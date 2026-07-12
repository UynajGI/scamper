//! Weighted composition of heterogeneous particle update kinds.

use crate::particle::ParticleError;
use rand::{Rng, RngExt};

/// One weighted move entry.
#[derive(Debug, Clone, PartialEq)]
pub struct WeightedMove<K> {
    pub kind: K,
    pub weight: f64,
}

/// Reusable weighted selector whose weights can be frozen for production.
#[derive(Debug, Clone, PartialEq)]
pub struct MoveMixture<K> {
    entries: Vec<WeightedMove<K>>,
    total_weight: f64,
    frozen: bool,
}

impl<K> Default for MoveMixture<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K> MoveMixture<K> {
    pub const fn new() -> Self {
        Self {
            entries: Vec::new(),
            total_weight: 0.0,
            frozen: false,
        }
    }

    /// Append a move kind and return its stable entry index.
    pub fn add(&mut self, kind: K, weight: f64) -> Result<usize, ParticleError> {
        if self.frozen {
            return Err(ParticleError::InvalidMoveMixture(
                "production move mixture is frozen".to_string(),
            ));
        }
        validate_weight(weight)?;
        let new_total = self.total_weight + weight;
        if !new_total.is_finite() || new_total <= 0.0 {
            return Err(ParticleError::InvalidMoveMixture(
                "move mixture must have a finite positive total weight".to_string(),
            ));
        }
        let index = self.entries.len();
        self.entries.push(WeightedMove { kind, weight });
        self.total_weight = new_total;
        Ok(index)
    }

    /// Change one move weight during setup or thermalization.
    pub fn set_weight(&mut self, index: usize, weight: f64) -> Result<(), ParticleError> {
        if self.frozen {
            return Err(ParticleError::InvalidMoveMixture(
                "production move mixture is frozen".to_string(),
            ));
        }
        validate_weight(weight)?;
        let old_weight = self
            .entries
            .get(index)
            .ok_or_else(|| {
                ParticleError::InvalidMoveMixture("move-mixture index is out of range".to_string())
            })?
            .weight;
        let new_total = self.total_weight - old_weight + weight;
        if !new_total.is_finite() || new_total <= 0.0 {
            return Err(ParticleError::InvalidMoveMixture(
                "move mixture must have a finite positive total weight".to_string(),
            ));
        }
        self.entries[index].weight = weight;
        self.total_weight = new_total;
        Ok(())
    }

    /// Irreversibly freeze weights before production measurements.
    pub fn freeze(&mut self) {
        self.frozen = true;
    }

    #[inline]
    pub const fn is_frozen(&self) -> bool {
        self.frozen
    }

    #[inline]
    pub fn entries(&self) -> &[WeightedMove<K>] {
        &self.entries
    }

    #[inline]
    pub const fn total_weight(&self) -> f64 {
        self.total_weight
    }

    /// Draw one entry according to normalized positive weights.
    pub fn select(&self, rng: &mut impl Rng) -> &K {
        assert!(
            !self.entries.is_empty() && self.total_weight.is_finite() && self.total_weight > 0.0,
            "cannot sample an empty or zero-weight move mixture"
        );
        let threshold = rng.random_range(0.0..self.total_weight);
        let mut cumulative = 0.0;
        for entry in &self.entries {
            cumulative += entry.weight;
            if threshold < cumulative {
                return &entry.kind;
            }
        }
        &self.entries.last().expect("non-empty mixture").kind
    }
}

fn validate_weight(weight: f64) -> Result<(), ParticleError> {
    if weight.is_finite() && weight >= 0.0 {
        Ok(())
    } else {
        Err(ParticleError::InvalidMoveMixture(
            "move weight must be finite and non-negative".to_string(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn production_freeze_rejects_weight_changes() {
        let mut mixture = MoveMixture::new();
        mixture.add("translation", 1.0).unwrap();
        mixture.add("rotation", 1.0).unwrap();
        mixture.freeze();
        assert!(mixture.set_weight(0, 2.0).is_err());
        assert!(mixture.add("volume", 1.0).is_err());
    }

    #[test]
    fn zero_weight_entry_is_never_selected() {
        let mut mixture = MoveMixture::new();
        mixture.add(0usize, 0.0).unwrap_err();

        let mut mixture = MoveMixture::new();
        mixture.add(0usize, 1.0).unwrap();
        mixture.add(1usize, 0.0).unwrap();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(9);
        for _ in 0..100 {
            assert_eq!(*mixture.select(&mut rng), 0);
        }
    }
}
