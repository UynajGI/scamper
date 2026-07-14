//! Straight event-chain Monte Carlo for identical periodic hard spheres.

use super::DynamicsError;
use crate::audit::should_audit_cache;
use crate::particle::{OrthorhombicCell, ParticleConfiguration, SimulationCell};
use rand::{Rng, RngExt};
use serde_json::{json, Value as Json};

const COLLISION_EPSILON: f64 = 1e-12;

/// One completed lifted event chain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EventChainOutcome {
    pub initial_particle: usize,
    pub final_particle: usize,
    pub axis: usize,
    pub direction: i8,
    pub distance: f64,
    pub collisions: u64,
}

/// Rejection-free straight event-chain kernel for identical hard spheres.
pub struct HardSphereEventChain<const D: usize> {
    configuration: ParticleConfiguration<D>,
    diameter: f64,
    chain_length: f64,
    chains: u64,
    collisions: u64,
    lifted_distance: f64,
    audit_interval: u64,
}

impl<const D: usize> HardSphereEventChain<D> {
    pub fn new(
        configuration: ParticleConfiguration<D>,
        diameter: f64,
        chain_length: f64,
        audit_interval: u64,
    ) -> Result<Self, DynamicsError> {
        let kernel = Self {
            configuration,
            diameter,
            chain_length,
            chains: 0,
            collisions: 0,
            lifted_distance: 0.0,
            audit_interval,
        };
        kernel.validate()?;
        Ok(kernel)
    }

    #[inline]
    pub const fn configuration(&self) -> &ParticleConfiguration<D> {
        &self.configuration
    }

    #[inline]
    pub const fn diameter(&self) -> f64 {
        self.diameter
    }

    #[inline]
    pub const fn chain_length(&self) -> f64 {
        self.chain_length
    }

    #[inline]
    pub const fn chains(&self) -> u64 {
        self.chains
    }

    #[inline]
    pub const fn collisions(&self) -> u64 {
        self.collisions
    }

    #[inline]
    pub const fn lifted_distance(&self) -> f64 {
        self.lifted_distance
    }

    pub fn validate(&self) -> Result<(), DynamicsError> {
        if D == 0 {
            return Err(DynamicsError::new(
                "hard-sphere event chain requires positive dimension",
            ));
        }
        if self.configuration.is_empty() {
            return Err(DynamicsError::new(
                "hard-sphere event chain requires at least one particle",
            ));
        }
        if !self.diameter.is_finite() || self.diameter <= 0.0 {
            return Err(DynamicsError::new(
                "hard-sphere diameter must be finite and positive",
            ));
        }
        if self.diameter > 0.5 * self.configuration.cell().minimum_length() {
            return Err(DynamicsError::new(
                "hard-sphere diameter must not exceed half the shortest cell length",
            ));
        }
        if !self.chain_length.is_finite() || self.chain_length <= 0.0 {
            return Err(DynamicsError::new(
                "event-chain length must be finite and positive",
            ));
        }
        if !self.lifted_distance.is_finite() || self.lifted_distance < 0.0 {
            return Err(DynamicsError::new("event-chain distance clock is invalid"));
        }
        let minimum_squared = self.diameter * self.diameter;
        for left in 0..self.configuration.len() {
            for right in left + 1..self.configuration.len() {
                let distance_squared = self.configuration.cell().distance_squared(
                    self.configuration.position(left),
                    self.configuration.position(right),
                );
                let tolerance = 1e-11 * (1.0 + minimum_squared);
                if distance_squared + tolerance < minimum_squared {
                    return Err(DynamicsError::new(format!(
                        "hard-sphere overlap between particles {left} and {right}"
                    )));
                }
            }
        }
        Ok(())
    }

    /// Execute one rejection-free chain with a random lifted particle and axis direction.
    pub fn step(&mut self, rng: &mut impl Rng) -> Result<EventChainOutcome, DynamicsError> {
        let initial_particle = rng.random_range(0..self.configuration.len());
        let axis = rng.random_range(0..D);
        let direction = if rng.random::<bool>() { 1 } else { -1 };
        self.step_with_lifting(initial_particle, axis, direction)
    }

    /// Deterministic lifted-chain entry point used by tests and custom schedulers.
    pub fn step_with_lifting(
        &mut self,
        initial_particle: usize,
        axis: usize,
        direction: i8,
    ) -> Result<EventChainOutcome, DynamicsError> {
        if initial_particle >= self.configuration.len() {
            return Err(DynamicsError::new(
                "event-chain active particle out of range",
            ));
        }
        if axis >= D {
            return Err(DynamicsError::new("event-chain axis out of range"));
        }
        if direction != -1 && direction != 1 {
            return Err(DynamicsError::new(
                "event-chain direction must be either -1 or +1",
            ));
        }

        let mut active = initial_particle;
        let mut remaining = self.chain_length;
        let mut collisions = 0u64;
        let mut guard = 0usize;
        let max_collisions = self
            .configuration
            .len()
            .saturating_mul(1_000_000)
            .max(1_000_000);
        while remaining > COLLISION_EPSILON * self.chain_length.max(1.0) {
            guard += 1;
            if guard > max_collisions {
                return Err(DynamicsError::new(
                    "event chain exceeded the collision safety bound",
                ));
            }
            let collision = self.next_collision(active, axis, direction, remaining)?;
            let (distance, next_active) = collision
                .map_or((remaining, None), |(distance, other)| {
                    (distance.min(remaining), Some(other))
                });
            self.translate(active, axis, direction, distance);
            remaining = (remaining - distance).max(0.0);
            if let Some(other) = next_active {
                active = other;
                collisions = collisions.saturating_add(1);
            }
            if distance == 0.0 && next_active.is_none() {
                break;
            }
        }
        if remaining > 0.0 {
            self.translate(active, axis, direction, remaining);
        }

        self.chains = self.chains.saturating_add(1);
        self.collisions = self.collisions.saturating_add(collisions);
        self.lifted_distance += self.chain_length;
        if should_audit_cache(self.chains, self.audit_interval) {
            self.validate()?;
        }
        Ok(EventChainOutcome {
            initial_particle,
            final_particle: active,
            axis,
            direction,
            distance: self.chain_length,
            collisions,
        })
    }

    fn translate(&mut self, particle: usize, axis: usize, direction: i8, distance: f64) {
        let mut position = *self.configuration.position(particle);
        position[axis] += f64::from(direction) * distance;
        self.configuration.cell().wrap(&mut position);
        self.configuration.set_position(particle, position);
    }

    fn next_collision(
        &self,
        active: usize,
        axis: usize,
        direction: i8,
        maximum_distance: f64,
    ) -> Result<Option<(f64, usize)>, DynamicsError> {
        let cell = self.configuration.cell();
        let active_position = self.configuration.position(active);
        let diameter_squared = self.diameter * self.diameter;
        let axis_length = cell.lengths()[axis];
        let sign = f64::from(direction);
        let mut best: Option<(f64, usize)> = None;

        for other in 0..self.configuration.len() {
            if other == active {
                continue;
            }
            let displacement =
                cell.displacement(active_position, self.configuration.position(other));
            let perpendicular_squared = displacement
                .iter()
                .enumerate()
                .filter_map(|(component_axis, component)| {
                    (component_axis != axis).then_some(component * component)
                })
                .sum::<f64>();
            if perpendicular_squared >= diameter_squared {
                continue;
            }
            let contact_projection = (diameter_squared - perpendicular_squared).sqrt();
            let oriented_parallel = sign * displacement[axis];
            let quotient =
                (COLLISION_EPSILON + contact_projection - oriented_parallel) / axis_length;
            let image_shift = quotient.floor() + 1.0;
            let distance = oriented_parallel + image_shift * axis_length - contact_projection;
            if !distance.is_finite() || distance <= COLLISION_EPSILON {
                continue;
            }
            if distance > maximum_distance + COLLISION_EPSILON {
                continue;
            }
            match best {
                None => best = Some((distance, other)),
                Some((best_distance, best_other)) => {
                    if distance < best_distance - COLLISION_EPSILON
                        || ((distance - best_distance).abs() <= COLLISION_EPSILON
                            && other < best_other)
                    {
                        best = Some((distance, other));
                    }
                }
            }
        }
        Ok(best)
    }

    pub fn save_snapshot(&self) -> Json {
        let cell_lengths = self.configuration.cell().lengths().to_vec();
        let positions = self
            .configuration
            .positions()
            .iter()
            .map(|position| position.to_vec())
            .collect::<Vec<_>>();
        json!({
            "format": "cmc-rs-hard-sphere-event-chain-v1",
            "dimension": D,
            "diameter": self.diameter,
            "chain_length": self.chain_length,
            "cell_lengths": cell_lengths,
            "positions": positions,
            "species": self.configuration.species(),
            "runtime": {
                "chains": self.chains,
                "collisions": self.collisions,
                "lifted_distance": self.lifted_distance,
                "audit_interval": self.audit_interval,
            }
        })
    }

    pub fn load_snapshot(&mut self, snapshot: &Json) -> Result<(), DynamicsError> {
        if snapshot["format"].as_str() != Some("cmc-rs-hard-sphere-event-chain-v1") {
            return Err(DynamicsError::new(
                "unknown hard-sphere event-chain snapshot format",
            ));
        }
        if snapshot["dimension"].as_u64() != Some(D as u64) {
            return Err(DynamicsError::new(
                "event-chain snapshot dimension mismatch",
            ));
        }
        require_same_f64(snapshot, "diameter", self.diameter)?;
        require_same_f64(snapshot, "chain_length", self.chain_length)?;
        let lengths = parse_array::<D>(&snapshot["cell_lengths"], "cell_lengths")?;
        if lengths
            .iter()
            .zip(self.configuration.cell().lengths())
            .any(|(left, right)| left.to_bits() != right.to_bits())
        {
            return Err(DynamicsError::new("event-chain snapshot cell mismatch"));
        }
        let positions = snapshot["positions"]
            .as_array()
            .ok_or_else(|| DynamicsError::new("event-chain snapshot positions must be an array"))?
            .iter()
            .map(|value| parse_array::<D>(value, "position"))
            .collect::<Result<Vec<_>, _>>()?;
        let species = snapshot["species"]
            .as_array()
            .ok_or_else(|| DynamicsError::new("event-chain snapshot species must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_u64()
                    .and_then(|number| u16::try_from(number).ok())
                    .ok_or_else(|| DynamicsError::new("event-chain snapshot species is invalid"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        if positions.len() != self.configuration.len() || species != self.configuration.species() {
            return Err(DynamicsError::new(
                "event-chain snapshot particle identity mismatch",
            ));
        }
        let cell = OrthorhombicCell::new(lengths)
            .map_err(|error| DynamicsError::new(error.to_string()))?;
        self.configuration = ParticleConfiguration::new(positions, species, cell)
            .map_err(|error| DynamicsError::new(error.to_string()))?;
        let runtime = &snapshot["runtime"];
        self.chains = required_u64(runtime, "chains")?;
        self.collisions = required_u64(runtime, "collisions")?;
        self.lifted_distance = required_f64(runtime, "lifted_distance")?;
        self.audit_interval = required_u64(runtime, "audit_interval")?;
        self.validate()
    }
}

fn parse_array<const D: usize>(value: &Json, field: &str) -> Result<[f64; D], DynamicsError> {
    let values = value
        .as_array()
        .ok_or_else(|| DynamicsError::new(format!("event-chain snapshot `{field}` is invalid")))?;
    if values.len() != D {
        return Err(DynamicsError::new(format!(
            "event-chain snapshot `{field}` dimension mismatch"
        )));
    }
    let mut result = [0.0; D];
    for (axis, value) in values.iter().enumerate() {
        result[axis] = value
            .as_f64()
            .filter(|number| number.is_finite())
            .ok_or_else(|| {
                DynamicsError::new(format!("event-chain snapshot `{field}` is invalid"))
            })?;
    }
    Ok(result)
}

fn required_u64(value: &Json, field: &str) -> Result<u64, DynamicsError> {
    value[field].as_u64().ok_or_else(|| {
        DynamicsError::new(format!("event-chain snapshot field `{field}` is invalid"))
    })
}

fn required_f64(value: &Json, field: &str) -> Result<f64, DynamicsError> {
    value[field]
        .as_f64()
        .filter(|number| number.is_finite())
        .ok_or_else(|| {
            DynamicsError::new(format!("event-chain snapshot field `{field}` is invalid"))
        })
}

fn require_same_f64(value: &Json, field: &str, expected: f64) -> Result<(), DynamicsError> {
    if value[field].as_f64().map(f64::to_bits) == Some(expected.to_bits()) {
        Ok(())
    } else {
        Err(DynamicsError::new(format!(
            "event-chain snapshot field `{field}` mismatch"
        )))
    }
}
