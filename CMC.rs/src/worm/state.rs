//! Generic physical/worm-sector state container.

use super::WormError;

/// Sector of an extended-configuration-space simulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WormSector {
    /// A physical configuration with no explicit defects.
    Physical,
    /// An extended configuration carrying a head and a tail defect.
    Worm,
}

impl WormSector {
    pub(crate) const fn as_str(self) -> &'static str {
        match self {
            Self::Physical => "physical",
            Self::Worm => "worm",
        }
    }

    pub(crate) fn parse(value: &str) -> Result<Self, WormError> {
        match value {
            "physical" => Ok(Self::Physical),
            "worm" => Ok(Self::Worm),
            _ => Err(WormError::new(format!("unknown worm sector `{value}`"))),
        }
    }
}

/// Configuration plus optional ordered worm endpoints.
///
/// The model-specific layer validates whether the defects agree with the
/// configuration constraints. This type enforces the structural invariant:
/// physical states have no endpoints and worm states have both endpoints.
#[derive(Debug, Clone, PartialEq)]
pub struct WormState<C, D> {
    configuration: C,
    sector: WormSector,
    head: Option<D>,
    tail: Option<D>,
}

impl<C, D> WormState<C, D> {
    pub const fn new(configuration: C) -> Self {
        Self {
            configuration,
            sector: WormSector::Physical,
            head: None,
            tail: None,
        }
    }

    pub fn from_parts(
        configuration: C,
        sector: WormSector,
        head: Option<D>,
        tail: Option<D>,
    ) -> Result<Self, WormError> {
        let state = Self {
            configuration,
            sector,
            head,
            tail,
        };
        state.validate_structure()?;
        Ok(state)
    }

    #[inline]
    pub const fn configuration(&self) -> &C {
        &self.configuration
    }

    #[inline]
    pub fn configuration_mut(&mut self) -> &mut C {
        &mut self.configuration
    }

    #[inline]
    pub const fn sector(&self) -> WormSector {
        self.sector
    }

    #[inline]
    pub const fn is_physical(&self) -> bool {
        matches!(self.sector, WormSector::Physical)
    }

    #[inline]
    pub const fn is_worm(&self) -> bool {
        matches!(self.sector, WormSector::Worm)
    }

    #[inline]
    pub const fn head(&self) -> Option<&D> {
        self.head.as_ref()
    }

    #[inline]
    pub const fn tail(&self) -> Option<&D> {
        self.tail.as_ref()
    }

    pub fn open(&mut self, defect: D) -> Result<(), WormError>
    where
        D: Clone,
    {
        if !self.is_physical() {
            return Err(WormError::new(
                "cannot open a worm outside the physical sector",
            ));
        }
        self.sector = WormSector::Worm;
        self.head = Some(defect.clone());
        self.tail = Some(defect);
        Ok(())
    }

    pub fn move_head(&mut self, new_head: D) -> Result<(), WormError> {
        if !self.is_worm() {
            return Err(WormError::new(
                "cannot move a worm head in the physical sector",
            ));
        }
        self.head = Some(new_head);
        Ok(())
    }

    pub fn close(&mut self) -> Result<(), WormError>
    where
        D: PartialEq,
    {
        if !self.is_worm() {
            return Err(WormError::new("cannot close a worm in the physical sector"));
        }
        if self.head != self.tail {
            return Err(WormError::new(
                "a worm can close only when its head and tail coincide",
            ));
        }
        self.sector = WormSector::Physical;
        self.head = None;
        self.tail = None;
        Ok(())
    }

    pub fn validate_structure(&self) -> Result<(), WormError> {
        match (self.sector, self.head.is_some(), self.tail.is_some()) {
            (WormSector::Physical, false, false) | (WormSector::Worm, true, true) => Ok(()),
            (WormSector::Physical, _, _) => Err(WormError::new(
                "physical-sector state must not contain worm endpoints",
            )),
            (WormSector::Worm, _, _) => Err(WormError::new(
                "worm-sector state must contain both head and tail",
            )),
        }
    }
}
