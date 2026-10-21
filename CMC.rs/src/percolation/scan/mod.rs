//! Newman-Ziff parameter scans (Z1).
//!
//! Site and bond scans activate one uniformly random permutation
//! incrementally and record the occupation-count curve; binomial reweighting
//! turns one curve into scalar estimates at any uniform probability without a
//! second sampling pass. Heterogeneous probabilities, mixed scans, deletion
//! processes, and Carlo.rs result-vector output are outside this layer.

mod newman_ziff;
mod reweight;

#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use newman_ziff::WorkspaceScanCapacityAudit;
pub use newman_ziff::{
    scan_bond, scan_site, NewmanZiffWorkspace, ScanCurve, ScanError, ScanMode, ScanPoint,
};
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use reweight::ReweightCapacityAudit;
pub use reweight::{
    reweight, ReweightError, ReweightWorkspace, ReweightedEstimate, ScanEstimator, ScanObservable,
};
