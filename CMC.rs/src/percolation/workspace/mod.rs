mod components;

#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use components::WorkspaceCapacityAudit;
pub use components::{ComponentLabel, ComponentWorkspace};
