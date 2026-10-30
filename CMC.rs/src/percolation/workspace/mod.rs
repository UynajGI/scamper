mod components;
mod winding;

#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use components::WorkspaceCapacityAudit;
pub use components::{ComponentLabel, ComponentWorkspace};
pub use winding::WindingWorkspace;
#[cfg(feature = "allocation-probe")]
#[doc(hidden)]
pub use winding::WorkspaceWindingCapacityAudit;
