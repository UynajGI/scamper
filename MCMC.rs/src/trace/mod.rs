mod memory;
mod view;

pub use memory::MemoryTrace;
pub use view::{ParameterIter, TraceView};

use crate::{EuclideanState, McmcError, TransitionReport};

/// Storage backend for retained posterior draws.
pub trait TraceStore: Send {
    fn record(
        &mut self,
        chain_id: usize,
        state: &EuclideanState,
        report: &TransitionReport,
    ) -> Result<bool, McmcError>;

    fn clear(&mut self);

    fn len(&self) -> usize;

    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
