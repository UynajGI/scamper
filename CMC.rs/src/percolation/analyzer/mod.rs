mod components;
mod winding;

pub use components::{analyze, analyze_with_labels, AnalysisError, AnalysisResult};
pub use winding::{analyze_winding, ComponentWinding, WindingError, WindingResult, WindingSummary};
