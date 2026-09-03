//! Zero-cost read-only topology capabilities.
//!
//! The same generic code accepts owned and borrowed topology without conversion:
//!
//! ```
//! use cmc_rs::{
//!     build_chain, BorrowedUndirectedCsr, GraphView, UndirectedGraphView, VertexId,
//! };
//!
//! fn degree<G: UndirectedGraphView>(graph: &G, vertex: VertexId) -> usize {
//!     graph.incidences(vertex).count()
//! }
//!
//! let owned = build_chain(3, false);
//! assert_eq!(degree(&owned, owned.vertex_id(1).unwrap()), 2);
//!
//! let offsets = [0, 1, 3, 4];
//! let neighbors = [1, 0, 2, 1];
//! let edge_ids = [0, 0, 1, 1];
//! let endpoints = [[0, 1], [1, 2]];
//! let borrowed = BorrowedUndirectedCsr::new(
//!     &offsets, &neighbors, &edge_ids, &endpoints,
//! ).unwrap();
//! assert_eq!(degree(&borrowed, borrowed.vertex_id(1).unwrap()), 2);
//! ```

mod borrowed_csr;
mod error;
mod id;
mod view;

pub use borrowed_csr::BorrowedUndirectedCsr;
pub use error::TopologyError;
pub use id::{EdgeId, VertexId};
pub use view::{GraphView, Incidence, UndirectedGraphView};
