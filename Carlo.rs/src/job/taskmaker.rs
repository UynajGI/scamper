//! Task builder for generating parameter sweeps.

use crate::job::{task_name, TaskInfo};
use crate::CarloError;
use std::collections::HashMap;

/// Builder for incrementally constructing a list of [`TaskInfo`] values.
///
/// Call [`set`](TaskMaker::set) to populate shared parameters, then
/// [`task`](TaskMaker::task) to snapshot them into a new task.
pub struct TaskMaker {
    tasks: Vec<TaskInfo>,
    current_params: HashMap<String, String>,
}

impl TaskMaker {
    /// Create an empty builder.
    pub fn new() -> Self {
        Self {
            tasks: Vec::new(),
            current_params: HashMap::new(),
        }
    }

    /// Set a parameter for subsequent tasks. Returns `&mut self` for chaining.
    pub fn set(&mut self, key: &str, value: &str) -> &mut Self {
        self.current_params
            .insert(key.to_string(), value.to_string());
        self
    }

    /// Snapshot the current parameters into a new task. Returns an error if
    /// required keys (`sweeps`, `thermalization`, `binsize`) are missing.
    pub fn task(&mut self) -> Result<&mut Self, CarloError> {
        let task_id = self.tasks.len() + 1;
        let name = task_name(task_id as u64);
        let task = TaskInfo::new(&name, self.current_params.clone())?;
        self.tasks.push(task);
        Ok(self)
    }

    /// Consume the builder and return the completed task list.
    pub fn make_tasks(self) -> Vec<TaskInfo> {
        self.tasks
    }

    /// Name the next task would receive.
    pub fn current_task_name(&self) -> String {
        task_name((self.tasks.len() + 1) as u64)
    }
}

impl Default for TaskMaker {
    fn default() -> Self {
        Self::new()
    }
}
