//! Task builder for generating parameter sweeps.

use crate::job::{task_name, TaskInfo};
use crate::CarloError;
use std::collections::HashMap;

pub struct TaskMaker {
    tasks: Vec<TaskInfo>,
    current_params: HashMap<String, String>,
}

impl TaskMaker {
    pub fn new() -> Self {
        Self {
            tasks: Vec::new(),
            current_params: HashMap::new(),
        }
    }

    pub fn set(&mut self, key: &str, value: &str) -> &mut Self {
        self.current_params
            .insert(key.to_string(), value.to_string());
        self
    }

    pub fn task(&mut self) -> Result<&mut Self, CarloError> {
        let task_id = self.tasks.len() + 1;
        let name = task_name(task_id as u64);
        let task = TaskInfo::new(&name, self.current_params.clone())?;
        self.tasks.push(task);
        Ok(self)
    }

    pub fn make_tasks(self) -> Vec<TaskInfo> {
        self.tasks
    }
    pub fn current_task_name(&self) -> String {
        task_name((self.tasks.len() + 1) as u64)
    }
}

impl Default for TaskMaker {
    fn default() -> Self {
        Self::new()
    }
}
