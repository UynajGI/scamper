//! Progress reporting for Monte Carlo simulations.
//!
//! Provides progress bars and status display using indicatif.

use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::time::{Duration, Instant};

/// Progress bar for a single simulation run.
pub struct SimProgress {
    pb: ProgressBar,
    thermalization: u64,
    thermalized: bool,
    start_time: Instant,
    last_rate_update: Instant,
    last_sweep_count: u64,
    /// Current sweep count.
    sweep_count: u64,
    /// Average sweeps per second.
    sweeps_per_sec: f64,
}

impl SimProgress {
    /// Create a new progress bar for a simulation.
    pub fn new(target_sweeps: u64, thermalization: u64, task_name: &str) -> Self {
        let pb = ProgressBar::new(target_sweeps);
        pb.set_style(
            ProgressStyle::with_template(
                "{msg:30} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({sweeps}/s)",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        pb.set_message(format!("{} [thermalizing] [0.0 s/s]", task_name));

        let now = Instant::now();
        Self {
            pb,
            thermalization,
            thermalized: false,
            start_time: now,
            last_rate_update: now,
            last_sweep_count: 0,
            sweep_count: 0,
            sweeps_per_sec: 0.0,
        }
    }

    /// Advance progress by one sweep.
    pub fn inc(&mut self) {
        self.pb.inc(1);
        self.sweep_count += 1;

        // Update rate every ~100ms
        let now = Instant::now();
        if now.duration_since(self.last_rate_update) >= Duration::from_millis(100) {
            let elapsed = now.duration_since(self.last_rate_update);
            let delta = self.sweep_count - self.last_sweep_count;
            if elapsed.as_secs_f64() > 0.0 {
                self.sweeps_per_sec = delta as f64 / elapsed.as_secs_f64();
            }
            self.last_rate_update = now;
            self.last_sweep_count = self.sweep_count;
        }

        // Update message with rate
        let task_name: String = self.pb.message().replace(" [thermalizing]", "");
        if !self.thermalized && self.pb.position() > self.thermalization {
            self.thermalized = true;
            self.pb.set_message(format!(
                "{} [measuring] [{:.1} s/s]",
                task_name, self.sweeps_per_sec
            ));
        } else {
            self.pb.set_message(format!(
                "{} [{:.1} s/s]",
                task_name.strip_suffix("]").unwrap_or(&task_name),
                self.sweeps_per_sec
            ));
        }
    }

    /// Update with current sweep count.
    pub fn set_position(&mut self, pos: u64) {
        self.pb.set_position(pos);
        self.sweep_count = pos;
    }

    /// Mark as finished.
    pub fn finish(&self) {
        self.pb.finish_with_message("done");
    }

    /// Finish with a custom message.
    pub fn finish_with_message(&self, msg: &str) {
        self.pb.finish_with_message(msg.to_string());
    }

    /// Get the inner progress bar for advanced customization.
    pub fn inner(&self) -> &ProgressBar {
        &self.pb
    }

    /// Get elapsed time.
    pub fn elapsed(&self) -> Duration {
        self.start_time.elapsed()
    }

    /// Get average sweeps per second.
    pub fn sweeps_per_sec(&self) -> f64 {
        if self.start_time.elapsed().as_secs_f64() > 0.0 {
            self.sweep_count as f64 / self.start_time.elapsed().as_secs_f64()
        } else {
            0.0
        }
    }

    /// Get current sweep count.
    pub fn sweep_count(&self) -> u64 {
        self.sweep_count
    }
}

/// Multi-task progress display.
pub struct MultiTaskProgress {
    multi: MultiProgress,
    bars: Vec<ProgressBar>,
    start_times: Vec<Instant>,
    last_sweeps: Vec<u64>,
    last_updates: Vec<Instant>,
    rates: Vec<f64>,
}

impl MultiTaskProgress {
    /// Create a new multi-task progress display.
    pub fn new() -> Self {
        Self {
            multi: MultiProgress::new(),
            bars: Vec::new(),
            start_times: Vec::new(),
            last_sweeps: Vec::new(),
            last_updates: Vec::new(),
            rates: Vec::new(),
        }
    }

    /// Add a task with its target sweeps.
    pub fn add_task(&mut self, task_name: &str, target_sweeps: u64) -> usize {
        let pb = self.multi.add(ProgressBar::new(target_sweeps));
        pb.set_style(
            ProgressStyle::with_template(
                "{msg:20} [{elapsed_precise}] {bar:30.cyan/blue} {pos}/{len} ({rate:.1}/s)",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        pb.set_message(task_name.to_string());

        let now = Instant::now();
        self.start_times.push(now);
        self.last_sweeps.push(0);
        self.last_updates.push(now);
        self.rates.push(0.0);

        self.bars.push(pb);
        self.bars.len() - 1
    }

    /// Update progress for a specific task.
    pub fn update(&mut self, task_idx: usize, sweeps: u64) {
        if let Some(pb) = self.bars.get(task_idx) {
            let now = Instant::now();
            let elapsed = now.duration_since(self.last_updates[task_idx]);
            let delta = sweeps.saturating_sub(self.last_sweeps[task_idx]);

            if elapsed.as_secs_f64() > 0.0 {
                self.rates[task_idx] = delta as f64 / elapsed.as_secs_f64();
            }

            self.last_sweeps[task_idx] = sweeps;
            self.last_updates[task_idx] = now;

            pb.set_position(sweeps);
        }
    }

    /// Get average rate for a task (sweeps/second).
    pub fn rate(&self, task_idx: usize) -> f64 {
        self.rates.get(task_idx).copied().unwrap_or(0.0)
    }

    /// Get elapsed time for a task.
    pub fn elapsed(&self, task_idx: usize) -> Duration {
        if let Some(&start) = self.start_times.get(task_idx) {
            start.elapsed()
        } else {
            Duration::ZERO
        }
    }

    /// Mark a task as complete.
    pub fn finish_task(&self, task_idx: usize) {
        if let Some(pb) = self.bars.get(task_idx) {
            pb.finish_with_message("done");
        }
    }

    /// Finish all tasks.
    pub fn finish(&self) {
        for pb in &self.bars {
            pb.finish();
        }
    }
}

impl Default for MultiTaskProgress {
    fn default() -> Self {
        Self::new()
    }
}

/// Display a status table for tasks.
pub fn print_status_table(tasks: &[super::job::TaskProgress]) {
    use std::io::Write;
    use termcolor::{Color, ColorChoice, ColorSpec, StandardStream, WriteColor};

    let stdout = StandardStream::stdout(ColorChoice::Auto);
    let mut handle = stdout.lock();

    // Calculate column widths
    let name_width = tasks
        .iter()
        .map(|t| {
            t.dir
                .file_name()
                .map(|n| n.to_string_lossy().len())
                .unwrap_or(4)
        })
        .max()
        .unwrap_or(10)
        .max(4);
    let sweeps_width = tasks
        .iter()
        .map(|t| t.sweeps.to_string().len())
        .max()
        .unwrap_or(1)
        .max(7);
    let target_width = tasks
        .iter()
        .map(|t| t.target_sweeps.to_string().len())
        .max()
        .unwrap_or(1)
        .max(6);

    // Format elapsed time for display
    fn format_duration(secs: f64) -> String {
        if secs < 60.0 {
            format!("{:.0}s", secs)
        } else if secs < 3600.0 {
            format!("{:.0}m", secs / 60.0)
        } else {
            format!("{:.1}h", secs / 3600.0)
        }
    }

    // Header
    let _ = writeln!(handle);
    let mut bold_spec = ColorSpec::new();
    bold_spec.set_bold(true);
    let _ = handle.set_color(&bold_spec);
    let _ = writeln!(
        handle,
        " {:name_width$}  {:>sweeps_width$}  {:>target_width$}  {:>4}  {:>5}  {:>8}  {:>6}",
        "Task",
        "Sweeps",
        "Target",
        "Runs",
        "Therm",
        "Elapsed",
        "Rate",
        name_width = name_width,
        sweeps_width = sweeps_width,
        target_width = target_width
    );
    let _ = handle.reset();
    let _ = writeln!(
        handle,
        " {}",
        "-".repeat(name_width + sweeps_width + target_width + 35)
    );

    // Rows
    for task in tasks {
        let name = task
            .dir
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "????".to_string());
        let therm_pct = (task.thermalization_fraction * 100.0) as u32;
        let elapsed = format_duration(task.elapsed_seconds);
        let rate = if task.sweeps_per_sec > 0.0 {
            format!("{:.1}/s", task.sweeps_per_sec)
        } else {
            "-".to_string()
        };

        let _ = writeln!(
            handle,
            " {:name_width$}  {:>sweeps_width$}  {:>target_width$}  {:>4}  {:>3}%  {:>8}  {:>6}",
            name,
            task.sweeps,
            task.target_sweeps,
            task.num_runs,
            therm_pct,
            elapsed,
            rate,
            name_width = name_width,
            sweeps_width = sweeps_width,
            target_width = target_width
        );
    }
    let _ = writeln!(handle);

    // Summary
    let all_done = tasks.iter().all(|t| t.sweeps >= t.target_sweeps);
    if all_done {
        let mut green_spec = ColorSpec::new();
        green_spec.set_fg(Some(Color::Green));
        let _ = handle.set_color(&green_spec);
        let _ = writeln!(handle, " ✓ All tasks complete");
        let _ = handle.reset();
    }

    // Total stats
    let total_sweeps: u64 = tasks.iter().map(|t| t.sweeps).sum();
    let total_runs: u64 = tasks.iter().map(|t| t.num_runs).sum();
    let avg_rate: f64 = tasks
        .iter()
        .map(|t| t.sweeps_per_sec)
        .sum::<f64>()
        / tasks.len().max(1) as f64;
    let _ = writeln!(
        handle,
        " Total: {} sweeps, {} runs, avg rate {:.1} sweeps/s\n",
        total_sweeps, total_runs, avg_rate
    );
}

/// Spinner for operations with unknown duration.
pub fn spinner(msg: &str) -> ProgressBar {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} {msg}")
            .unwrap()
            .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]),
    );
    pb.set_message(msg.to_string());
    pb.enable_steady_tick(Duration::from_millis(100));
    pb
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sim_progress_creation() {
        let progress = SimProgress::new(1000, 100, "task_0001");
        assert_eq!(progress.pb.length(), Some(1000));
    }
}
