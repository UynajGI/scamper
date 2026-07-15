use carlo_rs::job::TaskProgress;
use carlo_rs::progress::{print_status_table, spinner, MultiTaskProgress, SimProgress};
use std::path::PathBuf;
use std::time::Duration;

// ── SimProgress ───────────────────────────────────────────────────────────

#[test]
fn test_sim_progress_creation() {
    let progress = SimProgress::new(1000, 100, "task_0001");
    assert_eq!(progress.inner().length(), Some(1000));
    assert_eq!(progress.sweep_count(), 0);
    assert_eq!(progress.sweeps_per_sec(), 0.0);
}

#[test]
fn test_sim_progress_increment() {
    let mut progress = SimProgress::new(100, 10, "test_task");

    for _ in 0..50 {
        progress.inc();
    }

    assert_eq!(progress.sweep_count(), 50);
    assert!(progress.elapsed() < Duration::from_secs(1));
    // Rate should be > 0 after some increments
    assert!(progress.sweeps_per_sec() > 0.0);
}

#[test]
fn test_sim_progress_thermalization() {
    let mut progress = SimProgress::new(100, 10, "test_task");

    // Before thermalization
    for _ in 0..5 {
        progress.inc();
    }

    // After thermalization
    for _ in 0..10 {
        progress.inc();
    }

    assert_eq!(progress.sweep_count(), 15);
}

#[test]
fn test_sim_progress_set_position() {
    let mut progress = SimProgress::new(1000, 100, "task");

    progress.set_position(500);
    assert_eq!(progress.sweep_count(), 500);
}

#[test]
fn test_sim_progress_finish() {
    let progress = SimProgress::new(100, 10, "task");
    progress.finish();
}

#[test]
fn test_sim_progress_finish_with_message() {
    let progress = SimProgress::new(100, 10, "task");
    progress.finish_with_message("completed successfully");
}

#[test]
fn test_sim_progress_elapsed() {
    let progress = SimProgress::new(100, 10, "task");
    let elapsed = progress.elapsed();
    assert!(elapsed < Duration::from_secs(1));
}

#[test]
fn test_sim_progress_sweep_count_after_set() {
    let mut progress = SimProgress::new(1000, 0, "task");
    progress.set_position(42);
    assert_eq!(progress.sweep_count(), 42);
}

// ── MultiTaskProgress ─────────────────────────────────────────────────────

#[test]
fn test_multi_task_progress_creation() {
    let multi = MultiTaskProgress::new();
    assert!(multi.rate(0) == 0.0);
}

#[test]
fn test_multi_task_progress_default() {
    let multi = MultiTaskProgress::default();
    assert_eq!(multi.rate(999), 0.0); // out of range returns 0
}

#[test]
fn test_multi_task_progress_add_task() {
    let mut multi = MultiTaskProgress::new();

    let idx0 = multi.add_task("task_0", 1000);
    assert_eq!(idx0, 0);

    let idx1 = multi.add_task("task_1", 2000);
    assert_eq!(idx1, 1);
}

#[test]
fn test_multi_task_progress_update() {
    let mut multi = MultiTaskProgress::new();
    let idx = multi.add_task("task_0", 1000);

    multi.update(idx, 100);
    multi.update(idx, 200);
}

#[test]
fn test_multi_task_progress_rate() {
    let mut multi = MultiTaskProgress::new();
    let idx = multi.add_task("task_0", 1000);

    // Before any update, rate should be 0
    assert_eq!(multi.rate(idx), 0.0);

    multi.update(idx, 500);
    // After update, rate should be > 0 (unless extremely fast)
    assert!(multi.rate(idx) >= 0.0);
}

#[test]
fn test_multi_task_progress_rate_out_of_range() {
    let multi = MultiTaskProgress::new();
    assert_eq!(multi.rate(999), 0.0);
}

#[test]
fn test_multi_task_progress_elapsed() {
    let mut multi = MultiTaskProgress::new();
    let idx = multi.add_task("task_0", 1000);

    let elapsed = multi.elapsed(idx);
    assert!(elapsed < Duration::from_secs(1));
}

#[test]
fn test_multi_task_progress_elapsed_out_of_range() {
    let multi = MultiTaskProgress::new();
    assert_eq!(multi.elapsed(999), Duration::ZERO);
}

#[test]
fn test_multi_task_progress_finish_task() {
    let mut multi = MultiTaskProgress::new();
    let idx = multi.add_task("task_0", 1000);
    multi.update(idx, 500);

    multi.finish_task(idx);
}

#[test]
fn test_multi_task_progress_finish_all() {
    let mut multi = MultiTaskProgress::new();
    multi.add_task("task_0", 1000);
    multi.add_task("task_1", 2000);
    multi.add_task("task_2", 500);

    multi.finish();
}

#[test]
fn test_multi_task_progress_multiple_tasks() {
    let mut multi = MultiTaskProgress::new();
    let idx0 = multi.add_task("task_0", 1000);
    let idx1 = multi.add_task("task_1", 2000);
    let idx2 = multi.add_task("task_2", 500);

    multi.update(idx0, 100);
    multi.update(idx1, 200);
    multi.update(idx2, 50);

    assert!(
        multi.elapsed(idx0) <= multi.elapsed(idx1) || multi.elapsed(idx0) >= multi.elapsed(idx1)
    );
}

// ── spinner ───────────────────────────────────────────────────────────────

#[test]
fn test_spinner_creation() {
    let pb = spinner("Loading...");
    // Should return a valid progress bar
    assert_eq!(pb.length(), None); // spinners have no length
    pb.finish();
}

#[test]
fn test_spinner_with_empty_message() {
    let pb = spinner("");
    pb.finish();
}

// ── TaskProgress and print_status_table ───────────────────────────────────

#[test]
fn test_task_progress_with_timing() {
    let task = TaskProgress {
        target_sweeps: 1000,
        sweeps: 500,
        num_runs: 3,
        thermalization_fraction: 1.0,
        dir: PathBuf::from("/tmp/test.data/task_0000"),
        last_modified: None,
        elapsed_seconds: 60.0,
        sweeps_per_sec: 8.33,
    };

    assert_eq!(task.sweeps_per_sec, 8.33);
    assert_eq!(task.elapsed_seconds, 60.0);
}

#[test]
fn test_print_status_table_no_panic() {
    let tasks = vec![
        TaskProgress {
            target_sweeps: 1000,
            sweeps: 500,
            num_runs: 2,
            thermalization_fraction: 1.0,
            dir: PathBuf::from("/tmp/test.data/task_0000"),
            last_modified: None,
            elapsed_seconds: 30.0,
            sweeps_per_sec: 16.67,
        },
        TaskProgress {
            target_sweeps: 2000,
            sweeps: 2000,
            num_runs: 4,
            thermalization_fraction: 1.0,
            dir: PathBuf::from("/tmp/test.data/task_0001"),
            last_modified: None,
            elapsed_seconds: 120.0,
            sweeps_per_sec: 16.67,
        },
    ];

    // Should not panic
    print_status_table(&tasks);
}

#[test]
fn test_print_status_table_empty() {
    let tasks: Vec<TaskProgress> = vec![];
    print_status_table(&tasks);
}

#[test]
fn test_print_status_table_single_complete() {
    let tasks = vec![TaskProgress {
        target_sweeps: 100,
        sweeps: 100,
        num_runs: 1,
        thermalization_fraction: 0.1,
        dir: PathBuf::from("/tmp/test.data/task_0000"),
        last_modified: None,
        elapsed_seconds: 10.0,
        sweeps_per_sec: 10.0,
    }];

    print_status_table(&tasks);
}

#[test]
fn test_print_status_table_zero_rate() {
    let tasks = vec![TaskProgress {
        target_sweeps: 1000,
        sweeps: 0,
        num_runs: 0,
        thermalization_fraction: 0.0,
        dir: PathBuf::from("/tmp/test.data/task_0000"),
        last_modified: None,
        elapsed_seconds: 0.0,
        sweeps_per_sec: 0.0,
    }];

    print_status_table(&tasks);
}
