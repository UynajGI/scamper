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
    // Rate should be finite and > 0 after some increments
    let rate = progress.sweeps_per_sec();
    assert!(
        rate.is_finite(),
        "sweeps_per_sec should be finite, got {rate}"
    );
    assert!(
        rate > 0.0,
        "sweeps_per_sec should be > 0 after increments, got {rate}"
    );
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
    assert!(
        progress.inner().is_finished(),
        "progress bar should be finished after finish()"
    );
}

#[test]
fn test_sim_progress_finish_with_message() {
    let progress = SimProgress::new(100, 10, "task");
    progress.finish_with_message("completed successfully");
    assert!(
        progress.inner().is_finished(),
        "progress bar should be finished after finish_with_message()"
    );
}

#[test]
fn test_sim_progress_elapsed() {
    let progress = SimProgress::new(100, 10, "task");
    let elapsed = progress.elapsed();
    // Elapsed must be non-negative (Duration is unsigned) and reasonably small
    assert!(elapsed >= Duration::ZERO);
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

    // After updates, rate should be finite and non-negative
    let rate = multi.rate(idx);
    assert!(
        rate.is_finite(),
        "rate should be finite after updates, got {rate}"
    );
    assert!(rate >= 0.0, "rate should be non-negative, got {rate}");
}

#[test]
fn test_multi_task_progress_rate() {
    let mut multi = MultiTaskProgress::new();
    let idx = multi.add_task("task_0", 1000);

    // Before any update, rate should be 0
    assert_eq!(multi.rate(idx), 0.0);

    multi.update(idx, 500);
    // After update, rate should be finite and non-negative
    let rate = multi.rate(idx);
    assert!(rate.is_finite(), "rate should be finite, got {rate}");
    assert!(rate >= 0.0, "rate should be non-negative, got {rate}");
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
    assert!(elapsed >= Duration::ZERO, "elapsed should be non-negative");
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

    // finish_task marks the internal progress bar as done.
    // The bars field is private, so we verify indirectly: the call must not
    // panic and the task's elapsed time remains accessible afterward.
    multi.finish_task(idx);
    assert!(multi.elapsed(idx) >= Duration::ZERO);
}

#[test]
fn test_multi_task_progress_finish_all() {
    let mut multi = MultiTaskProgress::new();
    multi.add_task("task_0", 1000);
    multi.add_task("task_1", 2000);
    multi.add_task("task_2", 500);

    // finish() marks all internal progress bars as done.
    // The bars field is private, so we verify indirectly: the call must not
    // panic and all tasks' elapsed times remain accessible afterward.
    multi.finish();
    for idx in 0..3 {
        assert!(multi.elapsed(idx) >= Duration::ZERO);
    }
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

    // All three tasks should exist with non-negative elapsed times
    for &idx in &[idx0, idx1, idx2] {
        assert!(multi.elapsed(idx) >= Duration::ZERO);
        assert!(multi.rate(idx).is_finite());
    }
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
    // No-panic smoke test: empty message edge case for spinner
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
    // Intentional no-panic smoke test: print_status_table writes to stdout
    // via indicatif and cannot be captured in-process. Verifying no panic
    // across varied inputs (multiple tasks, different completion states) is
    // the meaningful assertion here.
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

    // Should not panic across varied inputs
    print_status_table(&tasks);
}

#[test]
fn test_print_status_table_empty() {
    // No-panic smoke test: empty task list edge case
    let tasks: Vec<TaskProgress> = vec![];
    print_status_table(&tasks);
}

#[test]
fn test_print_status_table_single_complete() {
    // No-panic smoke test: single completed task
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
    // No-panic smoke test: zero sweeps / zero rate edge case
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
