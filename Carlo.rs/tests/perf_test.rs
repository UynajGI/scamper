use carlo_rs::job::TaskProgress;
use carlo_rs::progress::print_status_table;
use carlo_rs::progress::SimProgress;
use std::path::PathBuf;
use std::time::Duration;

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
