use carlo_rs::{parse_duration, JobInfo, TaskMaker};
use std::time::Duration;

#[test]
fn test_full_workflow() {
    // Create tasks
    let mut tm = TaskMaker::new();
    tm.set("sweeps", "1000")
        .set("thermalization", "100")
        .set("binsize", "50")
        .task()
        .unwrap();

    let tasks = tm.make_tasks();

    // Create job
    let job = JobInfo::new(
        "/tmp/test_job",
        "TestMC",
        "Xoshiro256PlusPlus",
        tasks,
        Duration::from_secs(300),
        Duration::from_secs(3600),
        1,
    );

    assert_eq!(job.tasks.len(), 1);
}

#[test]
fn test_duration_parsing_integration() {
    let d1 = parse_duration("30").unwrap();
    assert_eq!(d1, Duration::from_secs(30));

    let d2 = parse_duration("1:30:45").unwrap();
    assert_eq!(d2, Duration::from_secs(5445));
}
