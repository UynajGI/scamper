//! Tests for the job module.

use carlo_rs::{job::parse_duration, TaskInfo, TaskMaker};
use std::collections::HashMap;
use std::time::Duration;

#[test]
fn test_task_info_validation_success() {
    let mut params = HashMap::new();
    params.insert("sweeps".to_string(), "1000".to_string());
    params.insert("thermalization".to_string(), "100".to_string());
    params.insert("binsize".to_string(), "10".to_string());

    let result = TaskInfo::new("test_task", params);
    assert!(result.is_ok());

    let task = result.unwrap();
    assert_eq!(task.name(), "test_task");
}

#[test]
fn test_task_info_validation_missing_sweeps() {
    let mut params = HashMap::new();
    params.insert("thermalization".to_string(), "100".to_string());
    params.insert("binsize".to_string(), "10".to_string());

    let result = TaskInfo::new("test_task", params);
    assert!(result.is_err());
}

#[test]
fn test_task_info_validation_missing_thermalization() {
    let mut params = HashMap::new();
    params.insert("sweeps".to_string(), "1000".to_string());
    params.insert("binsize".to_string(), "10".to_string());

    let result = TaskInfo::new("test_task", params);
    assert!(result.is_err());
}

#[test]
fn test_task_info_validation_missing_binsize() {
    let mut params = HashMap::new();
    params.insert("sweeps".to_string(), "1000".to_string());
    params.insert("thermalization".to_string(), "100".to_string());

    let result = TaskInfo::new("test_task", params);
    assert!(result.is_err());
}

#[test]
fn test_task_info_get_parameter() {
    let mut params = HashMap::new();
    params.insert("sweeps".to_string(), "1000".to_string());
    params.insert("thermalization".to_string(), "100".to_string());
    params.insert("binsize".to_string(), "10".to_string());
    params.insert("beta".to_string(), "2.5".to_string());

    let task = TaskInfo::new("test_task", params).unwrap();

    let sweeps: Option<u64> = task.get("sweeps");
    assert_eq!(sweeps, Some(1000));

    let beta: Option<f64> = task.get("beta");
    assert_eq!(beta, Some(2.5));

    let missing: Option<u64> = task.get("nonexistent");
    assert_eq!(missing, None);
}

#[test]
fn test_parse_duration_seconds_only() {
    let duration = parse_duration("30").unwrap();
    assert_eq!(duration, Duration::from_secs(30));
}

#[test]
fn test_parse_duration_minutes_seconds() {
    let duration = parse_duration("5:30").unwrap();
    assert_eq!(duration, Duration::from_secs(5 * 60 + 30));
}

#[test]
fn test_parse_duration_hours_minutes_seconds() {
    let duration = parse_duration("1:30:45").unwrap();
    assert_eq!(duration, Duration::from_secs(3600 + 30 * 60 + 45));
}

#[test]
fn test_parse_duration_days_hours_minutes_seconds() {
    let duration = parse_duration("2-12:30:15").unwrap();
    assert_eq!(
        duration,
        Duration::from_secs(2 * 86400 + 12 * 3600 + 30 * 60 + 15)
    );
}

#[test]
fn test_parse_duration_invalid_format() {
    let result = parse_duration("invalid");
    assert!(result.is_err());
}

#[test]
fn test_task_maker_basic_usage() {
    let mut maker = TaskMaker::new();

    maker
        .set("sweeps", "1000")
        .set("thermalization", "100")
        .set("binsize", "10")
        .set("beta", "2.0")
        .task()
        .unwrap();

    let tasks = maker.make_tasks();
    assert_eq!(tasks.len(), 1);
    assert_eq!(tasks[0].name(), "task0001");

    let beta: Option<f64> = tasks[0].get("beta");
    assert_eq!(beta, Some(2.0));
}

#[test]
fn test_task_maker_multiple_tasks() {
    let mut maker = TaskMaker::new();

    maker
        .set("sweeps", "1000")
        .set("thermalization", "100")
        .set("binsize", "10");

    maker
        .set("beta", "1.0")
        .task()
        .unwrap()
        .set("beta", "2.0")
        .task()
        .unwrap()
        .set("beta", "3.0")
        .task()
        .unwrap();

    let tasks = maker.make_tasks();
    assert_eq!(tasks.len(), 3);

    assert_eq!(tasks[0].name(), "task0001");
    assert_eq!(tasks[1].name(), "task0002");
    assert_eq!(tasks[2].name(), "task0003");

    let beta0: Option<f64> = tasks[0].get("beta");
    let beta1: Option<f64> = tasks[1].get("beta");
    let beta2: Option<f64> = tasks[2].get("beta");

    assert_eq!(beta0, Some(1.0));
    assert_eq!(beta1, Some(2.0));
    assert_eq!(beta2, Some(3.0));
}

#[test]
fn test_task_maker_missing_required_params() {
    let mut maker = TaskMaker::new();

    // Only set one required param
    maker.set("sweeps", "1000");

    let result = maker.task();
    assert!(result.is_err());
}

#[test]
fn test_task_maker_default() {
    let maker = TaskMaker::default();
    let tasks = maker.make_tasks();
    assert_eq!(tasks.len(), 0);
}
