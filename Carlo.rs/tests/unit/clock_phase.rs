use carlo_rs::{RunPhase, SimulationClock};

// ── SimulationClock ───────────────────────────────────────────────────────

#[test]
fn test_clock_sweeps_value() {
    let clock = SimulationClock::Sweeps(42);
    assert_eq!(clock.value(), 42.0);
}

#[test]
fn test_clock_attempts_value() {
    let clock = SimulationClock::Attempts(99);
    assert_eq!(clock.value(), 99.0);
}

#[test]
fn test_clock_accepted_moves_value() {
    let clock = SimulationClock::AcceptedMoves(7);
    assert_eq!(clock.value(), 7.0);
}

#[test]
fn test_clock_event_time_value() {
    let clock = SimulationClock::EventTime(2.5);
    assert!((clock.value() - 2.5).abs() < 1e-10);
}

#[test]
fn test_clock_zero_values() {
    assert_eq!(SimulationClock::Sweeps(0).value(), 0.0);
    assert_eq!(SimulationClock::Attempts(0).value(), 0.0);
    assert_eq!(SimulationClock::AcceptedMoves(0).value(), 0.0);
    assert_eq!(SimulationClock::EventTime(0.0).value(), 0.0);
}

#[test]
fn test_clock_equality() {
    assert_eq!(SimulationClock::Sweeps(5), SimulationClock::Sweeps(5));
    assert_ne!(SimulationClock::Sweeps(5), SimulationClock::Attempts(5));
    assert_ne!(SimulationClock::EventTime(5.0), SimulationClock::Sweeps(5));
}

#[test]
fn test_clock_serde_roundtrip() {
    let clocks = [
        SimulationClock::Sweeps(10),
        SimulationClock::Attempts(20),
        SimulationClock::AcceptedMoves(30),
        SimulationClock::EventTime(1.5),
    ];
    for clock in clocks {
        let json = serde_json::to_string(&clock).unwrap();
        let restored: SimulationClock = serde_json::from_str(&json).unwrap();
        assert_eq!(clock, restored);
    }
}

#[test]
fn test_clock_copy_clone() {
    let original = SimulationClock::Sweeps(5);
    let copied = original;
    let cloned = original; // Copy, not clone
    assert_eq!(original, copied);
    assert_eq!(original, cloned);
}

#[test]
fn test_clock_debug_format() {
    let clock = SimulationClock::EventTime(2.5);
    let debug = format!("{:?}", clock);
    assert!(debug.contains("EventTime"));
    assert!(debug.contains("2.5"));
}

// ── RunPhase ──────────────────────────────────────────────────────────────

#[test]
fn test_run_phase_default_is_initialization() {
    assert_eq!(RunPhase::default(), RunPhase::Initialization);
}

#[test]
fn test_phase_allows_adaptation() {
    assert!(!RunPhase::Initialization.allows_adaptation());
    assert!(RunPhase::Thermalization.allows_adaptation());
    assert!(!RunPhase::Measurement.allows_adaptation());
    assert!(!RunPhase::Finished.allows_adaptation());
}

#[test]
fn test_phase_collects_measurements() {
    assert!(!RunPhase::Initialization.collects_measurements());
    assert!(!RunPhase::Thermalization.collects_measurements());
    assert!(RunPhase::Measurement.collects_measurements());
    assert!(!RunPhase::Finished.collects_measurements());
}

#[test]
fn test_phase_equality_and_ordering() {
    let phases = [
        RunPhase::Initialization,
        RunPhase::Thermalization,
        RunPhase::Measurement,
        RunPhase::Finished,
    ];
    for (i, p) in phases.iter().enumerate() {
        for (j, q) in phases.iter().enumerate() {
            if i == j {
                assert_eq!(p, q);
            } else {
                assert_ne!(p, q);
            }
        }
    }
}

#[test]
fn test_phase_serde_roundtrip() {
    let phases = [
        RunPhase::Initialization,
        RunPhase::Thermalization,
        RunPhase::Measurement,
        RunPhase::Finished,
    ];
    for phase in phases {
        let json = serde_json::to_string(&phase).unwrap();
        let restored: RunPhase = serde_json::from_str(&json).unwrap();
        assert_eq!(phase, restored);
    }
}

#[test]
fn test_phase_copy_clone() {
    let original = RunPhase::Measurement;
    let copied = original;
    let cloned = original; // Copy, not clone
    assert_eq!(original, copied);
    assert_eq!(original, cloned);
}

#[test]
fn test_phase_debug_format() {
    let debug = format!("{:?}", RunPhase::Finished);
    assert!(debug.contains("Finished"));
}
