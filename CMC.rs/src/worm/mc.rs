//! Carlo.rs adapter and versioned JSON snapshot for the Ising graph worm.

use super::{
    EndpointPairHistogram, IsingGraphConfiguration, IsingGraphWormModel, WormConfig, WormError,
    WormKernel, WormSector, WormState, WormTransitionStatistics,
};
use crate::classical_mc::{build_lattice_from_params, parse_bool, parse_param};
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, Params};
use serde_json::{json, Value as Json};

/// Scheduler-ready persistent high-temperature Ising graph worm.
pub struct IsingGraphWormMC {
    kernel: WormKernel<IsingGraphWormModel>,
    endpoint_sample: Vec<f64>,
}

impl IsingGraphWormMC {
    pub fn new(model: IsingGraphWormModel, config: WormConfig) -> Result<Self, WormError> {
        let configuration = model.empty_configuration();
        let endpoint_sample = if config.track_endpoint_pairs {
            vec![
                0.0;
                model
                    .lattice()
                    .n_sites
                    .saturating_mul(model.lattice().n_sites)
            ]
        } else {
            Vec::new()
        };
        Ok(Self {
            kernel: WormKernel::new(model, configuration, config)?,
            endpoint_sample,
        })
    }

    #[inline]
    pub const fn kernel(&self) -> &WormKernel<IsingGraphWormModel> {
        &self.kernel
    }

    #[inline]
    pub fn kernel_mut(&mut self) -> &mut WormKernel<IsingGraphWormModel> {
        &mut self.kernel
    }

    /// Versioned model/kernel snapshot. RNG and measurement state remain owned
    /// by Carlo.rs's context checkpoint.
    pub fn save_snapshot(&self) -> Json {
        let state = self.kernel.state();
        let statistics = self.kernel.statistics();
        let endpoint_pairs = self.kernel.endpoint_pairs().map(|histogram| {
            json!({
                "bins": histogram.bins(),
                "counts": histogram.counts(),
                "samples": histogram.samples(),
            })
        });
        json!({
            "format": "cmc-rs-ising-worm-v1",
            "model": {
                "beta": self.kernel.model().beta(),
                "coupling": self.kernel.model().coupling(),
                "n_sites": self.kernel.model().lattice().n_sites,
                "n_edges": self.kernel.model().lattice().n_edges(),
                "edges": self.kernel.model().lattice().edges.iter().map(|edge| json!({
                    "source": edge.source,
                    "target": edge.target,
                    "kind": edge.kind.as_label(),
                    "weight": edge.weight,
                })).collect::<Vec<_>>(),
            },
            "config": {
                "local_updates_per_sweep": self.kernel.config().local_updates_per_sweep,
                "close_probability": self.kernel.config().close_probability,
                "log_worm_fugacity": self.kernel.config().log_worm_fugacity,
                "track_endpoint_pairs": self.kernel.config().track_endpoint_pairs,
                "cache_audit_interval": self.kernel.config().cache_audit_interval,
            },
            "state": {
                "sector": state.sector().as_str(),
                "head": state.head(),
                "tail": state.tail(),
                "occupied": state.configuration().occupied(),
            },
            "runtime": {
                "sweeps": self.kernel.sweeps(),
                "current_worm_steps": self.kernel.current_worm_steps(),
                "statistics": statistics_to_json(*statistics),
                "endpoint_pairs": endpoint_pairs,
            },
        })
    }

    pub fn load_snapshot(&mut self, snapshot: &Json) -> Result<(), CarloError> {
        if snapshot["format"].as_str() != Some("cmc-rs-ising-worm-v1") {
            return Err(checkpoint_error("unknown Ising worm snapshot format"));
        }
        self.validate_snapshot_model(&snapshot["model"])?;
        self.validate_snapshot_config(&snapshot["config"])?;

        let occupied = snapshot["state"]["occupied"]
            .as_array()
            .ok_or_else(|| checkpoint_error("state.occupied must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_bool()
                    .ok_or_else(|| checkpoint_error("state.occupied contains a non-boolean"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let configuration = IsingGraphConfiguration::from_occupied(self.kernel.model(), occupied)
            .map_err(worm_checkpoint_error)?;
        let sector = WormSector::parse(
            snapshot["state"]["sector"]
                .as_str()
                .ok_or_else(|| checkpoint_error("state.sector is missing"))?,
        )
        .map_err(worm_checkpoint_error)?;
        let head = optional_usize(&snapshot["state"]["head"], "state.head")?;
        let tail = optional_usize(&snapshot["state"]["tail"], "state.tail")?;
        let state = WormState::from_parts(configuration, sector, head, tail)
            .map_err(worm_checkpoint_error)?;

        let runtime = &snapshot["runtime"];
        let sweeps = required_u64(runtime, "sweeps")?;
        let current_worm_steps = required_u64(runtime, "current_worm_steps")?;
        let statistics = statistics_from_json(&runtime["statistics"])?;
        let endpoint_pairs = if runtime["endpoint_pairs"].is_null() {
            None
        } else {
            let value = &runtime["endpoint_pairs"];
            let bins = required_usize(value, "bins")?;
            let samples = required_u64(value, "samples")?;
            let counts = value["counts"]
                .as_array()
                .ok_or_else(|| checkpoint_error("endpoint_pairs.counts must be an array"))?
                .iter()
                .map(|entry| {
                    entry
                        .as_u64()
                        .ok_or_else(|| checkpoint_error("endpoint_pairs.counts contains a non-u64"))
                })
                .collect::<Result<Vec<_>, _>>()?;
            Some(
                EndpointPairHistogram::from_counts(bins, counts, samples)
                    .map_err(worm_checkpoint_error)?,
            )
        };
        self.kernel
            .restore_runtime(
                state,
                statistics,
                endpoint_pairs,
                sweeps,
                current_worm_steps,
            )
            .map_err(worm_checkpoint_error)
    }

    fn validate_snapshot_model(&self, model: &Json) -> Result<(), CarloError> {
        let expected = self.kernel.model();
        require_same_f64(model, "beta", expected.beta())?;
        require_same_f64(model, "coupling", expected.coupling())?;
        if required_usize(model, "n_sites")? != expected.lattice().n_sites
            || required_usize(model, "n_edges")? != expected.lattice().n_edges()
        {
            return Err(checkpoint_error(
                "Ising worm snapshot topology size mismatch",
            ));
        }
        let edges = model["edges"]
            .as_array()
            .ok_or_else(|| checkpoint_error("model.edges must be an array"))?;
        if edges.len() != expected.lattice().edges.len() {
            return Err(checkpoint_error("Ising worm snapshot edge-count mismatch"));
        }
        for (value, edge) in edges.iter().zip(&expected.lattice().edges) {
            if required_usize(value, "source")? != edge.source
                || required_usize(value, "target")? != edge.target
                || value["kind"].as_str() != Some(edge.kind.as_label())
                || value["weight"].as_f64().map(f64::to_bits) != Some(edge.weight.to_bits())
            {
                return Err(checkpoint_error("Ising worm snapshot edge mismatch"));
            }
        }
        Ok(())
    }

    fn validate_snapshot_config(&self, config: &Json) -> Result<(), CarloError> {
        let expected = self.kernel.config();
        if required_usize(config, "local_updates_per_sweep")? != expected.local_updates_per_sweep
            || config["close_probability"].as_f64().map(f64::to_bits)
                != Some(expected.close_probability.to_bits())
            || config["log_worm_fugacity"].as_f64().map(f64::to_bits)
                != Some(expected.log_worm_fugacity.to_bits())
            || config["track_endpoint_pairs"].as_bool() != Some(expected.track_endpoint_pairs)
            || required_u64(config, "cache_audit_interval")? != expected.cache_audit_interval
        {
            return Err(checkpoint_error(
                "Ising worm snapshot transition configuration mismatch",
            ));
        }
        Ok(())
    }
}

impl MonteCarlo for IsingGraphWormMC {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, context: &mut Context<Self::Rng>) {
        self.kernel
            .sweep(&mut context.rng)
            .expect("Ising graph worm sweep failed");
    }

    fn measure(&mut self, context: &mut Context<Self::Rng>) {
        let state = self.kernel.state();
        let last = *self.kernel.last_sweep_statistics();
        context.measure("WormSector", if state.is_worm() { 1.0 } else { 0.0 });
        context.measure(
            "PhysicalSector",
            if state.is_physical() { 1.0 } else { 0.0 },
        );
        context.measure("WormStepAcceptance", last.step_acceptance_fraction());
        context.measure("WormOpenAcceptance", last.open_acceptance_fraction());
        context.measure("WormCloseAcceptance", last.close_acceptance_fraction());
        context.measure(
            "CompletedWormLength",
            self.kernel.statistics().last_completed_worm_steps as f64,
        );

        if state.is_physical() {
            let occupied = state.configuration().occupied_edges() as f64;
            let n_edges = self.kernel.model().lattice().n_edges().max(1) as f64;
            context.measure("GraphOccupiedEdges", occupied);
            context.measure("GraphEdgeDensity", occupied / n_edges);
            context.measure(
                "Energy",
                self.kernel.model().energy_estimator(state.configuration()),
            );
        } else {
            let head = *state.head().expect("worm sector has a head");
            let tail = *state.tail().expect("worm sector has a tail");
            context.measure("WormHead", head as f64);
            context.measure("WormTail", tail as f64);
            if self.kernel.config().track_endpoint_pairs {
                self.endpoint_sample.fill(0.0);
                let n_sites = self.kernel.model().lattice().n_sites;
                self.endpoint_sample[tail * n_sites + head] = 1.0;
                context.measure_array("WormEndpointPairs", &self.endpoint_sample);
            }
        }
    }

    fn name(&self) -> &'static str {
        "IsingHighTemperatureGraphWorm"
    }
}

impl FromParams for IsingGraphWormMC {
    fn validate_params(params: &Params) -> Result<(), CarloError> {
        build_from_params(params).map(|_| ())
    }

    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let (model, config) = build_from_params(params)?;
        Self::new(model, config).map_err(|error| invalid("worm", error.to_string()))
    }
}

fn build_from_params(params: &Params) -> Result<(IsingGraphWormModel, WormConfig), CarloError> {
    let beta = parse_param::<f64>(params, "beta")?.unwrap_or(1.0);
    let coupling = parse_param::<f64>(params, "J")?.unwrap_or(1.0);
    let pbc = parse_bool(params, "pbc", true)?;
    let lattice = build_lattice_from_params(params, pbc)?;
    let model = IsingGraphWormModel::new(lattice, beta, coupling)
        .map_err(|error| invalid("worm_model", error.to_string()))?;

    let local_updates_per_sweep = parse_param::<usize>(params, "worm_updates_per_sweep")?
        .unwrap_or_else(|| model.lattice().n_edges().max(1));
    let close_probability = parse_param::<f64>(params, "worm_close_probability")?.unwrap_or(0.25);
    if params.contains("worm_fugacity") && params.contains("log_worm_fugacity") {
        return Err(invalid(
            "worm_fugacity",
            "specify either worm_fugacity or log_worm_fugacity, not both",
        ));
    }
    let log_worm_fugacity = if let Some(value) = parse_param::<f64>(params, "log_worm_fugacity")? {
        value
    } else {
        let fugacity = parse_param::<f64>(params, "worm_fugacity")?
            .unwrap_or_else(|| 1.0 / model.lattice().n_sites as f64);
        if !fugacity.is_finite() || fugacity <= 0.0 {
            return Err(invalid("worm_fugacity", "must be finite and positive"));
        }
        fugacity.ln()
    };
    let config = WormConfig {
        local_updates_per_sweep,
        close_probability,
        log_worm_fugacity,
        track_endpoint_pairs: parse_bool(params, "worm_track_endpoint_pairs", false)?,
        cache_audit_interval: parse_param::<u64>(params, "worm_cache_audit_interval")?.unwrap_or(0),
    };
    config
        .validate()
        .map_err(|error| invalid("worm", error.to_string()))?;
    Ok((model, config))
}

fn statistics_to_json(statistics: WormTransitionStatistics) -> Json {
    json!({
        "open_attempts": statistics.open_attempts,
        "open_accepts": statistics.open_accepts,
        "close_attempts": statistics.close_attempts,
        "close_accepts": statistics.close_accepts,
        "step_attempts": statistics.step_attempts,
        "step_accepts": statistics.step_accepts,
        "bounces": statistics.bounces,
        "physical_visits": statistics.physical_visits,
        "worm_visits": statistics.worm_visits,
        "completed_worms": statistics.completed_worms,
        "total_completed_worm_steps": statistics.total_completed_worm_steps,
        "last_completed_worm_steps": statistics.last_completed_worm_steps,
    })
}

fn statistics_from_json(value: &Json) -> Result<WormTransitionStatistics, CarloError> {
    Ok(WormTransitionStatistics {
        open_attempts: required_u64(value, "open_attempts")?,
        open_accepts: required_u64(value, "open_accepts")?,
        close_attempts: required_u64(value, "close_attempts")?,
        close_accepts: required_u64(value, "close_accepts")?,
        step_attempts: required_u64(value, "step_attempts")?,
        step_accepts: required_u64(value, "step_accepts")?,
        bounces: required_u64(value, "bounces")?,
        physical_visits: required_u64(value, "physical_visits")?,
        worm_visits: required_u64(value, "worm_visits")?,
        completed_worms: required_u64(value, "completed_worms")?,
        total_completed_worm_steps: required_u64(value, "total_completed_worm_steps")?,
        last_completed_worm_steps: required_u64(value, "last_completed_worm_steps")?,
    })
}

fn optional_usize(value: &Json, name: &str) -> Result<Option<usize>, CarloError> {
    if value.is_null() {
        Ok(None)
    } else {
        value
            .as_u64()
            .and_then(|number| usize::try_from(number).ok())
            .map(Some)
            .ok_or_else(|| checkpoint_error(format!("{name} must be null or usize")))
    }
}

fn required_u64(value: &Json, name: &str) -> Result<u64, CarloError> {
    value[name]
        .as_u64()
        .ok_or_else(|| checkpoint_error(format!("{name} must be u64")))
}

fn required_usize(value: &Json, name: &str) -> Result<usize, CarloError> {
    required_u64(value, name).and_then(|number| {
        usize::try_from(number).map_err(|_| checkpoint_error(format!("{name} does not fit usize")))
    })
}

fn require_same_f64(value: &Json, name: &str, expected: f64) -> Result<(), CarloError> {
    if value[name].as_f64().map(f64::to_bits) == Some(expected.to_bits()) {
        Ok(())
    } else {
        Err(checkpoint_error(format!("model.{name} mismatch")))
    }
}

fn worm_checkpoint_error(error: WormError) -> CarloError {
    checkpoint_error(error.to_string())
}

fn checkpoint_error(detail: impl Into<String>) -> CarloError {
    CarloError::CheckpointCorrupted {
        detail: detail.into(),
    }
}

fn invalid(field: impl Into<String>, reason: impl Into<String>) -> CarloError {
    CarloError::InvalidConfig {
        field: field.into(),
        reason: reason.into(),
    }
}
