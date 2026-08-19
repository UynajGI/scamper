//! Carlo.rs adapter and versioned JSON snapshot for the Ising graph worm.

use super::{
    EndpointPairHistogram, IsingGraphConfiguration, IsingGraphWormEnsemble, IsingGraphWormModel,
    WormConfig, WormError, WormKernel, WormSector, WormState, WormTransitionStatistics,
};
use crate::classical_mc::{build_lattice_from_params, parse_bool, parse_param};
use crate::lattice::graph::CsrLattice;
use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, Params};
use serde_json::{json, Value as Json};

/// Scheduler-ready persistent high-temperature Ising graph worm.
///
/// Multi-component (disconnected, possibly with isolated-site) lattices are
/// supported: the high-temperature graph ensemble factorizes over connected
/// components, so the adapter runs one independent two-defect worm per
/// component on domain-separated derived streams and combines observables
/// additively. A connected lattice is the single-component special case with
/// the historical single-kernel API ([`IsingGraphWormMC::kernel`]).
pub struct IsingGraphWormMC {
    ensemble: IsingGraphWormEnsemble,
    endpoint_sample: Vec<f64>,
}

impl IsingGraphWormMC {
    /// Build from a single-component [`IsingGraphWormModel`].
    ///
    /// The model constructor guarantees a connected lattice, so this always
    /// yields a one-component ensemble with the historical kernel accessors.
    pub fn new(model: IsingGraphWormModel, config: WormConfig) -> Result<Self, WormError> {
        let beta = model.beta();
        let coupling = model.coupling();
        let lattice = model.lattice().clone();
        Self::from_lattice(lattice, beta, coupling, config)
    }

    /// Build from an arbitrary (possibly disconnected) lattice.
    ///
    /// The lattice is decomposed into connected components and one
    /// independently driven worm is constructed per component.
    pub fn from_lattice(
        lattice: CsrLattice,
        beta: f64,
        coupling: f64,
        config: WormConfig,
    ) -> Result<Self, WormError> {
        let n_sites = lattice.n_sites;
        let ensemble = IsingGraphWormEnsemble::new(lattice, beta, coupling, config)?;
        let endpoint_sample = if ensemble.config().track_endpoint_pairs {
            vec![0.0; n_sites.saturating_mul(n_sites)]
        } else {
            Vec::new()
        };
        Ok(Self {
            ensemble,
            endpoint_sample,
        })
    }

    /// The kernel of the sole component (single-component lattices only).
    #[inline]
    pub fn kernel(&self) -> &WormKernel<IsingGraphWormModel> {
        self.ensemble.single_kernel()
    }

    /// Mutable counterpart of [`Self::kernel`].
    #[inline]
    pub fn kernel_mut(&mut self) -> &mut WormKernel<IsingGraphWormModel> {
        self.ensemble.single_kernel_mut()
    }

    /// Number of connected components being sampled.
    #[inline]
    pub fn n_components(&self) -> usize {
        self.ensemble.n_components()
    }

    #[inline]
    pub fn ensemble(&self) -> &IsingGraphWormEnsemble {
        &self.ensemble
    }

    #[inline]
    pub fn ensemble_mut(&mut self) -> &mut IsingGraphWormEnsemble {
        &mut self.ensemble
    }

    /// Worm-estimated two-point correlation `⟨s_tail s_head⟩`; both sites must
    /// belong to the same connected component.
    #[inline]
    pub fn endpoint_correlation(&self, tail: usize, head: usize) -> Option<f64> {
        self.ensemble.endpoint_correlation(tail, head)
    }

    /// Versioned model/kernel snapshot. RNG and measurement state remain owned
    /// by Carlo.rs's context checkpoint.
    ///
    /// Layout: single-component ensembles keep the historical `v1` layout
    /// (model/state/runtime at the top level). Multi-component ensembles use
    /// `v2` with a `components` array; the loader accepts both.
    pub fn save_snapshot(&self) -> Json {
        let worm_config = self.ensemble.config();
        let model = json!({
            "beta": self.ensemble.beta(),
            "coupling": self.ensemble.coupling(),
            "n_sites": self.ensemble.lattice().n_sites,
            "n_edges": self.ensemble.lattice().n_edges(),
            "edges": self.ensemble.lattice().edges.iter().map(|edge| json!({
                "source": edge.source,
                "target": edge.target,
                "kind": edge.kind.as_label(),
                "weight": edge.weight,
            })).collect::<Vec<_>>(),
        });
        let config = json!({
            "local_updates_per_sweep": worm_config.local_updates_per_sweep,
            "close_probability": worm_config.close_probability,
            "log_worm_fugacity": worm_config.log_worm_fugacity,
            "track_endpoint_pairs": worm_config.track_endpoint_pairs,
            "cache_audit_interval": worm_config.cache_audit_interval,
        });

        if self.ensemble.n_components() == 1 {
            let component = &self.ensemble.components()[0];
            let state = component.kernel().state();
            json!({
                "format": "cmc-rs-ising-worm-v1",
                "model": model,
                "config": config,
                "state": {
                    "sector": state.sector().as_str(),
                    "head": state.head(),
                    "tail": state.tail(),
                    "occupied": state.configuration().occupied(),
                },
                "runtime": {
                    "sweeps": component.kernel().sweeps(),
                    "current_worm_steps": component.kernel().current_worm_steps(),
                    "statistics": statistics_to_json(*component.kernel().statistics()),
                    "endpoint_pairs": component.kernel().endpoint_pairs().map(|histogram| {
                        json!({
                            "bins": histogram.bins(),
                            "counts": histogram.counts(),
                            "samples": histogram.samples(),
                        })
                    }),
                },
            })
        } else {
            let components = self
                .ensemble
                .components()
                .iter()
                .map(|component| {
                    let state = component.kernel().state();
                    json!({
                        "sites": component.global_sites(),
                        "state": {
                            "sector": state.sector().as_str(),
                            "head": state.head(),
                            "tail": state.tail(),
                            "occupied": state.configuration().occupied(),
                        },
                        "runtime": {
                            "sweeps": component.kernel().sweeps(),
                            "current_worm_steps": component.kernel().current_worm_steps(),
                            "statistics": statistics_to_json(*component.kernel().statistics()),
                            "endpoint_pairs": component.kernel().endpoint_pairs().map(|histogram| {
                                json!({
                                    "bins": histogram.bins(),
                                    "counts": histogram.counts(),
                                    "samples": histogram.samples(),
                                })
                            }),
                        },
                    })
                })
                .collect::<Vec<_>>();
            json!({
                "format": "cmc-rs-ising-worm-v2",
                "model": model,
                "config": config,
                "components": components,
            })
        }
    }

    pub fn load_snapshot(&mut self, snapshot: &Json) -> Result<(), CarloError> {
        let format = snapshot["format"].as_str();
        match format {
            Some("cmc-rs-ising-worm-v1") => self.load_snapshot_v1(snapshot),
            Some("cmc-rs-ising-worm-v2") => self.load_snapshot_v2(snapshot),
            _ => Err(checkpoint_error("unknown Ising worm snapshot format")),
        }
    }

    fn load_snapshot_v1(&mut self, snapshot: &Json) -> Result<(), CarloError> {
        if self.ensemble.n_components() != 1 {
            return Err(checkpoint_error(
                "single-component (v1) snapshots cannot restore a multi-component ensemble",
            ));
        }
        self.validate_snapshot_model_and_config(snapshot)?;

        let component_state =
            parse_component_state(self.ensemble.single_kernel().model(), &snapshot["state"])?;
        let runtime = &snapshot["runtime"];
        let statistics = statistics_from_json(&runtime["statistics"])?;
        let sweeps = required_u64(runtime, "sweeps")?;
        let current_worm_steps = required_u64(runtime, "current_worm_steps")?;
        let endpoint_pairs = parse_endpoint_pairs(&runtime["endpoint_pairs"])?;

        self.ensemble
            .single_kernel_mut()
            .restore_runtime(
                component_state,
                statistics,
                endpoint_pairs,
                sweeps,
                current_worm_steps,
            )
            .map_err(worm_checkpoint_error)
    }

    fn load_snapshot_v2(&mut self, snapshot: &Json) -> Result<(), CarloError> {
        self.validate_snapshot_model_and_config(snapshot)?;
        let components = snapshot["components"]
            .as_array()
            .ok_or_else(|| checkpoint_error("components must be an array"))?;
        if components.len() != self.ensemble.n_components() {
            return Err(checkpoint_error(
                "Ising worm snapshot component-count mismatch",
            ));
        }
        for (index, entry) in components.iter().enumerate() {
            let expected_sites = self.ensemble.components()[index].global_sites();
            let sites = entry["sites"]
                .as_array()
                .ok_or_else(|| checkpoint_error("components[].sites must be an array"))?;
            if sites.len() != expected_sites.len()
                || sites
                    .iter()
                    .zip(expected_sites.iter())
                    .any(|(value, &site)| value.as_u64() != Some(site as u64))
            {
                return Err(checkpoint_error(
                    "Ising worm snapshot component site partition mismatch",
                ));
            }
            let state =
                parse_component_state(self.ensemble.components()[index].model(), &entry["state"])?;
            let runtime = &entry["runtime"];
            let statistics = statistics_from_json(&runtime["statistics"])?;
            let sweeps = required_u64(runtime, "sweeps")?;
            let current_worm_steps = required_u64(runtime, "current_worm_steps")?;
            let endpoint_pairs = parse_endpoint_pairs(&runtime["endpoint_pairs"])?;
            self.ensemble.components_mut()[index]
                .kernel_mut()
                .restore_runtime(
                    state,
                    statistics,
                    endpoint_pairs,
                    sweeps,
                    current_worm_steps,
                )
                .map_err(worm_checkpoint_error)?;
        }
        Ok(())
    }

    /// Validate the global model/config sections against the current ensemble.
    fn validate_snapshot_model_and_config(&self, snapshot: &Json) -> Result<(), CarloError> {
        let model = &snapshot["model"];
        let config = &snapshot["config"];
        let lattice = self.ensemble.lattice();
        let worm_config = self.ensemble.config();
        require_same_f64(model, "beta", self.ensemble.beta())?;
        require_same_f64(model, "coupling", self.ensemble.coupling())?;
        if required_usize(model, "n_sites")? != lattice.n_sites
            || required_usize(model, "n_edges")? != lattice.n_edges()
        {
            return Err(checkpoint_error(
                "Ising worm snapshot topology size mismatch",
            ));
        }
        let edges = model["edges"]
            .as_array()
            .ok_or_else(|| checkpoint_error("model.edges must be an array"))?;
        if edges.len() != lattice.edges.len() {
            return Err(checkpoint_error("Ising worm snapshot edge-count mismatch"));
        }
        for (value, edge) in edges.iter().zip(&lattice.edges) {
            if required_usize(value, "source")? != edge.source
                || required_usize(value, "target")? != edge.target
                || value["kind"].as_str() != Some(edge.kind.as_label())
                || value["weight"].as_f64().map(f64::to_bits) != Some(edge.weight.to_bits())
            {
                return Err(checkpoint_error("Ising worm snapshot edge mismatch"));
            }
        }

        if required_usize(config, "local_updates_per_sweep")? != worm_config.local_updates_per_sweep
            || config["close_probability"].as_f64().map(f64::to_bits)
                != Some(worm_config.close_probability.to_bits())
            || config["log_worm_fugacity"].as_f64().map(f64::to_bits)
                != Some(worm_config.log_worm_fugacity.to_bits())
            || config["track_endpoint_pairs"].as_bool() != Some(worm_config.track_endpoint_pairs)
            || required_u64(config, "cache_audit_interval")? != worm_config.cache_audit_interval
        {
            return Err(checkpoint_error(
                "Ising worm snapshot transition configuration mismatch",
            ));
        }
        Ok(())
    }
}

/// Parse one component's worm state (occupied is indexed by component edges).
fn parse_component_state(
    model: &IsingGraphWormModel,
    state: &Json,
) -> Result<WormState<IsingGraphConfiguration, usize>, CarloError> {
    let occupied = state["occupied"]
        .as_array()
        .ok_or_else(|| checkpoint_error("state.occupied must be an array"))?
        .iter()
        .map(|value| {
            value
                .as_bool()
                .ok_or_else(|| checkpoint_error("state.occupied contains a non-boolean"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let configuration =
        IsingGraphConfiguration::from_occupied(model, occupied).map_err(worm_checkpoint_error)?;
    let sector = WormSector::parse(
        state["sector"]
            .as_str()
            .ok_or_else(|| checkpoint_error("state.sector is missing"))?,
    )
    .map_err(worm_checkpoint_error)?;
    let head = optional_usize(&state["head"], "state.head")?;
    let tail = optional_usize(&state["tail"], "state.tail")?;
    WormState::from_parts(configuration, sector, head, tail).map_err(worm_checkpoint_error)
}

fn parse_endpoint_pairs(value: &Json) -> Result<Option<EndpointPairHistogram>, CarloError> {
    if value.is_null() {
        return Ok(None);
    }
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
    let histogram =
        EndpointPairHistogram::from_counts(bins, counts, samples).map_err(worm_checkpoint_error)?;
    Ok(Some(histogram))
}

impl MonteCarlo for IsingGraphWormMC {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, context: &mut Context<Self::Rng>) {
        self.ensemble
            .sweep(&mut context.rng)
            .expect("Ising graph worm sweep failed");
    }

    fn measure(&mut self, context: &mut Context<Self::Rng>) {
        // Aggregate last-sweep transition diagnostics over all components.
        let mut step_attempts = 0u64;
        let mut step_accepts = 0u64;
        let mut open_attempts = 0u64;
        let mut open_accepts = 0u64;
        let mut close_attempts = 0u64;
        let mut close_accepts = 0u64;
        let mut completed_worm_length = 0u64;
        let mut worm_components = 0usize;
        for component in self.ensemble.components() {
            let last = *component.kernel().last_sweep_statistics();
            step_attempts += last.step_attempts;
            step_accepts += last.step_accepts;
            open_attempts += last.open_attempts;
            open_accepts += last.open_accepts;
            close_attempts += last.close_attempts;
            close_accepts += last.close_accepts;
            completed_worm_length += last.last_completed_worm_steps;
            worm_components += usize::from(component.kernel().state().is_worm());
        }
        let n_components = self.ensemble.n_components() as f64;
        context.measure("WormSector", worm_components as f64 / n_components);
        context.measure(
            "PhysicalSector",
            (self.ensemble.n_components() - worm_components) as f64 / n_components,
        );
        context.measure("WormStepAcceptance", ratio(step_accepts, step_attempts));
        context.measure("WormOpenAcceptance", ratio(open_accepts, open_attempts));
        context.measure("WormCloseAcceptance", ratio(close_accepts, close_attempts));
        context.measure("CompletedWormLength", completed_worm_length as f64);

        if let Some((energy, occupied_edges)) = self.ensemble.total_energy_and_occupied_edges() {
            // All components physical: the conditioned measure is exactly the
            // product of the per-component even-subgraph ensembles.
            let n_edges = self.ensemble.lattice().n_edges().max(1) as f64;
            context.measure("GraphOccupiedEdges", occupied_edges as f64);
            context.measure("GraphEdgeDensity", occupied_edges as f64 / n_edges);
            context.measure("Energy", energy);
        } else {
            // Endpoint observables from the first worm-sector component
            // (single-component ensembles: the only one, historical semantics).
            let mut head_tail: Option<(usize, usize)> = None;
            self.endpoint_sample.fill(0.0);
            let n_sites = self.ensemble.lattice().n_sites;
            let tracking = !self.endpoint_sample.is_empty();
            for component in self.ensemble.components() {
                let state = component.kernel().state();
                if !state.is_worm() {
                    continue;
                }
                let head = *state.head().expect("worm sector has a head");
                let tail = *state.tail().expect("worm sector has a tail");
                if tracking {
                    let global_tail = component.global_sites()[tail];
                    let global_head = component.global_sites()[head];
                    self.endpoint_sample[global_tail * n_sites + global_head] = 1.0;
                }
                if head_tail.is_none() {
                    head_tail = Some((
                        component.global_sites()[head],
                        component.global_sites()[tail],
                    ));
                }
            }
            if let Some((head, tail)) = head_tail {
                context.measure("WormHead", head as f64);
                context.measure("WormTail", tail as f64);
            }
            if tracking {
                context.measure_array("WormEndpointPairs", &self.endpoint_sample);
            }
        }
    }

    fn name(&self) -> &'static str {
        "IsingHighTemperatureGraphWorm"
    }
}

#[inline]
fn ratio(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

impl FromParams for IsingGraphWormMC {
    fn validate_params(params: &Params) -> Result<(), CarloError> {
        // Construct the full ensemble: the per-component model constructors
        // are the single source of truth for beta/coupling/edge validity.
        let (lattice, beta, coupling, config) = build_from_params(params)?;
        IsingGraphWormMC::from_lattice(lattice, beta, coupling, config)
            .map(|_| ())
            .map_err(|error| invalid("worm_model", error.to_string()))
    }

    fn from_params(params: &Params, _rng: &mut Self::Rng) -> Result<Self, CarloError> {
        let (lattice, beta, coupling, config) = build_from_params(params)?;
        Self::from_lattice(lattice, beta, coupling, config)
            .map_err(|error| invalid("worm_model", error.to_string()))
    }
}

fn build_from_params(params: &Params) -> Result<(CsrLattice, f64, f64, WormConfig), CarloError> {
    let beta = parse_param::<f64>(params, "beta")?.unwrap_or(1.0);
    let coupling = parse_param::<f64>(params, "J")?.unwrap_or(1.0);
    let pbc = parse_bool(params, "pbc", true)?;
    let lattice = build_lattice_from_params(params, pbc)?;

    let local_updates_per_sweep = parse_param::<usize>(params, "worm_updates_per_sweep")?
        .unwrap_or_else(|| {
            // Per component: one update stream per physical edge, matching the
            // connected-lattice default.
            lattice.n_edges().max(1)
        });
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
            .unwrap_or_else(|| 1.0 / lattice.n_sites as f64);
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
    Ok((lattice, beta, coupling, config))
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
