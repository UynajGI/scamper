//! Carlo.rs adapter for continuous-time spin-lattice QMC.

use carlo_rs::{CarloError, Context, FromParams, MonteCarlo, Params};

use crate::algorithm::{QmcKernel, UpdateSchedule};
use crate::graph::{CsrGraph, EdgeSpec};
use crate::local_space::SpinSpace;

use super::configuration::LatticeConfiguration;
use super::error::LatticeQmcError;
use super::model::{EdgeCoupling, GaugePolicy, SiteCoupling, SpinLatticeModel, SpinModelBuilder};
use super::observables::measure_observables;
use super::scattering::ScatteringPolicy;
use super::updates::ContinuousLatticeEngine;

/// Runnable continuous-time spin-lattice simulation.
pub struct LatticeSpinQmc {
    configuration: LatticeConfiguration,
    engine: ContinuousLatticeEngine,
    adaptive_schedule: bool,
    adaptation_interval: usize,
    adaptation_samples: usize,
    adaptation_order_sum: usize,
}

impl LatticeSpinQmc {
    /// Construct from a compiled model.
    pub fn new(
        model: SpinLatticeModel,
        beta: f64,
        schedule: UpdateSchedule,
        rng: &mut rand_xoshiro::Xoshiro256PlusPlus,
    ) -> Result<Self, LatticeQmcError> {
        let configuration = LatticeConfiguration::random(beta, &model, rng)?;
        Ok(Self {
            configuration,
            engine: ContinuousLatticeEngine::new(model, schedule),
            adaptive_schedule: false,
            adaptation_interval: 100,
            adaptation_samples: 0,
            adaptation_order_sum: 0,
        })
    }

    /// Current sampled configuration.
    pub fn configuration(&self) -> &LatticeConfiguration {
        &self.configuration
    }

    /// Update engine and compiled model.
    pub fn engine(&self) -> &ContinuousLatticeEngine {
        &self.engine
    }

    /// Enable warmup-only schedule adaptation.
    pub fn set_adaptive_schedule(&mut self, enabled: bool, interval: usize) {
        self.adaptive_schedule = enabled;
        self.adaptation_interval = interval.max(1);
    }

    fn adapt_during_warmup(&mut self, thermalized: bool) {
        if !self.adaptive_schedule || thermalized {
            return;
        }
        self.adaptation_samples += 1;
        self.adaptation_order_sum += self.configuration.expansion_order();
        if self.adaptation_samples < self.adaptation_interval {
            return;
        }
        let mean_order = self.adaptation_order_sum as f64 / self.adaptation_samples as f64;
        let previous = self.engine.schedule();
        self.engine.set_schedule(UpdateSchedule::new(
            mean_order.ceil().max(1.0) as usize,
            mean_order.sqrt().ceil().max(1.0) as usize,
            previous.max_loop_steps_factor,
        ));
        self.adaptation_samples = 0;
        self.adaptation_order_sum = 0;
    }
}

impl MonteCarlo for LatticeSpinQmc {
    type Rng = rand_xoshiro::Xoshiro256PlusPlus;

    fn sweep(&mut self, ctx: &mut Context<Self::Rng>) {
        self.engine
            .sweep(&mut self.configuration, &mut ctx.rng)
            .unwrap_or_else(|error| panic!("continuous-time lattice sweep failed: {error}"));
        self.adapt_during_warmup(ctx.is_thermalized());
    }

    fn measure(&mut self, ctx: &mut Context<Self::Rng>) {
        let observables = measure_observables(&self.configuration, self.engine.model())
            .unwrap_or_else(|error| panic!("lattice measurement failed: {error}"));
        ctx.measure("MagnetizationZ", observables.magnetization_z);
        ctx.measure("AbsMagnetizationZ", observables.abs_magnetization_z);
        ctx.measure("M2Z", observables.magnetization_z_squared);
        ctx.measure("M4Z", observables.magnetization_z_fourth);
        ctx.measure(
            "StaggeredMagnetizationZ",
            observables.staggered_magnetization_z,
        );
        ctx.measure(
            "StaggeredM2Z",
            observables.staggered_magnetization_z_squared,
        );
        ctx.measure("ChiZRaw", observables.susceptibility_z_raw);
        ctx.measure(
            "StaggeredChiZRaw",
            observables.staggered_susceptibility_z_raw,
        );
        ctx.measure("Energy", observables.energy_total);
        ctx.measure("EnergyPerSite", observables.energy_per_site);
        ctx.measure("ExpansionOrder", observables.expansion_order);
        ctx.measure("DiagonalOrder", observables.diagonal_order);
        ctx.measure("OffDiagonalOrder", observables.offdiagonal_order);
        ctx.measure("VertexDensity", observables.vertex_density);
        ctx.measure(
            "NearestNeighborSzCorrelation",
            observables.nearest_neighbor_sz_correlation,
        );
        ctx.measure("EdgeVertexFraction", observables.edge_vertex_fraction);
        ctx.measure("AverageSign", observables.average_sign);

        let stats = self.engine.stats();
        ctx.measure("DiagonalAcceptance", stats.diagonal_acceptance());
        ctx.measure("MeanLoopSteps", stats.mean_loop_steps());
        ctx.measure("BounceFraction", stats.bounce_fraction());
        ctx.measure("SpatialExitFraction", stats.spatial_exit_fraction());
        ctx.measure("AbortedLoops", stats.aborted_loops as f64);
    }

    fn name(&self) -> &'static str {
        "LatticeSpinQmc"
    }
}

impl FromParams for LatticeSpinQmc {
    fn from_params(params: &Params, rng: &mut Self::Rng) -> Result<Self, CarloError> {
        Self::validate_params(params)?;
        let graph = graph_from_params(params).map_err(to_carlo_error)?;
        let space = spin_space_from_params(params, graph.site_count()).map_err(to_carlo_error)?;
        let model_name = params
            .get::<String>("model")
            .unwrap_or_else(|| "heisenberg".into())
            .to_ascii_lowercase();
        let edge = edge_coupling_from_params(params, &model_name);
        let site = SiteCoupling {
            h_x: params.get::<f64>("h_x").unwrap_or(0.0),
            h_z: params.get::<f64>("h_z").unwrap_or(0.0),
            single_ion: params.get::<f64>("D").unwrap_or(0.0),
            shift: params.get::<f64>("site_shift"),
        };
        let gauge_policy = match params
            .get::<String>("gauge")
            .unwrap_or_else(|| "auto".into())
            .to_ascii_lowercase()
            .as_str()
        {
            "auto" | "marshall" => GaugePolicy::Auto,
            "identity" | "none" => GaugePolicy::Identity,
            other => {
                return Err(CarloError::InvalidConfig {
                    field: "gauge".into(),
                    reason: format!("unsupported gauge policy `{other}`"),
                })
            }
        };
        let scattering_policy = match params
            .get::<String>("scattering")
            .unwrap_or_else(|| "low_bounce".into())
            .to_ascii_lowercase()
            .as_str()
        {
            "low_bounce" | "low-bounce" => ScatteringPolicy::LowBounce,
            "metropolis" => ScatteringPolicy::Metropolis,
            other => {
                return Err(CarloError::InvalidConfig {
                    field: "scattering".into(),
                    reason: format!("unsupported scattering policy `{other}`"),
                })
            }
        };
        let mut builder = SpinModelBuilder::new(graph, space)
            .name(model_name)
            .uniform_edge(edge)
            .uniform_site(site)
            .gauge_policy(gauge_policy)
            .scattering_policy(scattering_policy);
        if let Some(margin) = params.get::<f64>("shift_margin") {
            builder = builder.shift_margin(margin);
        }
        let model = builder.build().map_err(to_carlo_error)?;
        let schedule = UpdateSchedule::new(
            params.get::<usize>("diagonal_proposals").unwrap_or(1),
            params.get::<usize>("directed_loops").unwrap_or(1),
            params.get::<usize>("max_loop_steps_factor").unwrap_or(32),
        );
        let beta = required::<f64>(params, "beta")?;
        let mut simulation = Self::new(model, beta, schedule, rng).map_err(to_carlo_error)?;
        simulation
            .engine
            .set_validate_each_sweep(params.get::<bool>("validate_each_sweep").unwrap_or(false));
        simulation
            .engine
            .set_strict_loop_limits(params.get::<bool>("strict_loop_limits").unwrap_or(false));
        simulation.set_adaptive_schedule(
            params.get::<bool>("adaptive_schedule").unwrap_or(true),
            params.get::<usize>("adaptation_interval").unwrap_or(100),
        );
        Ok(simulation)
    }

    fn validate_params(params: &Params) -> Result<(), CarloError> {
        let beta = required::<f64>(params, "beta")?;
        if !beta.is_finite() || beta <= 0.0 {
            return Err(CarloError::InvalidConfig {
                field: "beta".into(),
                reason: format!("must be finite and positive, got {beta}"),
            });
        }
        Ok(())
    }
}

fn graph_from_params(params: &Params) -> Result<CsrGraph, LatticeQmcError> {
    let topology = params
        .get::<String>("topology")
        .unwrap_or_else(|| "chain".into())
        .to_ascii_lowercase();
    let periodic = params.get::<bool>("pbc").unwrap_or(true);
    match topology.as_str() {
        "chain" => {
            CsrGraph::chain(params.get::<usize>("L").unwrap_or(16), periodic).map_err(Into::into)
        }
        "square" => {
            let width = params
                .get::<usize>("Lx")
                .or_else(|| params.get::<usize>("L"))
                .unwrap_or(8);
            let height = params
                .get::<usize>("Ly")
                .or_else(|| params.get::<usize>("L"))
                .unwrap_or(width);
            CsrGraph::square(width, height, periodic).map_err(Into::into)
        }
        "hypercubic" => {
            let dimensions = parse_usize_csv(params, "dims")?;
            CsrGraph::hypercubic(&dimensions, periodic).map_err(Into::into)
        }
        "edges" | "edge_list" | "edge-list" => {
            let n_sites = required_lattice::<usize>(params, "n_sites")?;
            let encoded = params
                .get::<String>("edges")
                .ok_or_else(|| LatticeQmcError::parameter("edges", "missing edge list"))?;
            CsrGraph::from_edges(n_sites, parse_edge_list(&encoded)?).map_err(Into::into)
        }
        "adjacency" => {
            let encoded = params
                .get::<String>("adjacency")
                .ok_or_else(|| LatticeQmcError::parameter("adjacency", "missing adjacency rows"))?;
            CsrGraph::from_adjacency(&parse_adjacency(&encoded)?).map_err(Into::into)
        }
        other => Err(LatticeQmcError::parameter(
            "topology",
            format!("unsupported topology `{other}`"),
        )),
    }
}

fn edge_coupling_from_params(params: &Params, model: &str) -> EdgeCoupling {
    let mut coupling = match model {
        "heisenberg" => EdgeCoupling::heisenberg(params.get::<f64>("J").unwrap_or(1.0)),
        "xy" => EdgeCoupling::xxz(
            params
                .get::<f64>("J_xy")
                .or_else(|| params.get::<f64>("J"))
                .unwrap_or(1.0),
            0.0,
        ),
        "xxz" => EdgeCoupling::xxz(
            params
                .get::<f64>("J_xy")
                .or_else(|| params.get::<f64>("J"))
                .unwrap_or(1.0),
            params
                .get::<f64>("J_z")
                .or_else(|| params.get::<f64>("J"))
                .unwrap_or(1.0),
        ),
        "xyz" => EdgeCoupling::xyz(
            params.get::<f64>("J_x").unwrap_or(1.0),
            params.get::<f64>("J_y").unwrap_or(1.0),
            params.get::<f64>("J_z").unwrap_or(1.0),
        ),
        "tfim" | "transverse_field_ising" | "transverse-field-ising" => {
            EdgeCoupling::xxz(0.0, params.get::<f64>("J_z").unwrap_or(1.0))
        }
        _ => EdgeCoupling::xyz(
            params.get::<f64>("J_x").unwrap_or(0.0),
            params.get::<f64>("J_y").unwrap_or(0.0),
            params.get::<f64>("J_z").unwrap_or(0.0),
        ),
    };
    coupling.shift = params.get::<f64>("bond_shift");
    coupling
}

fn spin_space_from_params(
    params: &Params,
    site_count: usize,
) -> Result<SpinSpace, LatticeQmcError> {
    if let Some(encoded) = params.get::<String>("two_s_by_site") {
        let values = parse_u16_values(&encoded, "two_s_by_site")?;
        if values.len() != site_count {
            return Err(LatticeQmcError::parameter(
                "two_s_by_site",
                format!("expected {site_count} entries, got {}", values.len()),
            ));
        }
        return SpinSpace::site_resolved(values).map_err(Into::into);
    }
    if let Some(encoded) = params.get::<String>("spins_by_site") {
        let values = encoded
            .split(',')
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| {
                value.parse::<f64>().map_err(|_| {
                    LatticeQmcError::parameter(
                        "spins_by_site",
                        format!("cannot parse `{value}` as a spin"),
                    )
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        if values.len() != site_count {
            return Err(LatticeQmcError::parameter(
                "spins_by_site",
                format!("expected {site_count} entries, got {}", values.len()),
            ));
        }
        let two_s = values
            .into_iter()
            .map(|spin| spin_to_two_s(spin, "spins_by_site"))
            .collect::<Result<Vec<_>, _>>()?;
        return SpinSpace::site_resolved(two_s).map_err(Into::into);
    }
    let two_s = if let Some(two_s) = params.get::<u16>("two_s") {
        if two_s == 0 {
            return Err(LatticeQmcError::parameter("two_s", "must be positive"));
        }
        two_s
    } else {
        spin_to_two_s(params.get::<f64>("spin").unwrap_or(0.5), "spin")?
    };
    SpinSpace::uniform(site_count, two_s).map_err(Into::into)
}

fn spin_to_two_s(spin: f64, field: &str) -> Result<u16, LatticeQmcError> {
    let two_s = 2.0 * spin;
    let rounded = two_s.round();
    if !spin.is_finite()
        || spin <= 0.0
        || (two_s - rounded).abs() > 1.0e-10
        || rounded > f64::from(u16::MAX)
    {
        return Err(LatticeQmcError::parameter(
            field,
            "must contain positive integer or half-integer spins",
        ));
    }
    Ok(rounded as u16)
}

fn parse_u16_values(encoded: &str, field: &str) -> Result<Vec<u16>, LatticeQmcError> {
    encoded
        .split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| {
            value.parse::<u16>().map_err(|_| {
                LatticeQmcError::parameter(field, format!("cannot parse `{value}` as u16"))
            })
        })
        .collect()
}

fn parse_usize_csv(params: &Params, field: &str) -> Result<Vec<usize>, LatticeQmcError> {
    let encoded = params
        .get::<String>(field)
        .ok_or_else(|| LatticeQmcError::parameter(field, "missing comma-separated values"))?;
    encoded
        .split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| {
            value.parse::<usize>().map_err(|_| {
                LatticeQmcError::parameter(field, format!("cannot parse `{value}` as usize"))
            })
        })
        .collect()
}

fn parse_edge_list(encoded: &str) -> Result<Vec<EdgeSpec>, LatticeQmcError> {
    encoded
        .split([',', ';'])
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .map(|token| {
            let normalized = token.replace('-', ":");
            let fields: Vec<_> = normalized.split(':').collect();
            if !(2..=4).contains(&fields.len()) {
                return Err(LatticeQmcError::parameter(
                    "edges",
                    format!("edge `{token}` must be u:v[:kind[:weight]]"),
                ));
            }
            let source = fields[0].parse::<usize>().map_err(|_| {
                LatticeQmcError::parameter("edges", format!("invalid source in `{token}`"))
            })?;
            let target = fields[1].parse::<usize>().map_err(|_| {
                LatticeQmcError::parameter("edges", format!("invalid target in `{token}`"))
            })?;
            let kind = fields
                .get(2)
                .map_or(Ok(0), |value| value.parse::<u16>())
                .map_err(|_| {
                    LatticeQmcError::parameter("edges", format!("invalid kind in `{token}`"))
                })?;
            let weight = fields
                .get(3)
                .map_or(Ok(1.0), |value| value.parse::<f64>())
                .map_err(|_| {
                    LatticeQmcError::parameter("edges", format!("invalid weight in `{token}`"))
                })?;
            Ok(EdgeSpec::typed(source, target, kind, weight))
        })
        .collect()
}

fn parse_adjacency(encoded: &str) -> Result<Vec<Vec<usize>>, LatticeQmcError> {
    encoded
        .split(';')
        .map(|row| {
            row.split(',')
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| {
                    value.parse::<usize>().map_err(|_| {
                        LatticeQmcError::parameter(
                            "adjacency",
                            format!("cannot parse neighbor `{value}`"),
                        )
                    })
                })
                .collect()
        })
        .collect()
}

fn required<T: std::str::FromStr>(params: &Params, field: &str) -> Result<T, CarloError> {
    params
        .get::<T>(field)
        .ok_or_else(|| CarloError::InvalidConfig {
            field: field.into(),
            reason: "missing or unparsable".into(),
        })
}

fn required_lattice<T: std::str::FromStr>(
    params: &Params,
    field: &str,
) -> Result<T, LatticeQmcError> {
    params
        .get::<T>(field)
        .ok_or_else(|| LatticeQmcError::parameter(field, "missing or unparsable"))
}

fn to_carlo_error(error: LatticeQmcError) -> CarloError {
    CarloError::InvalidConfig {
        field: "lattice_qmc".into(),
        reason: error.to_string(),
    }
}
