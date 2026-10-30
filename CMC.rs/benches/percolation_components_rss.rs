#[allow(dead_code)]
mod percolation_support;

use std::fs;

use cmc_rs::{
    analyze, BondBernoulli, BoundaryQuery, ComponentWorkspace, GraphView, MixedBernoulli,
    Probability, SiteBernoulli, StaticConfiguration, UndirectedGraphView,
};
use percolation_support::{case, P_BOND, P_SITE};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const SAMPLES: usize = 128;

fn kibibytes(field: &str) -> Result<usize, String> {
    let status = fs::read_to_string("/proc/self/status")
        .map_err(|error| format!("could not read /proc/self/status: {error}"))?;
    let line = status
        .lines()
        .find(|line| line.starts_with(field))
        .ok_or_else(|| format!("missing {field} in /proc/self/status"))?;
    let mut parts = line.split_whitespace();
    let _name = parts.next();
    parts
        .next()
        .ok_or_else(|| format!("missing value for {field}"))?
        .parse::<usize>()
        .map_err(|error| format!("invalid {field} value: {error}"))
}

fn probability(value: f64) -> Probability {
    Probability::new(value).expect("probe probability")
}

fn main() -> Result<(), String> {
    let case_name = std::env::var("SCUTTLE_RSS_CASE").unwrap_or_else(|_| "square".to_string());
    let law_name = std::env::var("SCUTTLE_RSS_LAW").unwrap_or_else(|_| "mixed".to_string());
    let law = match law_name.as_str() {
        "site" | "bond" | "mixed" => law_name.as_str(),
        _ => return Err(format!("unknown SCUTTLE_RSS_LAW={law_name}")),
    };
    let before_rss = kibibytes("VmRSS:")?;
    let before_hwm = kibibytes("VmHWM:")?;

    let case = case(&case_name).ok_or_else(|| format!("unknown SCUTTLE_RSS_CASE={case_name}"))?;
    let vertex_count = case.lattice.vertex_count();
    let edge_count = case.lattice.edge_count();
    let query = BoundaryQuery::new(&case.lattice, &case.from, &case.to)
        .map_err(|error| error.to_string())?;
    let queries = [query];
    let mut configuration = StaticConfiguration::new(vertex_count, edge_count);
    let mut workspace = ComponentWorkspace::new();
    workspace
        .prepare(vertex_count, queries.len(), false)
        .map_err(|error| error.to_string())?;
    let site = SiteBernoulli::new(probability(P_SITE).into());
    let bond = BondBernoulli::new(probability(P_BOND).into());
    let mixed = MixedBernoulli::new(probability(P_SITE).into(), probability(P_BOND).into());
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(0x4633_5253_5350);

    for _ in 0..8 + SAMPLES {
        match law {
            "site" => site.sample(&case.lattice, &mut configuration, &mut rng),
            "bond" => bond.sample(&case.lattice, &mut configuration, &mut rng),
            "mixed" => mixed.sample(&case.lattice, &mut configuration, &mut rng),
            _ => unreachable!("validated law"),
        }
        .map_err(|error| error.to_string())?;
        std::hint::black_box(
            analyze(
                &case.lattice,
                &configuration,
                &configuration,
                &queries,
                &mut workspace,
            )
            .map_err(|error| error.to_string())?,
        );
    }

    let after_rss = kibibytes("VmRSS:")?;
    let after_hwm = kibibytes("VmHWM:")?;
    eprintln!(
        "PERCOLATION_COMPONENT_RSS topology={} law={} V={} E={} measured_samples={} \
         before_vm_rss_kib={} after_vm_rss_kib={} vm_rss_delta_kib={} \
         before_vm_hwm_kib={} after_vm_hwm_kib={} vm_hwm_delta_kib={} \
         metric=linux_proc_process_lifetime_high_water",
        case.name,
        law_name,
        vertex_count,
        edge_count,
        SAMPLES,
        before_rss,
        after_rss,
        after_rss.saturating_sub(before_rss),
        before_hwm,
        after_hwm,
        after_hwm.saturating_sub(before_hwm),
    );
    Ok(())
}
