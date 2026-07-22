import fs from "fs";
import path from "path";

const DIR = "/home/jiangyuan/scuttle/research/rabi_finite_T/results/ed_fss";
const ETAS = [50, 200, 1000];
const BETAS = [0.5, 1.0, 2.0, 4.0, 16.0, 64.0, 256.0, 1024.0];

const C = {
  reset: "\x1b[0m", bold: "\x1b[1m", dim: "\x1b[2m",
  red: "\x1b[31m", green: "\x1b[32m", yellow: "\x1b[33m",
  blue: "\x1b[34m", cyan: "\x1b[36m", magenta: "\x1b[35m",
};

function loadData() {
  const data = {};
  for (const beta of BETAS) {
    const fname = `beta_${beta}.csv`;
    try {
      const csv = fs.readFileSync(path.join(DIR, fname), "utf8");
      const lines = csv.trim().split("\n");
      data[beta] = lines.slice(1).map(line => {
        const p = line.split(",").map(parseFloat);
        const row = { lambda: p[0] };
        ETAS.forEach((eta, i) => {
          const base = 1 + i * 6;
          row[`ntilde_${eta}`] = p[base];
          row[`sigmaz_${eta}`] = p[base + 1];
          row[`x2_${eta}`] = p[base + 2];
          row[`u4x_${eta}`] = p[base + 3];
          row[`gap_${eta}`] = p[base + 4];
          row[`cv_${eta}`] = p[base + 5];
        });
        return row;
      });
    } catch (e) { data[beta] = null; }
  }
  return data;
}

function findMax(rows, key, lo, hi) {
  let mx = -Infinity, ml = NaN;
  for (const r of rows) {
    if (!r || r.lambda < lo || r.lambda > hi) continue;
    if (r[key] > mx) { mx = r[key]; ml = r.lambda; }
  }
  return { lambda: ml, value: mx };
}

const data = loadData();

// 1. Gap closure — where does the QPT actually happen?
console.log(`${C.bold}═══ Energy gap ΔE₁₂ vs λ (β-independent) ═══${C.reset}`);
console.log("Gap = Δ×exp(-2λ²/η) approx. Closes when λ²/η >> 1.\n");
console.log("λ      η=50     η=200    η=1000");
const beta1024 = data[1024];
if (beta1024) {
  for (const targetL of [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]) {
    const row = beta1024.find(r => Math.abs(r.lambda - targetL) < 0.03);
    if (row) {
      const gaps = ETAS.map(eta => {
        const g = row[`gap_${eta}`];
        if (g < 0.01) return `${C.red}${g.toExponential(2)}${C.reset}`;
        if (g < 0.1) return `${C.yellow}${g.toExponential(2)}${C.reset}`;
        return g.toExponential(2);
      });
      console.log(`${targetL.toString().padStart(4)}   ${gaps.join("   ")}`);
    }
  }
}

// 2. U4(x) at large λ — does it finally become non-zero?
console.log(`\n${C.bold}═══ U4(x) at large λ (β=1024) ═══${C.reset}`);
console.log("λ      η=50        η=200       η=1000");
if (beta1024) {
  for (const targetL of [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]) {
    const row = beta1024.find(r => Math.abs(r.lambda - targetL) < 0.03);
    if (row) {
      const u4s = ETAS.map(eta => {
        const u = row[`u4x_${eta}`];
        const color = u > 0.1 ? C.green : u > 0.01 ? C.yellow : "";
        return `${color}${u.toFixed(6)}${C.reset}`;
      });
      console.log(`${targetL.toString().padStart(4)}   ${u4s.join("   ")}`);
    }
  }
}

// 3. U4 crossings between η=50 and η=200
console.log(`\n${C.bold}═══ U4(x) crossings: η=50 ∩ η=200 (β=1024) ═══${C.reset}`);
if (beta1024) {
  for (let i = 1; i < beta1024.length; i++) {
    const r0 = beta1024[i-1], r1 = beta1024[i];
    const d0 = r0.u4x_50 - r0.u4x_200;
    const d1 = r1.u4x_50 - r1.u4x_200;
    if (d0 * d1 < 0 && Math.abs(r0.u4x_50) > 1e-6) {
      const cross = r0.lambda + d0 / (d0 - d1) * (r1.lambda - r0.lambda);
      console.log(`  CROSSING at λ* ≈ ${cross.toFixed(4)}  (U4_50=${r0.u4x_50.toFixed(6)}, U4_200=${r0.u4x_200.toFixed(6)})`);
    }
  }
}

// 4. Specific heat peaks — THE finite-T diagnostic
console.log(`\n${C.bold}═══ C_V peak position λ*(β) — the finite-T shift ═══${C.reset}`);
console.log("η\\β     0.5       1.0       2.0       4.0       16        64        256       1024");
for (const eta of ETAS) {
  const peaks = BETAS.map(b => {
    if (!data[b]) return "N/A";
    return findMax(data[b], `cv_${eta}`, 0.5, 20).lambda.toFixed(2);
  });
  console.log(`${eta.toString().padStart(4)}   ${peaks.join("      ")}`);
}

// 5. C_V peak values
console.log(`\n${C.bold}═══ C_V peak heights ═══${C.reset}`);
console.log("η\\β     0.5       1.0       2.0       4.0       16        64        256       1024");
for (const eta of ETAS) {
  const peaks = BETAS.map(b => {
    if (!data[b]) return "N/A";
    return findMax(data[b], `cv_${eta}`, 0.5, 20).value.toExponential(2);
  });
  console.log(`${eta.toString().padStart(4)}   ${peaks.join("  ")}`);
}

// 6. ⟨σz⟩ — does the spin "flip" at large λ?
console.log(`\n${C.bold}═══ ⟨σz⟩ at large λ (β=1024) ═══${C.reset}`);
console.log("λ      η=50        η=200       η=1000");
if (beta1024) {
  for (const targetL of [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]) {
    const row = beta1024.find(r => Math.abs(r.lambda - targetL) < 0.03);
    if (row) {
      const szs = ETAS.map(eta => row[`sigmaz_${eta}`].toFixed(6));
      console.log(`${targetL.toString().padStart(4)}   ${szs.join("   ")}`);
    }
  }
}

// 7. ñ at large λ
console.log(`\n${C.bold}═══ ñ = η⟨n⟩ at large λ (β=1024) ═══${C.reset}`);
console.log("BO prediction: ñ → λ². λ      η=50/λ²    η=200/λ²   η=1000/λ²");
if (beta1024) {
  for (const targetL of [1, 2, 3, 5, 8, 10, 15, 20]) {
    const row = beta1024.find(r => Math.abs(r.lambda - targetL) < 0.03);
    if (row) {
      const ratios = ETAS.map(eta => (row[`ntilde_${eta}`] / (targetL * targetL)).toFixed(4));
      console.log(`${targetL.toString().padStart(4)}   ${ratios.join("      ")}`);
    }
  }
}
