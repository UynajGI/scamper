import fs from "fs";
import path from "path";

const DIR = "/home/jiangyuan/scuttle/research/rabi_finite_T/results/ed_fss";
const ETAS = [50, 200, 1000, 5000, 25000];
const BETAS = [16, 64, 256, 1024];

const C = {
  reset: "\x1b[0m", bold: "\x1b[1m", dim: "\x1b[2m",
  red: "\x1b[31m", green: "\x1b[32m", yellow: "\x1b[33m",
  blue: "\x1b[34m", cyan: "\x1b[36m",
};

function loadData() {
  const data = {};
  for (const beta of BETAS) {
    const csv = fs.readFileSync(path.join(DIR, `beta_${beta}.csv`), "utf8");
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
  }
  return data;
}

function findMax(rows, key, lambdaMin, lambdaMax) {
  let maxVal = -Infinity, maxLam = NaN;
  for (const r of rows) {
    if (r.lambda < lambdaMin || r.lambda > lambdaMax) continue;
    if (r[key] > maxVal) { maxVal = r[key]; maxLam = r.lambda; }
  }
  return { lambda: maxLam, value: maxVal };
}

function findMin(rows, key, lambdaMin, lambdaMax) {
  let minVal = Infinity, minLam = NaN;
  for (const r of rows) {
    if (r.lambda < lambdaMin || r.lambda > lambdaMax) continue;
    if (r[key] < minVal) { minVal = r[key]; minLam = r.lambda; }
  }
  return { lambda: minLam, value: minVal };
}

function findCrossing(rows, key1, key2, lambdaMin, lambdaMax) {
  for (let i = 1; i < rows.length; i++) {
    const r0 = rows[i-1], r1 = rows[i];
    if (r0.lambda < lambdaMin || r0.lambda > lambdaMax) continue;
    const d0 = r0[key1] - r0[key2];
    const d1 = r1[key1] - r1[key2];
    if (d0 * d1 < 0) {
      return r0.lambda + d0 / (d0 - d1) * (r1.lambda - r0.lambda);
    }
  }
  return NaN;
}

const data = loadData();

// 1. U4(x) — THE KEY QUESTION: do different η curves cross?
console.log(`${C.bold}═══ U4(x) = 1 - ⟨x⁴⟩/(3⟨x²⟩²) — Binder cumulant of position ═══${C.reset}`);
console.log(`${C.dim}Gaussian → 0, double-delta → 2/3. Crossings between η curves would indicate QPT.${C.reset}`);

for (const beta of [1024, 16]) {
  console.log(`\nβ=${beta}:`);
  console.log("λ       η=50     η=200    η=1000   η=5000   η=25000");
  for (const targetL of [0.15, 0.25, 0.35, 0.40, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60, 0.70, 0.90]) {
    const row = data[beta].find(r => Math.abs(r.lambda - targetL) < 0.004);
    if (row) {
      const vals = ETAS.map(eta => {
        const v = row[`u4x_${eta}`];
        const color = v < -0.1 ? C.red : v > 0.3 ? C.green : "";
        return `${color}${v.toFixed(6)}${C.reset}`;
      });
      console.log(`${targetL.toFixed(2)}    ${vals.join("  ")}`);
    }
  }
}

// 2. Check crossings between η=1000 and η=25000
console.log(`\n${C.bold}═══ U4(x) crossings between η pairs (β=1024) ═══${C.reset}`);
const pairs = [[50, 200], [200, 1000], [1000, 5000], [5000, 25000], [50, 25000]];
for (const [a, b] of pairs) {
  const cross = findCrossing(data[1024], `u4x_${a}`, `u4x_${b}`, 0.1, 1.0);
  console.log(`  η=${a} ∩ η=${b}: ${isNaN(cross) ? "NO CROSSING" : `λ* = ${cross.toFixed(4)}`}`);
}

// 3. Energy gap
console.log(`\n${C.bold}═══ Energy gap ΔE₁₂ (β-independent) ═══${C.reset}`);
console.log("λ       η=50     η=200    η=1000   η=5000   η=25000");
for (const targetL of [0.20, 0.30, 0.40, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60, 0.70]) {
  const row = data[1024].find(r => Math.abs(r.lambda - targetL) < 0.004);
  if (row) {
    const vals = ETAS.map(eta => row[`gap_${eta}`].toExponential(2));
    console.log(`${targetL.toFixed(2)}    ${vals.join("   ")}`);
  }
}
const minGap = {};
for (const eta of ETAS) {
  minGap[eta] = findMin(data[1024], `gap_${eta}`, 0.1, 1.0);
}
console.log(`\nMinimum gap location:`);
for (const eta of ETAS) {
  console.log(`  η=${eta.toString().padStart(5)}: λ*=${minGap[eta].lambda.toFixed(4)}, ΔE_min=${minGap[eta].value.toExponential(3)}`);
}

// 4. Specific heat — β-dependent peak
console.log(`\n${C.bold}═══ Specific heat C_V peak position λ*(β) ═══${C.reset}`);
console.log("This is THE finite-T diagnostic. Peak shifts with β.");
console.log("η\\β      16        64        256       1024      (b16-b1024)");
for (const eta of ETAS) {
  const peaks = BETAS.map(b => findMax(data[b], `cv_${eta}`, 0.1, 1.0).lambda);
  const shift = peaks[0] - peaks[3];
  const color = Math.abs(shift) > 0.02 ? C.cyan : C.dim;
  console.log(`${eta.toString().padStart(5)}    ${peaks.map(p => p.toFixed(4)).join("   ")}    ${color}${shift >= 0 ? "+" : ""}${shift.toFixed(4)}${C.reset}`);
}

// 5. ⟨σz⟩
console.log(`\n${C.bold}═══ ⟨σz⟩ spin polarization (β=1024) ═══${C.reset}`);
console.log("λ       η=50     η=200    η=1000   η=5000   η=25000");
for (const targetL of [0.20, 0.30, 0.40, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60, 0.70]) {
  const row = data[1024].find(r => Math.abs(r.lambda - targetL) < 0.004);
  if (row) {
    const vals = ETAS.map(eta => row[`sigmaz_${eta}`].toFixed(6));
    console.log(`${targetL.toFixed(2)}    ${vals.join("  ")}`);
  }
}

// 6. ⟨x²⟩
console.log(`\n${C.bold}═══ ⟨x²⟩ position fluctuation (β=1024) ═══${C.reset}`);
console.log("λ       η=50     η=200    η=1000   η=5000   η=25000");
for (const targetL of [0.20, 0.35, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60, 0.70]) {
  const row = data[1024].find(r => Math.abs(r.lambda - targetL) < 0.004);
  if (row) {
    const vals = ETAS.map(eta => row[`x2_${eta}`].toExponential(3));
    console.log(`${targetL.toFixed(2)}    ${vals.join("  ")}`);
  }
}
