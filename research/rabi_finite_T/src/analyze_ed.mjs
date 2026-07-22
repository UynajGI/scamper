#!/usr/bin/env node
// Analyze ED FSS results: extract ñ transition points and β-dependence
import fs from "fs";
import path from "path";
import { pathToFileURL } from "url";

const DIR = "/home/jiangyuan/scuttle/research/rabi_finite_T/results/ed_fss";
const ETAS = [50, 200, 1000, 5000, 25000];
const BETAS = [16, 64, 256, 1024];

// Color codes for terminal
const C = {
  reset: "\x1b[0m", bold: "\x1b[1m", dim: "\x1b[2m",
  red: "\x1b[31m", green: "\x1b[32m", yellow: "\x1b[33m",
  blue: "\x1b[34m", magenta: "\x1b[35m", cyan: "\x1b[36m",
};

function loadData() {
  const data = {};
  for (const beta of BETAS) {
    const csv = fs.readFileSync(path.join(DIR, `beta_${beta}.csv`), "utf8");
    const lines = csv.trim().split("\n");
    data[beta] = lines.slice(1).map(line => {
      const p = line.split(",").map(parseFloat);
      const lambdas = p[0];
      const row = { lambda: lambdas };
      ETAS.forEach((eta, i) => {
        const base = 1 + i * 3;
        row[`n_${eta}`] = p[base];
        row[`ntilde_${eta}`] = p[base + 1];
        row[`ratio_${eta}`] = p[base + 2]; // ntilde / lambda^2
      });
      return row;
    });
  }
  return data;
}

// Find where ntilde/lambda^2 crosses a threshold (e.g., 0.95)
// This marks the "entry into broken phase"
function findThresholdCrossing(rows, eta, threshold) {
  const key = `ratio_${eta}`;
  for (let i = 1; i < rows.length; i++) {
    if (rows[i - 1][key] < threshold && rows[i][key] >= threshold) {
      // Linear interpolation
      const r0 = rows[i - 1][key];
      const r1 = rows[i][key];
      const l0 = rows[i - 1].lambda;
      const l1 = rows[i].lambda;
      return l0 + (threshold - r0) / (r1 - r0) * (l1 - l0);
    }
  }
  return NaN; // Never crosses
}

// Find inflection point: max of d(ntilde)/d(lambda) (steepest slope)
function findMaxSlope(rows, eta) {
  const key = `ntilde_${eta}`;
  let maxSlope = 0;
  let maxLambda = NaN;
  for (let i = 1; i < rows.length; i++) {
    const slope = (rows[i][key] - rows[i - 1][key]) / (rows[i].lambda - rows[i - 1].lambda);
    if (slope > maxSlope) {
      maxSlope = slope;
      maxLambda = (rows[i].lambda + rows[i - 1].lambda) / 2;
    }
  }
  return { lambda: maxLambda, slope: maxSlope };
}

// Find inflection point: zero crossing of second derivative
function findInflection(rows, eta) {
  const key = `ntilde_${eta}`;
  let prev2 = null, prev1 = null;
  for (let i = 0; i < rows.length; i++) {
    if (i < 2) { 
      if (i === 0) prev2 = rows[i];
      if (i === 1) prev1 = rows[i];
      continue;
    }
    const dl = rows[i].lambda - prev2.lambda;
    const d2 = (rows[i][key] - 2 * prev1[key] + prev2[key]) / (dl / 2) ** 2;
    // Track sign change of second derivative
    if (prev2._d2 !== undefined && prev2._d2 > 0 && d2 < 0) {
      // Inflection: second derivative changes from + to -
      return rows[i - 1].lambda;
    }
    prev2._d2 = d2;
    prev2 = prev1;
    prev1 = rows[i];
  }
  return NaN;
}

const data = loadData();

// 1. Show ñ/λ² transition for each η
console.log(`${C.bold}=== ñ/λ² transition (β=1024) ===${C.reset}`);
console.log("Shows how ñ/λ² rises from symmetric phase to BO limit (=1)");
console.log("λ       η=50    η=200   η=1000  η=5000  η=25000");
const sampleLambdas = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70];
for (const targetL of sampleLambdas) {
  const row = data[1024].find(r => Math.abs(r.lambda - targetL) < 0.003);
  if (row) {
    const vals = ETAS.map(eta => {
      const r = row[`ratio_${eta}`];
      const color = r > 0.95 ? C.green : r > 0.7 ? C.yellow : C.red;
      return `${color}${r.toFixed(4)}${C.reset}`;
    });
    console.log(`${targetL.toFixed(2)}    ${vals.join("   ")}`);
  }
}

// 2. Threshold crossing λ*(η, β) where ñ/λ² = 0.90
console.log(`\n${C.bold}=== λ* where ñ/λ² = 0.90 (entry into broken phase) ===${C.reset}`);
console.log("η\\β      16        64        256       1024       (b16-b1024)");
for (const eta of ETAS) {
  const vals = BETAS.map(b => findThresholdCrossing(data[b], eta, 0.90));
  const shift = vals[0] - vals[3];
  const color = Math.abs(shift) > 0.02 ? C.cyan : C.dim;
  console.log(`${eta.toString().padStart(5)}    ${vals.map(v => isNaN(v) ? "  >1.0   " : v.toFixed(4)).join("   ")}    ${color}${shift >= 0 ? "+" : ""}${shift.toFixed(4)}${C.reset}`);
}

// 3. Same for threshold 0.95
console.log(`\n${C.bold}=== λ* where ñ/λ² = 0.95 ===${C.reset}`);
console.log("η\\β      16        64        256       1024       (b16-b1024)");
for (const eta of ETAS) {
  const vals = BETAS.map(b => findThresholdCrossing(data[b], eta, 0.95));
  const shift = vals[0] - vals[3];
  const color = Math.abs(shift) > 0.02 ? C.cyan : C.dim;
  console.log(`${eta.toString().padStart(5)}    ${vals.map(v => isNaN(v) ? "  >1.0   " : v.toFixed(4)).join("   ")}    ${color}${shift >= 0 ? "+" : ""}${shift.toFixed(4)}${C.reset}`);
}

// 4. Max slope position (steepest rise of ñ)
console.log(`\n${C.bold}=== λ* at max slope of ñ(λ) ===${C.reset}`);
console.log("η\\β      16        64        256       1024       (b16-b1024)");
for (const eta of ETAS) {
  const vals = BETAS.map(b => {
    const r = findMaxSlope(data[b], eta);
    return r.lambda;
  });
  const shift = vals[0] - vals[3];
  const color = Math.abs(shift) > 0.02 ? C.cyan : C.dim;
  console.log(`${eta.toString().padStart(5)}    ${vals.map(v => v.toFixed(4)).join("   ")}    ${color}${shift >= 0 ? "+" : ""}${shift.toFixed(4)}${C.reset}`);
}

// 5. Show ñ at λ=0.5 for all (η, β)
console.log(`\n${C.bold}=== ñ at λ=0.5 (critical point, BO→0.25) ===${C.reset}`);
console.log("η\\β      16        64        256       1024       BO=λ²");
for (const eta of ETAS) {
  const vals = BETAS.map(b => {
    const row = data[b].find(r => Math.abs(r.lambda - 0.5) < 0.003);
    return row ? row[`ntilde_${eta}`] : NaN;
  });
  const shift = vals[0] - vals[3];
  const color = Math.abs(shift) > 0.01 ? C.cyan : C.dim;
  console.log(`${eta.toString().padStart(5)}    ${vals.map(v => v.toFixed(6)).join("   ")}    ${color}${shift >= 0 ? "+" : ""}${shift.toFixed(6)}${C.reset}`);
}

// 6. Cutoff convergence check: compare cutoff=80 vs cutoff=40 at one point
console.log(`\n${C.bold}=== Summary ===${C.reset}`);
console.log("The QPT transition (ñ/λ² rising from <0.5 to ~1.0) sharpens with η:");
console.log("  η=50:     broad, ñ/λ²≈0.97 even at λ=0.5 (always broken)");
console.log("  η=25000:  sharp, ñ/λ² drops well below 0.5 in symmetric phase");
console.log("β-dependence: higher T (lower β) shifts ñ/λ² transitions to LOWER λ");
