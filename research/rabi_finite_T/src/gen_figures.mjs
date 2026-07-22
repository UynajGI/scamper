#!/usr/bin/env node
// Generate figure data for the Rabi QPT report.
// Output: report/figs/*.csv

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

const DIR = "/home/jiangyuan/scuttle/research/rabi_finite_T";
const OUT = path.join(DIR, "report/figs");

// ===== Figure 1: QMC vs ED validation (n_tilde) =====
{
  const mcCsv = fs.readFileSync(path.join(DIR, "results/qmc_verify/verify.csv"), "utf8");
  const edCsv = fs.readFileSync(path.join(DIR, "results/ed_fss/beta_1024.csv"), "utf8");
  const mcData = mcCsv.trim().split("\n").slice(1).map(l => l.split(",").map(parseFloat));
  const edData = edCsv.trim().split("\n").slice(1).map(l => l.split(",").map(parseFloat));
  
  let out = "lambda,ntilde_qmc,ntilde_ed,eta\n";
  for (const r of mcData) {
    const lam = r[0];
    if (lam < 0.2 || lam > 2.0) continue;
    // QMC: ntilde at col 3, 6, 9, 12 for eta 10,50,200,500
    // ED: ntilde at col 2, 5, 8, 11 for eta 10,50,200,500
    const etas = [10, 50, 200, 500];
    for (let i = 0; i < 4; i++) {
      const edRow = edData.find(e => Math.abs(e[0] - lam) < 0.005);
      if (edRow) {
        out += `${lam},${r[3 + i*3]},${edRow[2 + i*3]},${etas[i]}\n`;
      }
    }
  }
  fs.writeFileSync(path.join(OUT, "qmc_vs_ed_ntilde.csv"), out);
}

// ===== Figure 2: U4(x) vs lambda =====
{
  const csv = fs.readFileSync(path.join(DIR, "results/qmc_verify/verify_u4.csv"), "utf8");
  const data = csv.trim().split("\n").slice(1).map(l => l.split(",").map(parseFloat));
  let out = "lambda,u4_qmc,u4_ed,eta\n";
  for (const r of data) {
    const lam = r[0];
    if (lam > 2.0) continue;
    for (let i = 0; i < 4; i++) {
      const eta = [10, 50, 200, 500][i];
      out += `${lam},${r[1+i*4]},${r[2+i*4]},${eta}\n`;
    }
  }
  fs.writeFileSync(path.join(OUT, "u4x_vs_lambda.csv"), out);
}

// ===== Figure 3: U4(x) extended (ED only, larger λ range) =====
{
  const csv = fs.readFileSync(path.join(DIR, "results/ed_fss/beta_1024.csv"), "utf8");
  const data = csv.trim().split("\n").slice(1).map(l => l.split(",").map(parseFloat));
  // Columns: lambda, ntilde_eta50, sigmaz, x2, u4x, gap, cv, ...(rep for each eta)
  let out = "lambda,u4x,eta\n";
  for (const r of data) {
    const lam = r[0];
    for (let i = 0; i < 5; i++) {
      const eta = [50, 200, 1000, 5000, 25000][i];
      const u4 = r[1 + i*6 + 3];
      if (lam <= 5) out += `${lam},${u4},${eta}\n`;
    }
  }
  fs.writeFileSync(path.join(OUT, "u4x_ed_extended.csv"), out);
}

// ===== Figure 4: Energy gap vs lambda =====
{
  const csv = fs.readFileSync(path.join(DIR, "results/ed_fss/beta_1024.csv"), "utf8");
  const data = csv.trim().split("\n").slice(1).map(l => l.split(",").map(parseFloat));
  let out = "lambda,gap,eta\n";
  for (const r of data) {
    const lam = r[0];
    for (let i = 0; i < 5; i++) {
      const eta = [50, 200, 1000, 5000, 25000][i];
      out += `${lam},${r[1 + i*6 + 4]},${eta}\n`;
    }
  }
  fs.writeFileSync(path.join(OUT, "gap_vs_lambda.csv"), out);
}

// ===== Figure 5: C_V vs lambda for different beta (η=50) =====
{
  let out = "lambda,cv,beta\n";
  for (const beta of [4, 16, 64, 256, 1024]) {
    try {
      const csv = fs.readFileSync(path.join(DIR, `results/ed_fss/beta_${beta}.csv`), "utf8");
      const data = csv.trim().split("\n").slice(1).map(l => l.split(",").map(parseFloat));
      for (const r of data) {
        if (r[0] <= 20) out += `${r[0]},${r[1 + 0*6 + 5]},${beta}\n`;
      }
    } catch(e) {}
  }
  fs.writeFileSync(path.join(OUT, "cv_vs_lambda_eta50.csv"), out);
}

// ===== Figure 6: C_V peak position scaling =====
{
  let out = "predicted,observed,eta,beta\n";
  const peaks = {
    // eta=50: beta -> lambda*
    50: {4:3.57, 16:6.89, 64:9.06, 256:10.81, 1024:12.31},
    200: {4:7.15, 16:13.77, 64:18.12, 256:21.61, 1024:24.61},
    1000: {4:15.99, 16:30.80, 64:40.52, 256:48.33},
  };
  for (const [eta, betas] of Object.entries(peaks)) {
    for (const [beta, lam] of Object.entries(betas)) {
      const predicted = (eta/2) * Math.log(beta/2.4);
      out += `${predicted},${lam*lam},${eta},${beta}\n`;
    }
  }
  fs.writeFileSync(path.join(OUT, "cv_peak_scaling.csv"), out);
}

// ===== Figure 7: U4(n) vs lambda =====
{
  const csv = fs.readFileSync(path.join(DIR, "results/qmc_verify/u4n.csv"), "utf8");
  const data = csv.trim().split("\n").slice(1).map(l => l.split(",").map(parseFloat));
  let out = "lambda,u4n_ed,eta\n";
  for (const r of data) {
    for (let i = 0; i < 4; i++) {
      const eta = [10, 50, 200, 500][i];
      out += `${r[0]},${r[2+i*3]},${eta}\n`;
    }
  }
  fs.writeFileSync(path.join(OUT, "u4n_vs_lambda.csv"), out);
}

// ===== Figure 8: Wormhole U4 vs lambda =====
{
  let out_w = "lambda,u4,eta\n";
  for (const beta of [1024]) {
    try {
      const csv = fs.readFileSync(path.join(DIR, `results/qmc_wormhole/beta_${beta}.csv`), "utf8");
      const data = csv.trim().split("\n").slice(1).map(l => l.split(",").map(parseFloat));
      for (const r of data) {
        out_w += `${r[0]},${r[4]},50\n`;
        out_w += `${r[0]},${r[10]},200\n`;
        out_w += `${r[0]},${r[16]},1000\n`;
      }
    } catch(e) {}
  }
  fs.writeFileSync(path.join(OUT, "wormhole_u4.csv"), out_w);
}

// ===== Figure 9: ntilde over lambda^2 vs lambda =====
{
  const csv = fs.readFileSync(path.join(DIR, "results/ed_fss/beta_1024.csv"), "utf8");
  const data = csv.trim().split("\n").slice(1).map(l => l.split(",").map(parseFloat));
  let out = "lambda,ntilde_over_lam2,eta\n";
  for (const r of data) {
    const lam = r[0];
    if (lam > 5) continue;
    for (let i = 0; i < 5; i++) {
      const eta = [50, 200, 1000, 5000, 25000][i];
      out += `${lam},${r[1 + i*6 + 2]},${eta}\n`;
    }
  }
  fs.writeFileSync(path.join(OUT, "ntilde_over_lam2.csv"), out);
}

console.log("Generated figure data in", OUT);
