// Find the correct FSS observable for the Rabi QPT.
// Try: U4(n), <n>/η, susceptibility d<n>/dg
// At fixed r = g/g_c, scan many η values to see which observable crosses.

const DELTA = 1.0;

function buildH(omega, g, cutoff) {
    const dim = 2 * cutoff;
    const H = Array(dim).fill(0).map(() => Array(dim).fill(0));
    for (let n = 0; n < cutoff; n++) {
        H[2*n][2*n] = omega * n - DELTA/2;
        H[2*n+1][2*n+1] = omega * n + DELTA/2;
        if (n + 1 < cutoff) {
            const me = g * Math.sqrt(n + 1);
            H[2*n][2*(n+1)+1] = me; H[2*(n+1)+1][2*n] = me;
            H[2*n+1][2*(n+1)] = me; H[2*(n+1)][2*n+1] = me;
        }
    }
    return H;
}

function jacobi(A) {
    const n = A.length;
    const M = A.map(r => [...r]);
    const V = Array(n).fill(0).map((_, i) =>
        Array(n).fill(0).map((_, j) => i === j ? 1 : 0));
    for (let iter = 0; iter < 500 * n * n; iter++) {
        let p = 0, q = 1, largest = 0;
        for (let i = 0; i < n; i++)
            for (let j = i + 1; j < n; j++)
                if (Math.abs(M[i][j]) > largest) { largest = Math.abs(M[i][j]); p = i; q = j; }
        if (largest < 1e-15) break;
        const theta = 0.5 * Math.atan2(2 * M[p][q], M[q][q] - M[p][p]);
        const c = Math.cos(theta), s = Math.sin(theta);
        for (let i = 0; i < n; i++) {
            if (i !== p && i !== q) {
                const ip = M[i][p], iq = M[i][q];
                M[i][p] = c * ip - s * iq; M[p][i] = M[i][p];
                M[i][q] = s * ip + c * iq; M[q][i] = M[i][q];
            }
        }
        const pp = M[p][p], qq = M[q][q], pq = M[p][q];
        M[p][p] = c*c*pp - 2*s*c*pq + s*s*qq;
        M[q][q] = s*s*pp + 2*s*c*pq + c*c*qq;
        M[p][q] = 0; M[q][p] = 0;
        for (let i = 0; i < n; i++) {
            const vp = V[i][p], vq = V[i][q];
            V[i][p] = c * vp - s * vq;
            V[i][q] = s * vp + c * vq;
        }
    }
    const idx = Array(n).fill(0).map((_, i) => i);
    idx.sort((a, b) => M[a][a] - M[b][b]);
    return { eigenvalues: idx.map(i => M[i][i]), eigenvectors: V };
}

function buildOp(cutoff, kind) {
    const dim = 2 * cutoff;
    const op = Array(dim).fill(0).map(() => Array(dim).fill(0));
    if (kind === 'n') {
        for (let n = 0; n < cutoff; n++)
            for (let s = 0; s < 2; s++)
                op[2*n+s][2*n+s] = n;
    } else if (kind === 'n2') {
        for (let n = 0; n < cutoff; n++)
            for (let s = 0; s < 2; s++)
                op[2*n+s][2*n+s] = n * n;
    } else if (kind === 'n4') {
        for (let n = 0; n < cutoff; n++)
            for (let s = 0; s < 2; s++)
                op[2*n+s][2*n+s] = n*n*n*n;
    } else if (kind === 'x2' || kind === 'x4') {
        const power = kind === 'x2' ? 2 : 4;
        const x = Array(cutoff).fill(0).map(() => Array(cutoff).fill(0));
        for (let n = 0; n < cutoff - 1; n++) { const sq = Math.sqrt(n+1); x[n][n+1] = sq; x[n+1][n] = sq; }
        let xk = Array(cutoff).fill(0).map((_, i) => Array(cutoff).fill(0).map((_, j) => i===j?1:0));
        for (let p = 0; p < power; p++) {
            const r = Array(cutoff).fill(0).map(() => Array(cutoff).fill(0));
            for (let i = 0; i < cutoff; i++) for (let k = 0; k < cutoff; k++) {
                if (xk[i][k] === 0) continue;
                for (let j = 0; j < cutoff; j++) r[i][j] += xk[i][k] * x[k][j];
            }
            xk = r;
        }
        for (let i = 0; i < cutoff; i++) for (let j = 0; j < cutoff; j++)
            for (let s = 0; s < 2; s++) op[2*i+s][2*j+s] = xk[i][j];
    }
    return op;
}

function thermalExp(eigen, op, beta) {
    const dim = eigen.eigenvalues.length;
    const ground = eigen.eigenvalues[0];
    let Z = 0, res = 0;
    for (let k = 0; k < dim; k++) {
        const w = Math.exp(-beta * (eigen.eigenvalues[k] - ground));
        if (w < 1e-300) continue;
        let d = 0;
        for (let i = 0; i < dim; i++) {
            const vik = eigen.eigenvectors[i][k];
            if (vik === 0) continue;
            for (let j = 0; j < dim; j++)
                d += vik * op[i][j] * eigen.eigenvectors[j][k];
        }
        res += w * d; Z += w;
    }
    return res / Z;
}

// ─── Main: scan multiple observables at fixed r ──────────────────────
console.log("=== Observable comparison at fixed r ===\n");
console.log("Looking for which observable shows η-crossings\n");

const etas = [100, 1000, 10000];
const r_values = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0];
const beta = 1024.0;

console.log("r     η      <n>      n/η     <n²>    U4(n)   <x²>    U4(x)   <n>/η   χ_n");
for (const r of r_values) {
    let prev_n = {};
    for (const eta of etas) {
        const omega = eta * DELTA;
        const gc = Math.sqrt(omega * DELTA) / 2;
        const g = r * gc;
        // cutoff: need ~10*<n>_max. <n> ~ r²Δ/(4Ω) = r²/(4η). Even r=10 gives <n>=0.025 at η=1000.
        // But need larger cutoff for η=100 (⟨n⟩ up to 0.25 at r=5).
        const cutoff = Math.max(30, Math.ceil(20 * r * r / eta) + 10);
        const H = buildH(omega, g, cutoff);
        const eigen = jacobi(H);
        const beta_val = beta;

        const n_mean = thermalExp(eigen, buildOp(cutoff, 'n'), beta_val);
        const n2 = thermalExp(eigen, buildOp(cutoff, 'n2'), beta_val);
        const n4 = thermalExp(eigen, buildOp(cutoff, 'n4'), beta_val);
        const x2 = thermalExp(eigen, buildOp(cutoff, 'x2'), beta_val);
        const x4 = thermalExp(eigen, buildOp(cutoff, 'x4'), beta_val);

        const u4_n = n2 > 1e-14 ? 1 - n4 / (3 * n2 * n2) : 0;
        const u4_x = x2 > 1e-14 ? 1 - x4 / (3 * x2 * x2) : 0;
        const n_over_eta = n_mean / eta;

        console.log(
            `${r.toFixed(2).padStart(5)} ${eta.toString().padStart(6)} ` +
            `${n_mean.toFixed(8).padStart(10)} ${n_over_eta.toFixed(8).padStart(10)} ` +
            `${n2.toFixed(8).padStart(10)} ${u4_n.toFixed(6).padStart(8)} ` +
            `${x2.toFixed(8).padStart(10)} ${u4_x.toFixed(6).padStart(8)}`
        );
    }
    console.log("");
}

// Now check: does <n>/η cross at fixed r?
console.log("\n=== <n>/η at fixed r for many η ===");
console.log("If <n>/η → 0 below gc and → finite above gc, curves should cross near r=1");
console.log("\nr     <n>/η(100)  <n>/η(1000)  <n>/η(10000)");
for (const r of [0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0]) {
    const vals = [];
    for (const eta of [100, 1000, 10000]) {
        const omega = eta * DELTA;
        const gc = Math.sqrt(omega * DELTA) / 2;
        const g = r * gc;
        const cutoff = Math.max(30, Math.ceil(20 * r * r / eta) + 10);
        const H = buildH(omega, g, cutoff);
        const eigen = jacobi(H);
        const n_mean = thermalExp(eigen, buildOp(cutoff, 'n'), 1024.0);
        vals.push(n_mean / eta);
    }
    console.log(
        `${r.toFixed(2).padStart(5)} ` +
        vals.map(v => v.toExponential(3).padStart(13)).join(" ")
    );
}
