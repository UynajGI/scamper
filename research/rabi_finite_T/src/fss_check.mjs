// Standalone Rabi model FSS analysis
// Build H from scratch, verify against known results, compute U4(x) and U4(n)
//
// H = (Delta/2) sigma_z + Omega * a^dag a + g * sigma_x * (a + a^dag)
//
// Basis: |n, s> where n=0..cutoff-1, s=0(up)/1(down)
// Index: state = 2*n + s  (s is fast variable)
//
// sigma_z|up> = +1, sigma_z|down> = -1
// sigma_x|up> = |down>, sigma_x|down> = |up>
// a|n> = sqrt(n)|n-1>, a^dag|n> = sqrt(n+1)|n+1>

const DELTA = 1.0;

// Build Hamiltonian matrix
function buildH(omega, g, cutoff) {
    const dim = 2 * cutoff;
    const H = Array(dim).fill(0).map(() => Array(dim).fill(0));

    for (let n = 0; n < cutoff; n++) {
        // Diagonal: Omega*n + (Delta/2)*sigma_z
        H[2*n][2*n] = omega * n + DELTA/2;       // |n,up>: sigma_z = +1
        H[2*n+1][2*n+1] = omega * n - DELTA/2;    // |n,down>: sigma_z = -1

        // Off-diagonal: g * sigma_x * (a + a^dag)
        // sigma_x flips spin, (a+a^dag) changes n by +/-1
        // Matrix element: g * <n',s'| sigma_x * (a+a^dag) |n,s>
        //              = g * <s'|sigma_x|s> * <n'|(a+a^dag)|n>
        // sigma_x: <up|sigma_x|down> = 1, <down|sigma_x|up> = 1
        // (a+a^dag): <n+1|(a+a^dag)|n> = sqrt(n+1), <n-1|(a+a^dag)|n> = sqrt(n)

        // |n,up> -> |n+1,down>: g * sqrt(n+1)
        if (n + 1 < cutoff) {
            const me = g * Math.sqrt(n + 1);
            H[2*(n+1)+1][2*n] = me;
            H[2*n][2*(n+1)+1] = me;
        }
        // |n,down> -> |n+1,up>: g * sqrt(n+1)
        if (n + 1 < cutoff) {
            const me = g * Math.sqrt(n + 1);
            H[2*(n+1)][2*n+1] = me;
            H[2*n+1][2*(n+1)] = me;
        }
    }

    return H;
}

// Build (a + a^dag)^power operator
function buildDisplacement(cutoff, power) {
    const dim = 2 * cutoff;

    // Build (a+a^dag) in boson space
    const x = Array(cutoff).fill(0).map(() => Array(cutoff).fill(0));
    for (let n = 0; n < cutoff - 1; n++) {
        const sq = Math.sqrt(n + 1);
        x[n][n+1] = sq;
        x[n+1][n] = sq;
    }

    // x^power
    let xk = Array(cutoff).fill(0).map((_, i) =>
        Array(cutoff).fill(0).map((_, j) => i === j ? 1 : 0));
    for (let p = 0; p < power; p++) {
        const result = Array(cutoff).fill(0).map(() => Array(cutoff).fill(0));
        for (let i = 0; i < cutoff; i++)
            for (let k = 0; k < cutoff; k++) {
                if (xk[i][k] === 0) continue;
                for (let j = 0; j < cutoff; j++)
                    result[i][j] += xk[i][k] * x[k][j];
            }
        xk = result;
    }

    // Embed in full space: I_spin tensor x^power
    const op = Array(dim).fill(0).map(() => Array(dim).fill(0));
    for (let i = 0; i < cutoff; i++)
        for (let j = 0; j < cutoff; j++)
            for (let s = 0; s < 2; s++)
                op[2*i+s][2*j+s] = xk[i][j];

    return op;
}

// Build n = a^dag a operator (photon number)
function buildPhotonNumber(cutoff) {
    const dim = 2 * cutoff;
    const op = Array(dim).fill(0).map(() => Array(dim).fill(0));
    for (let n = 0; n < cutoff; n++)
        for (let s = 0; s < 2; s++)
            op[2*n+s][2*n+s] = n;
    return op;
}

// Matrix-vector multiply
function matVec(A, v) {
    const n = A.length;
    return v.map((_, i) => {
        let s = 0;
        for (let j = 0; j < n; j++) s += A[i][j] * v[j];
        return s;
    });
}

// Power iteration to find ground state (for small matrices, use full diagonalization)
// For simplicity, use Jacobi rotation
function jacobi(A, maxIter = 200) {
    const n = A.length;
    const M = A.map(r => [...r]);
    const V = Array(n).fill(0).map((_, i) =>
        Array(n).fill(0).map((_, j) => i === j ? 1 : 0));

    for (let iter = 0; iter < maxIter * n * n; iter++) {
        // Find largest off-diagonal
        let p = 0, q = 1, largest = 0;
        for (let i = 0; i < n; i++)
            for (let j = i + 1; j < n; j++)
                if (Math.abs(M[i][j]) > largest) {
                    largest = Math.abs(M[i][j]);
                    p = i; q = j;
                }
        if (largest < 1e-14) break;

        const theta = 0.5 * Math.atan2(2 * M[p][q], M[q][q] - M[p][p]);
        const c = Math.cos(theta), s = Math.sin(theta);

        for (let i = 0; i < n; i++) {
            if (i !== p && i !== q) {
                const ip = M[i][p], iq = M[i][q];
                M[i][p] = c * ip - s * iq;
                M[p][i] = M[i][p];
                M[i][q] = s * ip + c * iq;
                M[q][i] = M[i][q];
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

    // Sort eigenvalues
    const idx = Array(n).fill(0).map((_, i) => i);
    idx.sort((a, b) => M[a][a] - M[b][b]);
    const eigenvalues = idx.map(i => M[i][i]);
    const eigenvectors = Array(n).fill(0).map((_, row) =>
        idx.map(i => V[row][i]));

    return { eigenvalues, eigenvectors };
}

// Thermal expectation <O>_beta
function thermalExp(eigen, op, beta) {
    const dim = eigen.eigenvalues.length;
    const ground = eigen.eigenvalues[0];
    let Z = 0, result = 0;

    for (let k = 0; k < dim; k++) {
        const w = Math.exp(-beta * (eigen.eigenvalues[k] - ground));
        if (w < 1e-300) continue;
        // <k|O|k> = sum_{i,j} V[i][k] * O[i][j] * V[j][k]
        let diag = 0;
        for (let i = 0; i < dim; i++) {
            const vik = eigen.eigenvectors[i][k];
            if (vik === 0) continue;
            for (let j = 0; j < dim; j++)
                diag += vik * op[i][j] * eigen.eigenvectors[j][k];
        }
        result += w * diag;
        Z += w;
    }
    return result / Z;
}

// ─── Sanity check ─────────────────────────────────────────────────────

console.log("=== Sanity check: g=0 (uncoupled) ===");
{
    const omega = 5.0, g = 0.0, cutoff = 30;
    const H = buildH(omega, g, cutoff);
    const eigen = jacobi(H);

    // Ground state: |0, down>, E = -Delta/2
    console.log("E0 =", eigen.eigenvalues[0].toFixed(6), "(expected", (-DELTA/2).toFixed(6), ")");

    const x2 = buildDisplacement(cutoff, 2);
    const n_op = buildPhotonNumber(cutoff);
    const beta = 100;
    console.log("<x^2> =", thermalExp(eigen, x2, beta).toFixed(6), "(expected 1.0)");
    console.log("<n> =", thermalExp(eigen, n_op, beta).toFixed(6), "(expected 0.0)");
    console.log("U4(x) =", (1 - thermalExp(eigen, buildDisplacement(cutoff, 4), beta) / (3 * Math.pow(thermalExp(eigen, x2, beta), 2))).toFixed(6), "(expected ~0, Gaussian)");
}

// ─── FSS analysis ─────────────────────────────────────────────────────
console.log("\n=== FSS: U4(x) and U4(n) vs r for different eta ===\n");

const omegas = [50.0, 100.0, 200.0];
const beta = 1024.0;
const cutoff = 50;
const r_values = [];
for (let i = 0; i < 30; i++) r_values.push(0.3 + 4.7 * i / 29);

console.log("r      g_c=3.54  U4x(50)  U4x(100) U4x(200)  U4n(50)  U4n(100) U4n(200)");

for (const r of r_values) {
    let row = r.toFixed(3).padStart(6) + "  ";
    const u4x = [], u4n = [];

    for (const omega of omegas) {
        const gc = Math.sqrt(omega * DELTA) / 2;
        const g = r * gc;
        const H = buildH(omega, g, cutoff);
        const eigen = jacobi(H);

        const x2 = thermalExp(eigen, buildDisplacement(cutoff, 2), beta);
        const x4 = thermalExp(eigen, buildDisplacement(cutoff, 4), beta);
        const u4_x = x2 > 1e-14 ? 1 - x4 / (3 * x2 * x2) : 0;

        const n2 = thermalExp(eigen, buildPhotonNumber(cutoff), beta);
        // n^2 = a^dag a a^dag a = a^dag(a^dag a + 1)a = (a^dag)^2 a^2 + a^dag a
        // Simpler: build n^2 and n^4 directly
        // n^2 op: n*n in diagonal
        const dim = 2 * cutoff;
        const n2_op = Array(dim).fill(0).map(() => Array(dim).fill(0));
        const n4_op = Array(dim).fill(0).map(() => Array(dim).fill(0));
        for (let n = 0; n < cutoff; n++) {
            for (let s = 0; s < 2; s++) {
                n2_op[2*n+s][2*n+s] = n*n;
                n4_op[2*n+s][2*n+s] = n*n*n*n;
            }
        }
        const nn2 = thermalExp(eigen, n2_op, beta);
        const nn4 = thermalExp(eigen, n4_op, beta);
        const u4_n = nn2 > 1e-14 ? 1 - nn4 / (3 * nn2 * nn2) : 0;

        u4x.push(u4_x);
        u4n.push(u4_n);
    }

    // Check for crossings
    let cross_x = "";
    if (r > 0.5) {
        if ((u4x[0] - u4x[1]) * (u4x[0] - u4x[1]) < 0.001) cross_x = " <-close50-100";
    }

    console.log(
        r.toFixed(3).padStart(6) + "  " +
        u4x.map(v => v.toFixed(6).padStart(9)).join(" ") + "  " +
        u4n.map(v => v.toFixed(6).padStart(9)).join(" ") + cross_x
    );
}
