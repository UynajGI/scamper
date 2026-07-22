// Compare U4 for +g vs -g to test sign independence

const DELTA = 1.0;

function buildH(omega, g, cutoff) {
    const dim = 2 * cutoff;
    const H = Array(dim).fill(0).map(() => Array(dim).fill(0));
    for (let n = 0; n < cutoff; n++) {
        H[2*n][2*n] = omega * n - DELTA/2;     // Down
        H[2*n+1][2*n+1] = omega * n + DELTA/2;  // Up
        if (n + 1 < cutoff) {
            const me = g * Math.sqrt(n + 1);
            H[2*n][2*(n+1)+1] = me;
            H[2*(n+1)+1][2*n] = me;
            H[2*n+1][2*(n+1)] = me;
            H[2*(n+1)][2*n+1] = me;
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
                if (Math.abs(M[i][j]) > largest) {
                    largest = Math.abs(M[i][j]); p = i; q = j;
                }
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
    return {
        eigenvalues: idx.map(i => M[i][i]),
        eigenvectors: Array(n).fill(0).map((_, row) => idx.map(i => V[row][i]))
    };
}

function buildDisp(cutoff, power) {
    const dim = 2 * cutoff;
    const x = Array(cutoff).fill(0).map(() => Array(cutoff).fill(0));
    for (let n = 0; n < cutoff - 1; n++) {
        const sq = Math.sqrt(n + 1);
        x[n][n+1] = sq; x[n+1][n] = sq;
    }
    let xk = Array(cutoff).fill(0).map((_, i) =>
        Array(cutoff).fill(0).map((_, j) => i === j ? 1 : 0));
    for (let p = 0; p < power; p++) {
        const r = Array(cutoff).fill(0).map(() => Array(cutoff).fill(0));
        for (let i = 0; i < cutoff; i++)
            for (let k = 0; k < cutoff; k++) {
                if (xk[i][k] === 0) continue;
                for (let j = 0; j < cutoff; j++)
                    r[i][j] += xk[i][k] * x[k][j];
            }
        xk = r;
    }
    const op = Array(dim).fill(0).map(() => Array(dim).fill(0));
    for (let i = 0; i < cutoff; i++)
        for (let j = 0; j < cutoff; j++)
            for (let s = 0; s < 2; s++)
                op[2*i+s][2*j+s] = xk[i][j];
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
        res += w * d;
        Z += w;
    }
    return res / Z;
}

// Compare +g vs -g
const omega = 50.0, g_abs = 5.0, cutoff = 50, beta = 1024;

for (const sign of [+1, -1]) {
    const g = sign * g_abs;
    const H = buildH(omega, g, cutoff);
    const eigen = jacobi(H);
    const x2 = thermalExp(eigen, buildDisp(cutoff, 2), beta);
    const x4 = thermalExp(eigen, buildDisp(cutoff, 4), beta);
    const u4 = x2 > 1e-14 ? 1 - x4 / (3 * x2 * x2) : 0;
    console.log(`g=${g.toFixed(2)}: E0=${eigen.eigenvalues[0].toFixed(6)} <x2>=${x2.toFixed(8)} <x4>=${x4.toFixed(8)} U4=${u4.toFixed(8)}`);
}
