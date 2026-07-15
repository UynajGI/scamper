//! Exact finite-basis transfer matrix and local estimators.

use crate::impurity::ImpurityError;

#[derive(Debug, Clone)]
pub struct SymmetricEigensystem {
    pub values: Vec<f64>,
    /// Eigenvectors stored by rows: `vectors[basis][eigenstate]`.
    pub vectors: Vec<Vec<f64>>,
}

impl SymmetricEigensystem {
    pub fn diagonalize(matrix: Vec<Vec<f64>>) -> Result<Self, ImpurityError> {
        let dimension = matrix.len();
        if dimension == 0 || matrix.iter().any(|row| row.len() != dimension) {
            return Err(ImpurityError::InvalidConfiguration(
                "Hamiltonian must be non-empty and square".into(),
            ));
        }
        let (values, vectors) = jacobi_eigensystem(matrix);
        Ok(Self { values, vectors })
    }

    pub fn matrix_function(&self, function: impl Fn(f64) -> f64) -> Vec<Vec<f64>> {
        let n = self.values.len();
        let mut result = vec![vec![0.0; n]; n];
        for (i, row) in result.iter_mut().enumerate() {
            for (j, element) in row.iter_mut().enumerate() {
                *element = (0..n)
                    .map(|state| {
                        self.vectors[i][state]
                            * function(self.values[state])
                            * self.vectors[j][state]
                    })
                    .sum();
            }
        }
        result
    }

    pub fn thermal_density_matrix(&self, beta: f64) -> Vec<Vec<f64>> {
        let ground = self.values[0];
        let weights = self
            .values
            .iter()
            .map(|&value| (-beta * (value - ground)).exp())
            .collect::<Vec<_>>();
        let partition: f64 = weights.iter().sum();
        let n = self.values.len();
        let mut rho = vec![vec![0.0; n]; n];
        for (i, row) in rho.iter_mut().enumerate() {
            for (j, element) in row.iter_mut().enumerate() {
                *element = (0..n)
                    .map(|state| self.vectors[i][state] * weights[state] * self.vectors[j][state])
                    .sum::<f64>()
                    / partition;
            }
        }
        rho
    }
}

pub fn multiply(left: &[Vec<f64>], right: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = left.len();
    let mut result = vec![vec![0.0; n]; n];
    for i in 0..n {
        for k in 0..n {
            let value = left[i][k];
            if value == 0.0 {
                continue;
            }
            for j in 0..n {
                result[i][j] += value * right[k][j];
            }
        }
    }
    result
}

#[allow(clippy::needless_range_loop)]
fn jacobi_eigensystem(mut matrix: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = matrix.len();
    let mut vectors = vec![vec![0.0; n]; n];
    for i in 0..n {
        vectors[i][i] = 1.0;
    }
    for _ in 0..(160 * n * n) {
        let mut p = 0usize;
        let mut q = 1usize.min(n.saturating_sub(1));
        let mut largest = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                if matrix[i][j].abs() > largest {
                    largest = matrix[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if largest < 1e-13 {
            break;
        }
        let angle = 0.5 * (2.0 * matrix[p][q]).atan2(matrix[q][q] - matrix[p][p]);
        let (s, c) = angle.sin_cos();
        for i in 0..n {
            if i != p && i != q {
                let ip = matrix[i][p];
                let iq = matrix[i][q];
                matrix[i][p] = c * ip - s * iq;
                matrix[p][i] = matrix[i][p];
                matrix[i][q] = s * ip + c * iq;
                matrix[q][i] = matrix[i][q];
            }
        }
        let pp = matrix[p][p];
        let qq = matrix[q][q];
        let pq = matrix[p][q];
        matrix[p][p] = c * c * pp - 2.0 * s * c * pq + s * s * qq;
        matrix[q][q] = s * s * pp + 2.0 * s * c * pq + c * c * qq;
        matrix[p][q] = 0.0;
        matrix[q][p] = 0.0;
        for row in &mut vectors {
            let vp = row[p];
            let vq = row[q];
            row[p] = c * vp - s * vq;
            row[q] = s * vp + c * vq;
        }
    }
    let mut order = (0..n).collect::<Vec<_>>();
    order.sort_by(|&a, &b| matrix[a][a].total_cmp(&matrix[b][b]));
    let values = order.iter().map(|&i| matrix[i][i]).collect();
    let sorted = (0..n)
        .map(|row| order.iter().map(|&col| vectors[row][col]).collect())
        .collect();
    (values, sorted)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonalize_recovers_known_eigenvalues_of_diagonal_matrix() {
        // diag(3, -1, 2) has eigenvalues {-1, 2, 3}.
        let matrix = vec![
            vec![3.0, 0.0, 0.0],
            vec![0.0, -1.0, 0.0],
            vec![0.0, 0.0, 2.0],
        ];
        let eigen = SymmetricEigensystem::diagonalize(matrix).expect("diagonalize");
        assert!((eigen.values[0] - (-1.0)).abs() < 1e-10);
        assert!((eigen.values[1] - 2.0).abs() < 1e-10);
        assert!((eigen.values[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn diagonalize_recovers_eigenvalues_of_symmetric_2x2() {
        // [[2, 1], [1, 2]] has eigenvalues {1, 3}.
        let matrix = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
        let eigen = SymmetricEigensystem::diagonalize(matrix).expect("diagonalize");
        assert!((eigen.values[0] - 1.0).abs() < 1e-10);
        assert!((eigen.values[1] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn eigenvectors_are_orthonormal() {
        let matrix = vec![vec![2.0, 1.0], vec![1.0, 2.0]];
        let eigen = SymmetricEigensystem::diagonalize(matrix).expect("diagonalize");
        let n = eigen.values.len();
        for i in 0..n {
            for j in 0..n {
                let overlap: f64 = (0..n)
                    .map(|k| eigen.vectors[k][i] * eigen.vectors[k][j])
                    .sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (overlap - expected).abs() < 1e-10,
                    "overlap <{i},{j}> = {overlap}"
                );
            }
        }
    }

    #[test]
    fn rejects_empty_or_non_square_matrix() {
        assert!(SymmetricEigensystem::diagonalize(vec![]).is_err());
        let ragged = vec![vec![1.0, 0.0], vec![0.0]];
        assert!(SymmetricEigensystem::diagonalize(ragged).is_err());
    }

    #[test]
    fn matrix_function_applies_to_eigenvalues() {
        // For diag(1, 4): exp returns diag(exp(1), exp(4)).
        let matrix = vec![vec![1.0, 0.0], vec![0.0, 4.0]];
        let eigen = SymmetricEigensystem::diagonalize(matrix).expect("diagonalize");
        let result = eigen.matrix_function(f64::exp);
        assert!((result[0][0] - 1.0_f64.exp()).abs() < 1e-10);
        assert!((result[1][1] - 4.0_f64.exp()).abs() < 1e-10);
        assert!(result[0][1].abs() < 1e-10);
    }

    #[test]
    fn thermal_density_matrix_has_unit_trace_and_is_symmetric() {
        let matrix = vec![vec![0.0, 0.5], vec![0.5, 1.0]];
        let eigen = SymmetricEigensystem::diagonalize(matrix).expect("diagonalize");
        let rho = eigen.thermal_density_matrix(2.0);
        let trace = rho[0][0] + rho[1][1];
        assert!(
            (trace - 1.0).abs() < 1e-10,
            "trace should be 1, got {trace}"
        );
        assert!((rho[0][1] - rho[1][0]).abs() < 1e-10);
    }

    #[test]
    fn thermal_density_matrix_approaches_ground_state_at_large_beta() {
        // eigenvalues {~−0.3, ~1.3}; at large beta the density matrix projects
        // onto the ground eigenstate.
        let matrix = vec![vec![0.0, 0.5], vec![0.5, 1.0]];
        let eigen = SymmetricEigensystem::diagonalize(matrix).expect("diagonalize");
        let rho = eigen.thermal_density_matrix(40.0);
        let ground_weight = eigen.vectors.iter().map(|row| row[0] * row[0]).sum::<f64>();
        // Diagonal of rho in the original basis approaches |v_0><v_0|.
        let expected00 = eigen.vectors[0][0].powi(2) * ground_weight / ground_weight;
        assert!((rho[0][0] - expected00).abs() < 1e-6);
    }

    #[test]
    fn multiply_is_associative_with_identity() {
        let identity = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let matrix = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let product = multiply(&identity, &matrix);
        for i in 0..2 {
            for j in 0..2 {
                assert!((product[i][j] - matrix[i][j]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn multiply_matches_manual_computation() {
        let a = vec![vec![1.0, 2.0], vec![0.0, 1.0]];
        let b = vec![vec![3.0, 0.0], vec![1.0, 4.0]];
        let c = multiply(&a, &b);
        // Manual: [[1*3+2*1, 1*0+2*4],[0*3+1*1, 0*0+1*4]] = [[5,8],[1,4]]
        assert!((c[0][0] - 5.0).abs() < 1e-14);
        assert!((c[0][1] - 8.0).abs() < 1e-14);
        assert!((c[1][0] - 1.0).abs() < 1e-14);
        assert!((c[1][1] - 4.0).abs() < 1e-14);
    }
}
