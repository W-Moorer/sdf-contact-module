import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


def solve_mass(M, x):
    if isinstance(M, csr_matrix):
        return spsolve(M, x)
    return np.linalg.solve(M, x)


def solve_rank_deficient_kkt(M, C, Q, J, b_c, rank_tol=1e-9):
    """
    Solve:
        M a - J^T λ = Q - C   (dynamics)
        J a = b_c              (acceleration-level constraints)

    Handles rank-deficient J using scaled Schur complement + SVD pseudo-inverse.

    Returns: (acc, lambda, diagnostics)
    """
    nc = J.shape[0]
    nv = M.shape[0]

    if nv == 0:
        return np.zeros(0), b_c.copy(), {'rank': 0, 'num_constraints': nc, 'redundant': nc}

    if nc == 0:
        rhs = Q - C
        a = solve_mass(M, rhs)
        return a, np.zeros(0), {'rank': 0, 'num_constraints': 0, 'redundant': 0}

    J_dense = J.toarray() if isinstance(J, csr_matrix) else np.asarray(J)
    M_dense = M.toarray() if isinstance(M, csr_matrix) else np.asarray(M)

    # 1. MinvQ = M^{-1} (Q - C)
    MinvQ = solve_mass(M, Q - C)

    # 2. Row scaling using approximate row norm
    W_raw = J_dense @ solve_mass(M, J_dense.T)
    row_norm = np.sqrt(np.maximum(np.abs(np.diag(W_raw)), 1e-30))
    D = 1.0 / np.maximum(row_norm, 1e-12)

    Js = D[:, None] * J_dense
    bs = D * b_c

    # 3. Scaled Schur complement
    Minv_JsT = solve_mass(M, Js.T)
    W = Js @ Minv_JsT
    rhs = bs - Js @ MinvQ

    # 4. SVD pseudo-inverse
    U, S, Vt = np.linalg.svd(W, full_matrices=False)

    smax = S[0] if len(S) > 0 else 0.0
    keep = S > rank_tol * max(smax, 1.0)

    # 5. Consistency check
    rhs_in_range = U[:, keep] @ (U[:, keep].T @ rhs)
    inconsistency = np.linalg.norm(rhs - rhs_in_range)

    # 6. Solve λ = W⁺ rhs
    alpha = Vt.T[:, keep] @ ((U[:, keep].T @ rhs) / S[keep])

    # 7. Acceleration
    a = MinvQ + Minv_JsT @ alpha

    # 8. Diagnostics
    constraint_residual = np.linalg.norm(J_dense @ a - b_c)
    rank = int(np.sum(keep))

    diag = {
        'rank': rank,
        'num_constraints': nc,
        'redundant': nc - rank,
        'constraint_residual': float(constraint_residual),
        'smallest_singular_value': float(S[-1]) if len(S) else 0.0,
        'inconsistency': float(inconsistency),
    }

    return a, alpha, diag


def solve_kkt(M, C, Q, J, b_c):
    """Wrapper: calls solve_rank_deficient_kkt and extracts (acc, lam)."""
    a, lam, diag = solve_rank_deficient_kkt(M, C, Q, J, b_c)
    if diag.get('redundant', 0) > 0 and diag.get('inconsistency', 0) > 1e-6:
        print(f"  [KKT] {diag['redundant']} redundant constraints, "
              f"inconsistency={diag['inconsistency']:.2e}")
    return a, lam
