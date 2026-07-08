import numpy as np


class AABB:
    def __init__(self, center, half_extents):
        self.center = np.asarray(center, dtype=float)
        self.half = np.asarray(half_extents, dtype=float)

    @classmethod
    def from_points(cls, points):
        if len(points) == 0:
            return cls([0, 0, 0], [0, 0, 0])
        lo = points.min(axis=0)
        hi = points.max(axis=0)
        return cls((lo + hi) / 2, (hi - lo) / 2)


class SDFContactEngine:
    def __init__(self):
        self._pairs = {}
        self._aabb_a = {}
        self._aabb_b = {}

    def register_pair(self, cp_id, quadrature_mesh, sdf_grid,
                      aabb_b_half=None, narrow_margin=0.01):
        self._pairs[cp_id] = (quadrature_mesh, sdf_grid)
        self._aabb_a[cp_id] = AABB.from_points(quadrature_mesh.X_q)
        self._aabb_b[cp_id] = aabb_b_half
        self._narrow_margin = narrow_margin

    def _check_broad_phase(self, cp, state, model):
        body_a = model.bodies[cp.body_a_id]
        body_b = model.bodies[cp.body_b_id]
        idx_a = state.body_idx(body_a.id)
        idx_b = state.body_idx(body_b.id)

        rA, RA = state.r[idx_a], state.R[idx_a]
        rB = state.r[idx_b]

        aabb_a_local = self._aabb_a.get(cp.id)
        if aabb_a_local is None:
            return True

        aabb_a_world_center = rA + RA @ aabb_a_local.center
        aabb_a_world_half = np.sum(np.abs(RA) * aabb_a_local.half[None, :], axis=1)

        aabb_b_half = self._aabb_b.get(cp.id)
        if aabb_b_half is None:
            return True

        diff = aabb_a_world_center - rB
        return all(abs(diff[i]) <= aabb_a_world_half[i] + aabb_b_half[i] for i in range(3))

    def _narrow_filter(self, X_w_local, rB, RB):
        margin = self._narrow_margin
        if RB is None:
            return slice(None)
        rB_flat = rB.ravel()
        dz = X_w_local[:, 2] - rB_flat[2]
        mask = np.abs(dz) < margin
        return mask

    def evaluate(self, cp, state, model):
        pair_key = cp.id
        if pair_key not in self._pairs:
            return np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

        if not self._check_broad_phase(cp, state, model):
            return np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

        quad_mesh, sdf = self._pairs[pair_key]

        body_a = model.bodies[cp.body_a_id]
        body_b = model.bodies[cp.body_b_id]
        idx_a = state.body_idx(body_a.id)
        idx_b = state.body_idx(body_b.id)

        rA, RA = state.r[idx_a], state.R[idx_a]
        rB, RB = state.r[idx_b], state.R[idx_b]
        vA, wA = state.v[idx_a], state.omega[idx_a]
        vB, wB = state.v[idx_b], state.omega[idx_b]

        X_w = rA[:, None] + RA @ quad_mesh.X_q.T
        X_w_local = RB.T @ (X_w - rB[:, None])

        narrow_mask = self._narrow_filter(X_w_local.T, rB, RB)
        if np.isscalar(narrow_mask) or (isinstance(narrow_mask, slice) and narrow_mask == slice(None)):
            active_indices = np.arange(quad_mesh.num_points)
        elif isinstance(narrow_mask, np.ndarray) and narrow_mask.dtype == bool:
            active_indices = np.where(narrow_mask)[0]
        else:
            active_indices = np.arange(quad_mesh.num_points)

        if len(active_indices) == 0:
            return np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

        k_n = cp.normal_stiffness
        mu = cp.friction_coefficient
        epsilon = cp.activation_distance
        regularizer = cp.quadrature_settings.get('regularizer', 1e-6)
        c_n = cp.quadrature_settings.get('damping', 0.0)

        Y_local = X_w_local[:, active_indices].T
        g, raw_grad, grad_norm, unit_normal, valid = sdf.query_batch(Y_local)

        valid_mask = valid & (k_n * np.maximum(epsilon - g, 0.0) > 1e-30)
        if not np.any(valid_mask):
            return np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

        idx_active = active_indices[valid_mask]
        gv = g[valid_mask]
        rg = raw_grad[valid_mask]
        gn = grad_norm[valid_mask]
        p = k_n * np.maximum(epsilon - gv, 0.0)

        xw_v = X_w[:, idx_active].T
        wv = quad_mesh.w_q[idx_active]
        off_A = xw_v - rA[None, :]

        # Normal penalty force
        fn = p[:, None] * rg
        p_hat = p * gn

        # Normal damping force: c_n * v_n * gradient
        if c_n > 0:
            xw_v_T = xw_v.T
            uA_damp = vA[:, None] + np.cross(wA, xw_v_T - rA[:, None], axis=0)
            uB_damp = vB[:, None] + np.cross(wB, xw_v_T - rB[:, None], axis=0)
            vn = np.sum((uA_damp - uB_damp).T * rg, axis=1)
            # Viscous damping: opposes motion in both directions (approach AND separation)
            # fd = c * (-vn) * grad, where vn = normal component of relative velocity
            # Positive vn = separating, negative vn = approaching
            # -vn is positive when approaching (compression) and negative when separating (tension)
            fd = c_n * (-vn)[:, None] * rg * (p > 0)[:, None]
            fn += fd

        uA = vA[None, :] + np.cross(wA[None, :], off_A)
        off_B = xw_v - rB[None, :]
        uB = vB[None, :] + np.cross(wB[None, :], off_B)
        u_rel = uA - uB

        un = rg / np.maximum(gn[:, None], 1e-30)
        udotn = np.sum(u_rel * un, axis=1, keepdims=True)
        u_t = u_rel - udotn * un
        u_t_norm = np.linalg.norm(u_t, axis=1)
        denom = np.sqrt(u_t_norm**2 + regularizer**2)
        ft = np.zeros_like(fn)
        mu_active = mu * p_hat > 0
        if np.any(mu_active):
            ft[mu_active] = -mu * p_hat[mu_active, None] * u_t[mu_active] / denom[mu_active, None]

        fq = fn + ft
        F_A = np.sum(wv[:, None] * fq, axis=0)
        tau_A = np.sum(wv[:, None] * np.cross(off_A, fq), axis=0)

        F_B = -F_A
        tau_B = -(tau_A + np.cross(rA - rB, F_A))
        return F_A, tau_A, F_B, tau_B

    def compute_Q_contact(self, model, state):
        nb_m = model.num_movable
        Q = np.zeros(6 * nb_m)
        for cp in model.contacts.values():
            F_A, tau_A, F_B, tau_B = self.evaluate(cp, state, model)
            self._add_to_Q(Q, model, cp.body_a_id, F_A, tau_A)
            self._add_to_Q(Q, model, cp.body_b_id, F_B, tau_B)
        return Q

    def _add_to_Q(self, Q, model, body_id, F, tau):
        for i, body in enumerate(model.movable_bodies):
            if body.id == body_id:
                Q[6*i:6*i+3] += F
                Q[6*i+3:6*i+6] += tau
                return
