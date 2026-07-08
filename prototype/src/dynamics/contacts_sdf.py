import numpy as np
from dataclasses import dataclass
from ..sdf_grid import AnalyticPlaneSDF
from ..geometry.bvh import AABB_BVH


@dataclass
class ContactMeasure:
    """Geometric/contact-kinematic measure for one registered contact pair."""
    cp_id: int
    gap: float
    penetration: float
    active_count: int
    valid_count: int
    normal_velocity: float = 0.0
    min_point_index: int = -1


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
        self._narrow_margin = 0.01

    def register_pair(self, cp_id, quadrature_mesh, sdf_grid,
                      aabb_b_half=None, narrow_margin=0.01):
        self._pairs[cp_id] = (quadrature_mesh, sdf_grid)
        self._aabb_a[cp_id] = AABB.from_points(quadrature_mesh.X_q)
        self._aabb_b[cp_id] = aabb_b_half
        self._narrow_margin = narrow_margin
        # Build BVH over quadrature points in body-A local frame
        self._bvh[cp_id] = AABB_BVH(quadrature_mesh.X_q, leaf_size=16)

    @property
    def _bvh(self):
        if not hasattr(self, '__bvh'):
            self.__bvh = {}
        return self.__bvh

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

    def _narrow_filter(self, X_w_local, rB=None, RB=None):
        # X_w_local is already expressed in body-B local coordinates; do not subtract rB again.
        margin = self._narrow_margin
        if RB is None:
            return slice(None)
        dz = X_w_local[:, 2]
        mask = np.abs(dz) < margin
        return mask

    def _pair_query(self, cp, state, model):
        """Return quadrature world points and SDF query arrays for exactly one pair."""
        pair_key = cp.id
        if pair_key not in self._pairs:
            return None

        quad_mesh, sdf = self._pairs[pair_key]
        body_a = model.bodies[cp.body_a_id]
        body_b = model.bodies[cp.body_b_id]
        idx_a = state.body_idx(body_a.id)
        idx_b = state.body_idx(body_b.id)

        rA, RA = state.r[idx_a], state.R[idx_a]
        rB, RB = state.r[idx_b], state.R[idx_b]

        # BVH culling
        bvh = self._bvh.get(pair_key)
        if bvh is not None and hasattr(sdf, 'bmin') and hasattr(sdf, 'bmax'):
            margin = max(self._narrow_margin, 0.05)
            bmin_local = RA.T @ (sdf.bmin - rA) - margin
            bmax_local = RA.T @ (sdf.bmax - rA) + margin
            cull_idx = bvh.query(bmin_local, bmax_local)
        else:
            cull_idx = np.arange(quad_mesh.num_points)

        if len(cull_idx) == 0:
            return None

        X_w = rA[:, None] + RA @ quad_mesh.X_q[cull_idx].T
        Y_local = RB.T @ (X_w - rB[:, None])
        g, raw_grad, grad_norm, unit_normal, valid = sdf.query_batch(Y_local.T)
        return quad_mesh, sdf, X_w, Y_local, g, raw_grad, grad_norm, unit_normal, valid, cull_idx

    def measure_pair(self, cp, state, model, activation_distance=None):
        """Measure signed gap and normal relative velocity for one contact pair.

        Negative gap means penetration. The result is per ContactPair, which is the
        current code's natural representation of a contact region.
        """
        q = self._pair_query(cp, state, model)
        if q is None:
            return ContactMeasure(cp.id, np.inf, 0.0, 0, 0, 0.0, -1)

        quad_mesh, sdf, X_w, Y_local, g, raw_grad, grad_norm, unit_normal, valid = q[:9]
        if not np.any(valid):
            return ContactMeasure(cp.id, np.inf, 0.0, 0, 0, 0.0, -1)

        valid_idx = np.where(valid)[0]
        g_valid = g[valid]
        local_argmin = int(np.argmin(g_valid))
        min_idx = int(valid_idx[local_argmin])
        g_min = float(g[min_idx])
        penetration = max(0.0, -g_min)

        eps = cp.activation_distance if activation_distance is None else activation_distance
        active_count = int(np.sum(valid & (g <= eps)))

        # normal relative velocity at the minimum-gap quadrature point
        body_a = model.bodies[cp.body_a_id]
        body_b = model.bodies[cp.body_b_id]
        idx_a = state.body_idx(body_a.id)
        idx_b = state.body_idx(body_b.id)
        rA, rB = state.r[idx_a], state.r[idx_b]
        RB = state.R[idx_b]
        vA, wA = state.v[idx_a], state.omega[idx_a]
        vB, wB = state.v[idx_b], state.omega[idx_b]
        xw = X_w[:, min_idx]
        n_world = RB @ unit_normal[min_idx]
        n_norm = np.linalg.norm(n_world)
        if n_norm > 1e-30:
            n_world = n_world / n_norm
        uA = vA + np.cross(wA, xw - rA)
        uB = vB + np.cross(wB, xw - rB)
        normal_velocity = float(np.dot(uA - uB, n_world))

        return ContactMeasure(cp.id, g_min, penetration, active_count,
                              int(np.sum(valid)), normal_velocity, min_idx)

    def measure_all(self, model, state, cp_ids=None):
        wanted = None if cp_ids is None else set(cp_ids)
        out = {}
        for cp in model.contacts.values():
            if wanted is not None and cp.id not in wanted:
                continue
            out[cp.id] = self.measure_pair(cp, state, model)
        return out

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

        # --- BVH broad-phase culling ---
        # Get cube's SDF bounding box in world (cube body frame is world)
        if hasattr(sdf, 'bmin') and hasattr(sdf, 'bmax'):
            bmin_w, bmax_w = sdf.bmin, sdf.bmax
        else:
            bmin_w, bmax_w = np.full(3, -1.0), np.full(3, 1.0)

        margin = max(self._narrow_margin, 0.05)
        # Transform cube bbox to body-A local frame
        bmin_local = RA.T @ (bmin_w - rA) - margin
        bmax_local = RA.T @ (bmax_w - rA) + margin

        bvh = self._bvh.get(pair_key)
        if bvh is not None:
            bvh_indices = bvh.query(bmin_local, bmax_local)
        else:
            bvh_indices = np.arange(quad_mesh.num_points)

        if len(bvh_indices) == 0:
            return np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

        # --- Transform only BVH-culled points ---
        X_w = rA[:, None] + RA @ quad_mesh.X_q[bvh_indices].T
        X_w_local = RB.T @ (X_w - rB[:, None])

        narrow_mask = self._narrow_filter(X_w_local.T, rB, RB)
        if np.isscalar(narrow_mask) or (isinstance(narrow_mask, slice) and narrow_mask == slice(None)):
            local_active = np.arange(len(bvh_indices))
        elif isinstance(narrow_mask, np.ndarray) and narrow_mask.dtype == bool:
            local_active = np.where(narrow_mask)[0]
        else:
            local_active = np.arange(len(bvh_indices))

        if len(local_active) == 0:
            return np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

        # Map back to global quad mesh indices
        active_indices = bvh_indices[local_active]
        k_n = cp.normal_stiffness
        mu = cp.friction_coefficient
        epsilon = cp.activation_distance
        regularizer = cp.quadrature_settings.get('regularizer', 1e-6)
        c_n = cp.quadrature_settings.get('damping', 0.0)

        Y_local = X_w_local[:, local_active].T
        g, raw_grad, grad_norm, unit_normal, valid = sdf.query_batch(Y_local)

        valid_mask = valid & (k_n * np.maximum(epsilon - g, 0.0) > 1e-30)
        if not np.any(valid_mask):
            return np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)

        # idx_active_global: global quad mesh indices for active points
        idx_active_local = np.where(valid_mask)[0]
        idx_active_global = active_indices[valid_mask]
        gv = g[valid_mask]
        rg_local = raw_grad[valid_mask]
        gn = grad_norm[valid_mask]
        p = k_n * np.maximum(epsilon - gv, 0.0)

        # Forces are accumulated in world coordinates. SDF gradients are local to body B.
        rg = (RB @ rg_local.T).T

        # Compute world positions for active points only
        xw_active = rA[:, None] + RA @ quad_mesh.X_q[idx_active_global].T
        xw_v = xw_active.T
        wv = quad_mesh.w_q[idx_active_global]
        off_A = xw_v - rA[None, :]

        # Normal penalty force
        fn = p[:, None] * rg
        p_hat = p * gn

        # Normal damping force: active as long as gap < epsilon (even after penalty force ends)
        if c_n > 0:
            uA_damp = vA[None, :] + np.cross(wA[None, :], off_A)
            off_B_damp = xw_v - rB[None, :]
            uB_damp = vB[None, :] + np.cross(wB[None, :], off_B_damp)
            vn = np.sum((uA_damp - uB_damp) * rg, axis=1)
            damp_active = gv < epsilon
            fd = c_n * (-vn)[:, None] * rg * damp_active[:, None]
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

    def min_gap(self, model, state, cp_id=None):
        """Minimum signed distance. Negative = penetration."""
        if cp_id is not None:
            cp = model.contacts.get(cp_id)
            if cp is None:
                return np.inf
            return self.measure_pair(cp, state, model).gap

        g_min = np.inf
        for cp in model.contacts.values():
            g_min = min(g_min, self.measure_pair(cp, state, model).gap)
        return g_min if np.isfinite(g_min) else np.inf

    def compute_Q_contact(self, model, state, active_ids=None):
        nb_m = model.num_movable
        Q = np.zeros(6 * nb_m)
        wanted = None if active_ids is None else set(active_ids)
        for cp in model.contacts.values():
            if wanted is not None and cp.id not in wanted:
                continue
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
