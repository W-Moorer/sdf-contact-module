import numpy as np
from enum import Enum
from dataclasses import dataclass
from ..dynamics.assembler import DynamicsAssembler
from ..dynamics.constraints import (_all_residuals, num_joint_constraints,
                                    num_drive_constraints)
from ..dynamics.mass_matrix import build_mass_matrix, build_gyroscopic
from ..geometry.box_sdf import BoxSDF
from ..sdf_grid import AnalyticPlaneSDF


class ContactPhase(Enum):
    """Per-contact-region contact phase state machine."""
    INACTIVE = 0       # Far from contact, can use full nominal dt
    APPROACHING = 1    # A trial step predicts first contact; TOI localization needed
    IMPACTING = 2      # First-contact neighbourhood; very small dt
    ACTIVE = 3         # Contact force is active; use reduced dt
    STABLE = 4         # Contact is persistent; dt may grow gradually
    LEAVING = 5        # Contact is separating; grow dt / reset after debounce


@dataclass
class ContactRegionState:
    """State attached to one contact region.

    In the current IR the natural region key is ContactPair.id. If a single
    ContactPair later contains several disconnected patches, this dataclass can
    be keyed by a derived region id such as (cp_id, component_id).
    """
    cp_id: int
    phase: ContactPhase = ContactPhase.INACTIVE
    dt_current: float = 1e-3
    stable_count: int = 0
    release_count: int = 0
    min_dt_hits: int = 0
    last_gap: float = np.inf
    last_penetration: float = 0.0
    last_normal_velocity: float = 0.0
    last_toi: float = np.inf


@dataclass
class ContactEvent:
    cp_id: int
    toi: float
    penetration: float
    min_dt_hit: bool = False


def _pseudo_inverse_solve(A, b, tol=1e-9):
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    smax = S[0] if len(S) > 0 else 0.0
    keep = S > tol * max(smax, 1.0)
    return Vt.T[:, keep] @ ((U[:, keep].T @ b) / S[keep])


class BackwardEulerIntegrator:
    def __init__(self, model, contact_engine=None,
                 max_iter=20, tol=1e-8, line_search=True,
                 contact_penetration_tol=1e-6,
                 contact_dt_min=1e-8,
                 contact_dt_impact=1e-7,
                 contact_dt_active=1e-5,
                 contact_dt_max=1e-4,
                 contact_dt_normal=1e-3,
                 contact_dt_growth=1.5,
                 contact_stable_steps=5,
                 contact_release_steps=3):
        self.model = model
        self.contact_engine = contact_engine
        self.assembler = DynamicsAssembler(model, contact_engine=contact_engine)
        self.max_iter = max_iter
        self.tol = tol
        self.line_search = line_search
        self.history = []
        self._temp_state = None

        # Contact-time stepping parameters.
        self._dt_min = float(contact_dt_min)
        self._dt_impact = float(contact_dt_impact)
        self._dt_contact = float(contact_dt_active)
        self._dt_max_contact = float(contact_dt_max)
        self._dt_normal = float(contact_dt_normal)
        self._dt_growth = float(contact_dt_growth)
        self._pen_max = float(contact_penetration_tol)
        self._release_gap = max(5.0 * self._pen_max, 1e-9)
        self._stable_vn_tol = 1e-5
        self._stable_rel_pen_change = 5e-2
        self._stable_abs_pen_change = self._pen_max
        self._contact_stable_steps = int(contact_stable_steps)
        self._contact_release_steps = int(contact_release_steps)
        self._max_toi_bisect_iter = 60

        # Residual evaluation uses this temporary filter. None means all contact
        # pairs are allowed; an empty set means no contact forces are allowed.
        self._active_contact_ids_for_residual = None

        # Contact state machine. Kept in the old names too for compatibility
        # with any external debugging code that reads these attributes.
        self._contact_states = {}
        self._contact_phases = {}
        self._contact_dt = {}
        if contact_engine is not None:
            for cp in model.contacts.values():
                self._contact_states[cp.id] = ContactRegionState(cp_id=cp.id,
                                                                  dt_current=self._dt_normal)
                self._contact_phases[cp.id] = ContactPhase.INACTIVE
                self._contact_dt[cp.id] = self._dt_normal

    # ------------------------------------------------------------------
    # Public stepping API
    # ------------------------------------------------------------------
    def integrate(self, state, T, dt, callback=None):
        """Integrate until state.t + T using variable accepted substeps.

        The old implementation used int(T / dt) and could silently simulate less
        physical time when contact shortened a step. This loop uses actual time.
        """
        t_end = state.t + T
        i = 0
        n_nominal = max(1, int(np.ceil(T / dt)))
        while state.t < t_end - 1e-15:
            h_req = min(dt, t_end - state.t)
            h_done = self.step(state, h_req)
            if h_done <= 0:
                # Failsafe to avoid a zero-progress loop.
                h_done = min(self._dt_min, t_end - state.t)
                self._step_be(state, h_done, record=True, active_contact_ids=self._active_contact_ids())
            if callback:
                callback(state, i, n_nominal)
            i += 1
        return self.history

    def step(self, state, dt):
        """Advance by at most dt and return the accepted physical time step."""
        if self.model.num_movable == 0:
            state.t += dt
            self._record(state)
            return dt

        if self.contact_engine is None or len(self.model.contacts) == 0:
            self._step_be(state, dt, record=True, active_contact_ids=None)
            return dt

        # If the simulation starts with penetration, attach the state machine to
        # those contact regions immediately rather than waiting for TOI logic.
        self._activate_initially_penetrating_regions(state)

        # Existing active/stable contacts impose a local cap. New inactive
        # regions are pre-detected inside that cap; otherwise a new patch could
        # tunnel while another patch is using a reduced dt.
        h_cap = min(dt, self._active_contact_dt_cap())

        event = self._find_earliest_new_contact(state, h_cap)
        if event is not None:
            active_ids_before_event = self._active_contact_ids()
            # During the approach-to-TOI leg, the new pair is deliberately not in
            # the contact force set; otherwise the penalty force can hide the
            # first crossing in the trial solve.
            self._step_be(state, event.toi, record=False,
                          active_contact_ids=active_ids_before_event)
            st = self._contact_states[event.cp_id]
            st.phase = ContactPhase.IMPACTING
            st.dt_current = min(self._dt_impact, self._dt_contact)
            st.last_toi = event.toi
            if event.min_dt_hit:
                st.min_dt_hits += 1
            self._update_last_contact_measurements(state)
            self._sync_legacy_contact_maps()
            self._record(state)
            return event.toi

        # No new first contact inside h_cap. Advance with forces only from the
        # currently active contact regions. If all regions are inactive, this is
        # an ordinary BE step without contact forces.
        active_ids = self._active_contact_ids()
        self._step_be(state, h_cap, record=False, active_contact_ids=active_ids)
        self._refresh_contact_states(state)
        self._sync_legacy_contact_maps()
        self._record(state)
        return h_cap

    # ------------------------------------------------------------------
    # Contact state queries and transitions
    # ------------------------------------------------------------------
    def _contact_all_inactive(self):
        return all(st.phase == ContactPhase.INACTIVE
                   for st in self._contact_states.values())

    def _active_contact_ids(self):
        return {cid for cid, st in self._contact_states.items()
                if st.phase in (ContactPhase.IMPACTING,
                                ContactPhase.ACTIVE,
                                ContactPhase.STABLE,
                                ContactPhase.LEAVING)}

    def _active_contact_dt_cap(self):
        active = [st.dt_current for st in self._contact_states.values()
                  if st.phase in (ContactPhase.IMPACTING,
                                  ContactPhase.ACTIVE,
                                  ContactPhase.STABLE,
                                  ContactPhase.LEAVING)]
        return min(active) if active else self._dt_normal

    def _sync_legacy_contact_maps(self):
        for cid, st in self._contact_states.items():
            self._contact_phases[cid] = st.phase
            self._contact_dt[cid] = st.dt_current

    def _measure_contacts(self, state, cp_ids=None):
        if self.contact_engine is None:
            return {}
        if hasattr(self.contact_engine, 'measure_all'):
            return self.contact_engine.measure_all(self.model, state, cp_ids=cp_ids)
        # Backward-compatible fallback for older contact engines.
        out = {}
        wanted = None if cp_ids is None else set(cp_ids)
        for cp in self.model.contacts.values():
            if wanted is not None and cp.id not in wanted:
                continue
            g = self.contact_engine.min_gap(self.model, state)
            out[cp.id] = type('Measure', (), {
                'cp_id': cp.id, 'gap': float(g), 'penetration': max(0.0, -float(g)),
                'active_count': int(g <= cp.activation_distance),
                'valid_count': 1, 'normal_velocity': 0.0, 'min_point_index': -1
            })()
        return out

    def _activate_initially_penetrating_regions(self, state):
        measures = self._measure_contacts(state)
        for cid, m in measures.items():
            st = self._contact_states[cid]
            st.last_gap = m.gap
            st.last_penetration = m.penetration
            st.last_normal_velocity = m.normal_velocity
            if st.phase == ContactPhase.INACTIVE and m.gap <= 0.0:
                st.phase = ContactPhase.IMPACTING
                st.dt_current = min(self._dt_impact, self._dt_contact)

    def _update_last_contact_measurements(self, state):
        measures = self._measure_contacts(state)
        for cid, m in measures.items():
            st = self._contact_states[cid]
            st.last_gap = m.gap
            st.last_penetration = m.penetration
            st.last_normal_velocity = m.normal_velocity

    def _refresh_contact_states(self, state):
        measures = self._measure_contacts(state)
        for cid, m in measures.items():
            st = self._contact_states[cid]
            old_pen = st.last_penetration
            st.last_gap = m.gap
            st.last_penetration = m.penetration
            st.last_normal_velocity = m.normal_velocity

            if st.phase == ContactPhase.INACTIVE:
                continue

            # Contact has ended only after a small debounce window. This prevents
            # reset during tiny rebound/oscillation around g=0.
            if m.gap > self._release_gap and m.normal_velocity >= -self._stable_vn_tol:
                st.release_count += 1
                st.phase = ContactPhase.LEAVING
                st.dt_current = min(self._dt_normal,
                                    max(st.dt_current * self._dt_growth, self._dt_contact))
                if st.release_count >= self._contact_release_steps:
                    self._reset_contact_state(st)
                continue

            st.release_count = 0

            if st.phase == ContactPhase.IMPACTING:
                st.phase = ContactPhase.ACTIVE
                st.dt_current = max(self._dt_contact, self._dt_min)
                continue

            # Stable contact is persistent contact with small normal relative
            # motion and small penetration change. The penetration magnitude
            # itself is not forced below 1e-7 because penalty contact may settle
            # at mg/kA, which is usually larger than that.
            pen_change = abs(m.penetration - old_pen)
            rel_base = max(abs(old_pen), self._pen_max)
            stable_pen = (pen_change <= self._stable_abs_pen_change or
                          pen_change / rel_base <= self._stable_rel_pen_change)
            stable_vn = abs(m.normal_velocity) <= self._stable_vn_tol
            if m.active_count > 0 and stable_pen and stable_vn:
                st.stable_count += 1
            else:
                st.stable_count = 0

            if st.stable_count >= self._contact_stable_steps:
                st.phase = ContactPhase.STABLE
                st.dt_current = min(self._dt_max_contact,
                                    self._dt_normal,
                                    max(st.dt_current * self._dt_growth, self._dt_contact))
            else:
                st.phase = ContactPhase.ACTIVE
                st.dt_current = max(self._dt_contact, self._dt_min)

    def _reset_contact_state(self, st):
        st.phase = ContactPhase.INACTIVE
        st.dt_current = self._dt_normal
        st.stable_count = 0
        st.release_count = 0
        st.last_toi = np.inf
        st.min_dt_hits = 0

    # ------------------------------------------------------------------
    # Trial stepping and time-of-impact localization
    # ------------------------------------------------------------------
    def _trial_state(self, state, dt, active_contact_ids):
        trial = state.copy()
        self._step_be(trial, dt, record=False, active_contact_ids=active_contact_ids)
        return trial

    def _find_earliest_new_contact(self, state, dt):
        inactive_ids = [cid for cid, st in self._contact_states.items()
                        if st.phase == ContactPhase.INACTIVE]
        if not inactive_ids or dt <= 0:
            return None

        active_ids = self._active_contact_ids()
        start_measures = self._measure_contacts(state, cp_ids=inactive_ids)

        # Full-step pre-detection: run the candidate step, then inspect whether
        # any inactive contact region would cross into penetration.
        trial = self._trial_state(state, dt, active_ids)
        end_measures = self._measure_contacts(trial, cp_ids=inactive_ids)

        events = []
        for cid in inactive_ids:
            m0 = start_measures.get(cid)
            m1 = end_measures.get(cid)
            if m0 is None or m1 is None:
                continue
            # Already penetrating is handled before this function. Here we only
            # need first crossings inside the proposed step.
            crosses = (m0.gap > 0.0 and m1.gap <= 0.0)
            too_deep_after_trial = (m0.gap > 0.0 and m1.penetration > self._pen_max)
            if crosses or too_deep_after_trial:
                event = self._bisect_toi_pair(state, cid, dt, active_ids,
                                              high_measure=m1)
                if event is not None:
                    events.append(event)

        if not events:
            return None
        return min(events, key=lambda e: e.toi)

    def _bisect_toi_pair(self, state, cp_id, dt, active_contact_ids, high_measure=None):
        """Locate first contact for one pair using BE trial solves.

        We keep a bracket [lo, hi] where lo is non-contact and hi is contact or
        predicted excessive penetration. The accepted first-contact step is hi,
        unless the bracket reaches contact_dt_min, in which case the minimum-step
        rule intentionally accepts the tiny step and lets the next step continue.
        """
        lo = 0.0
        hi = dt
        hi_measure = high_measure
        min_dt_hit = False

        for _ in range(self._max_toi_bisect_iter):
            if hi - lo <= self._dt_min:
                min_dt_hit = True
                break
            mid = 0.5 * (lo + hi)
            trial = self._trial_state(state, mid, active_contact_ids)
            m = self._measure_contacts(trial, cp_ids=[cp_id])[cp_id]
            if m.gap <= 0.0 or m.penetration > self._pen_max:
                hi = mid
                hi_measure = m
            else:
                lo = mid

        # Guarantee forward progress under the user-specified minimum-step rule.
        if hi < self._dt_min and dt >= self._dt_min:
            hi = self._dt_min
            min_dt_hit = True
            trial = self._trial_state(state, hi, active_contact_ids)
            hi_measure = self._measure_contacts(trial, cp_ids=[cp_id])[cp_id]

        if hi <= 0.0:
            return None

        penetration = 0.0 if hi_measure is None else hi_measure.penetration
        # In the normal case, bisection should keep first-contact penetration <=
        # tolerance. If min_dt_hit is true, a larger value is explicitly allowed
        # by the minimum-step escape rule.
        if penetration > self._pen_max and not min_dt_hit:
            # One more shrink attempt; if still too deep, keep the last safe
            # positive substep and let the next outer step re-detect.
            hi = max(lo, min(hi, self._dt_min))
            min_dt_hit = True
        return ContactEvent(cp_id=cp_id, toi=min(hi, dt),
                            penetration=penetration, min_dt_hit=min_dt_hit)

    # Backward-compatible helper kept for older scripts. The new implementation
    # does not use this linearized vertical estimate for TOI decisions.
    def _bisect_toi(self, state, dt_guess, tol=1e-6):
        gap0 = self.contact_engine.min_gap(self.model, state)
        vz_min = min((state.v[state.body_idx(b.id), 2]
                      for b in self.model.movable_bodies), default=0.0)
        if vz_min >= 0:
            return 0.0

        t_low, t_high = 0.0, min(dt_guess, self._dt_normal)
        for _ in range(40):
            t_mid = (t_low + t_high) / 2.0
            if t_mid - t_low < self._dt_min:
                break
            gap_mid = gap0 + vz_min * t_mid
            if gap_mid < tol:
                t_high = t_mid
            else:
                t_low = t_mid
        return t_low

    def _get_contact_state(self, state):
        results = {}
        for cid, m in self._measure_contacts(state).items():
            results[cid] = {'gap': m.gap, 'toi_gap': m.gap, 'vz': m.normal_velocity,
                            'penetration': m.penetration,
                            'phase': self._contact_states[cid].phase.name}
        return results

    def _estimate_vertical_gap_single(self, state, cp):
        ba = self.model.bodies.get(cp.body_a_id)
        bb = self.model.bodies.get(cp.body_b_id)
        if ba is None or bb is None:
            return None
        try:
            idx_a = state.body_idx(ba.id)
            idx_b = state.body_idx(bb.id)
        except KeyError:
            return None

        z_a_low = state.r[idx_a, 2]
        z_b_high = state.r[idx_b, 2]

        pair = self.contact_engine._pairs.get(cp.id) if self.contact_engine else None
        if pair:
            qm, sdf = pair
            z_a_low = state.r[idx_a, 2] + np.min(qm.X_q[:, 2])
            if isinstance(sdf, AnalyticPlaneSDF):
                z_b_high = state.r[idx_b, 2] + sdf.z_top
            elif isinstance(sdf, BoxSDF):
                z_b_high = state.r[idx_b, 2] + sdf.center[2] + sdf.half[2]
        return z_a_low - z_b_high

    def _estimate_vertical_gap(self, state):
        if not self.contact_engine:
            return None
        g_min = np.inf
        for cp in self.model.contacts.values():
            g = self._estimate_vertical_gap_single(state, cp)
            if g is not None:
                g_min = min(g_min, g)
        return g_min if np.isfinite(g_min) else None

    # ------------------------------------------------------------------
    # Backward Euler nonlinear solve
    # ------------------------------------------------------------------
    def _step_be(self, state, dt, record=True, active_contact_ids=None):
        q0 = state.copy()
        V0 = q0.pack_V()
        nc = num_joint_constraints(self.model) + num_drive_constraints(self.model)
        nb_m = self.model.num_movable
        n_v = 6 * nb_m

        old_filter = self._active_contact_ids_for_residual
        self._active_contact_ids_for_residual = active_contact_ids
        try:
            z = np.concatenate([V0, np.zeros(nc)])
            R = self._residual(z, q0, dt)
            for iteration in range(self.max_iter):
                err = np.max(np.abs(R)) if len(R) else 0.0
                if err < self.tol:
                    break
                K = self._tangent(z, q0, dt, R_cache=R, use_fd=False)
                dz_raw = self._solve_kkt(K, -R)
                alpha = self._line_search(z, dz_raw, q0, dt, R0=R) if self.line_search else 1.0
                z = z + alpha * dz_raw
                R = self._residual(z, q0, dt)

            V = z[:n_v]
            self._kinematic_update(q0, V, dt, state)
            state.t += dt
            if record:
                self._record(state)
        finally:
            self._active_contact_ids_for_residual = old_filter

    def _kinematic_update(self, q0, V, dt, state):
        state.r[:] = q0.r.copy()
        state.R[:] = q0.R.copy()
        state.v[:] = q0.v.copy()
        state.omega[:] = q0.omega.copy()
        state.t = q0.t

        tmp = 0
        for body in self.model.movable_bodies:
            idx = state.body_idx(body.id)
            state.v[idx] = V[6*tmp:6*tmp+3]
            state.omega[idx] = V[6*tmp+3:6*tmp+6]
            state.r[idx] = q0.r[idx] + dt * state.v[idx]

            dtheta = dt * state.omega[idx]
            dtheta_norm = np.linalg.norm(dtheta)
            if dtheta_norm > 1e-30:
                axis = dtheta / dtheta_norm
                c = np.cos(dtheta_norm)
                s = np.sin(dtheta_norm)
                dR = (c * np.eye(3)
                      + s * np.array([[0, -axis[2], axis[1]],
                                      [axis[2], 0, -axis[0]],
                                      [-axis[1], axis[0], 0]])
                      + (1 - c) * np.outer(axis, axis))
                state.R[idx] = dR @ q0.R[idx]
            else:
                state.R[idx] = q0.R[idx].copy()
            tmp += 1

    def _residual(self, z, q0, dt):
        nb_m = self.model.num_movable
        V = z[:6*nb_m]
        lam = z[6*nb_m:]

        if self._temp_state is None:
            self._temp_state = q0.copy()
        self._kinematic_update(q0, V, dt, self._temp_state)
        state = self._temp_state

        M = build_mass_matrix(self.model, state).toarray()
        C = build_gyroscopic(self.model, state)

        Q = np.zeros(6 * nb_m)
        g = np.asarray(self.model.gravity, dtype=float)
        tmp = 0
        for body in self.model.movable_bodies:
            Q[6*tmp:6*tmp+3] = body.mass * g
            tmp += 1

        if self.contact_engine is not None:
            active_ids = self._active_contact_ids_for_residual
            if hasattr(self.contact_engine, 'compute_Q_contact'):
                try:
                    Q += self.contact_engine.compute_Q_contact(self.model, state,
                                                               active_ids=active_ids)
                except TypeError:
                    # Older engine without active-id filter.
                    if active_ids is None:
                        Q += self.contact_engine.compute_Q_contact(self.model, state)
                    elif len(active_ids) > 0:
                        Q += self.contact_engine.compute_Q_contact(self.model, state)

        V0 = q0.pack_V()
        R_dyn = M @ (V - V0) / dt + C - Q

        nc = num_joint_constraints(self.model) + num_drive_constraints(self.model)
        J = None
        if nc > 0:
            from ..dynamics.constraints import eval_constraint_jacobian
            J_joint = eval_constraint_jacobian(self.model, state)
            nd = num_drive_constraints(self.model)
            J_drive = np.zeros((nd, 6 * nb_m))
            if nd > 0:
                from ..dynamics.drives import motion_drive_jacobian
                row = 0
                for drive in self.model.drives.values():
                    Jd = motion_drive_jacobian(drive, self.model, state)
                    J_drive[row:row+1] = Jd
                    row += 1
            J = np.vstack([J_joint, J_drive]) if len(J_joint) > 0 and len(J_drive) > 0 else (J_joint if len(J_joint) > 0 else J_drive)
            R_dyn -= J.T @ lam

        Phi = _all_residuals(self.model, state)
        Phi_scaled = Phi / dt

        self._cache_M = M
        self._cache_J = J
        self._cache_nb_m = nb_m
        self._cache_nc = nc

        any_contact = False
        if self.contact_engine is not None:
            ids = self._active_contact_ids_for_residual
            if ids is None:
                # All pairs allowed.
                any_contact = self.contact_engine.min_gap(self.model, state) < 0.001
            elif len(ids) > 0:
                measures = self._measure_contacts(state, cp_ids=ids)
                any_contact = any(m.gap < 0.001 for m in measures.values())
        self._cache_any_contact = any_contact

        return np.concatenate([R_dyn, Phi_scaled])

    def _contact_tangent(self, state, dt):
        """Estimate diagonal contact stiffness contribution to KKT tangent.
        
        For KORDER=n: ∂F/∂V = k_n * n * δ^(n-1) * area * dt  (where δ = ε - g)
        For n=1: ∂F/∂V = k_n * area * dt (δ-independent)
        For n=2: ∂F/∂V = 2*k_n*δ*area*dt (δ-dependent, zero when δ=0)
        """
        if self.contact_engine is None:
            return None
        nb_m = self.model.num_movable
        Kc = np.zeros(6 * nb_m)
        active_ids = self._active_contact_ids_for_residual
        for cp in self.model.contacts.values():
            if active_ids is not None and cp.id not in active_ids:
                continue
            k_n = cp.normal_stiffness
            k_order = getattr(cp, 'k_order', 1)
            c_n = cp.quadrature_settings.get('damping', 0.0)
            pair = self.contact_engine._pairs.get(cp.id)
            if pair is None:
                continue
            qm, _ = pair
            area = np.sum(qm.w_q)
            # Estimate penetration δ from mesh velocity
            ba = self.model.bodies.get(cp.body_a_id)
            if ba is None:
                continue
            try:
                idx_all = state.body_idx(ba.id)
                vel = np.linalg.norm(state.v[idx_all])
            except (KeyError, IndexError):
                vel = 1.0
            delta_est = max(vel * dt, 1e-10)
            # With F = (k_n/A)*δ^k_order*A = k_n*δ^k_order:
            # For any k_order: ∂F/∂V = k_n * k_order * δ^(k_order-1) * dt
            if k_order == 1:
                stiff = k_n * dt
            else:
                stiff = k_n * k_order * (delta_est ** (k_order - 1)) * dt
                stiff += 5e2  # regularization for δ→0 where ∂F/∂V→0
            damp = c_n * area
            total = stiff + damp
            if ba is None:
                continue
            try:
                idx_all = state.body_idx(ba.id)
            except KeyError:
                continue
            for movable_i, body in enumerate(self.model.movable_bodies):
                if body.id == ba.id:
                    for d in range(3):
                        Kc[6*movable_i + d] += total
                    break
        return Kc

    def _tangent(self, z, q0, dt, R_cache=None, use_fd=False):
        nb_m = self._cache_nb_m
        nc = self._cache_nc
        nv = 6 * nb_m
        n = nv + nc

        if use_fd:
            eps = 1e-7
            R0 = self._residual(z, q0, dt) if R_cache is None else R_cache
            K = np.zeros((n, n))
            for i in range(n):
                zp = z.copy()
                zp[i] += eps
                Rp = self._residual(zp, q0, dt)
                K[:, i] = (Rp - R0) / eps
            return K

        K = np.zeros((n, n))
        K[:nv, :nv] = self._cache_M / dt
        Kc = self._contact_tangent(self._temp_state, dt)
        if Kc is not None:
            for i in range(nv):
                K[i, i] += Kc[i]
        if nc > 0 and self._cache_J is not None:
            K[:nv, nv:] = -self._cache_J.T
            K[nv:, :nv] = self._cache_J
        return K

    def _line_search(self, z, dz, q0, dt, R0=None):
        R0_norm = np.linalg.norm(R0) if R0 is not None else np.linalg.norm(self._residual(z, q0, dt))
        for alpha in [1.0, 0.5, 0.25, 0.125, 0.0625]:
            z_new = z + alpha * dz
            R_norm = np.linalg.norm(self._residual(z_new, q0, dt))
            if R_norm < R0_norm:
                return alpha
        return 1.0

    @staticmethod
    def _solve_kkt(K, rhs):
        return np.linalg.lstsq(K, rhs, rcond=1e-12)[0]

    def _record(self, state):
        from ..dynamics.joints_revolute import revolute_joint_coordinate
        record = {'t': state.t}
        for joint in self.model.joints.values():
            if joint.type.name == 'REVOLUTE':
                record[f'joint_{joint.name}_theta'] = revolute_joint_coordinate(
                    joint, self.model, state)
        for body in self.model.bodies.values():
            idx = state.body_idx(body.id)
            record[f'body_{body.name}_r'] = state.r[idx].copy()
            record[f'body_{body.name}_v'] = state.v[idx].copy()
        record['contact_phases'] = {cid: st.phase.name for cid, st in self._contact_states.items()}
        record['contact_dt'] = {cid: st.dt_current for cid, st in self._contact_states.items()}
        record['contact_gap'] = {cid: st.last_gap for cid, st in self._contact_states.items()}
        record['contact_penetration'] = {cid: st.last_penetration for cid, st in self._contact_states.items()}
        self.history.append(record)
