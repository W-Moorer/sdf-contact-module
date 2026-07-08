import numpy as np
from ..dynamics.assembler import DynamicsAssembler
from ..dynamics.constraints import (_all_residuals, num_joint_constraints,
                                    num_drive_constraints)
from ..dynamics.mass_matrix import build_mass_matrix, build_gyroscopic


class BackwardEulerIntegrator:
    def __init__(self, model, contact_engine=None,
                 max_iter=20, tol=1e-8, line_search=True):
        self.model = model
        self.contact_engine = contact_engine
        self.assembler = DynamicsAssembler(model, contact_engine=contact_engine)
        self.max_iter = max_iter
        self.tol = tol
        self.line_search = line_search
        self.history = []

    def integrate(self, state, T, dt, callback=None):
        n_steps = int(T / dt)
        for i in range(n_steps):
            self.step(state, dt)
            if callback:
                callback(state, i, n_steps)
        return self.history

    def step(self, state, dt):
        if self.model.num_movable == 0:
            state.t += dt
            self._record(state)
            return

        q0 = state.copy()
        V0 = q0.pack_V()
        nc = num_joint_constraints(self.model) + num_drive_constraints(self.model)
        nb_m = self.model.num_movable

        z = np.concatenate([q0.pack_V(), np.zeros(nc)])

        for iteration in range(self.max_iter):
            R = self._residual(z, q0, dt)
            err = np.max(np.abs(R))
            if err < self.tol:
                break

            K = self._tangent(z, q0, dt)
            dz = np.linalg.solve(K, -R)

            if self.line_search:
                alpha = self._line_search(z, dz, q0, dt)
            else:
                alpha = 1.0

            z = z + alpha * dz

        V = z[:6*nb_m]
        self._kinematic_update(q0, V, dt, state)
        state.t += dt
        self._record(state)

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

        state = q0.copy()
        self._kinematic_update(q0, V, dt, state)

        M = build_mass_matrix(self.model, state).toarray()
        C = build_gyroscopic(self.model, state)

        Q = np.zeros(6 * nb_m)
        g = np.asarray(self.model.gravity, dtype=float)
        tmp = 0
        for body in self.model.movable_bodies:
            Q[6*tmp:6*tmp+3] = body.mass * g
            tmp += 1

        if self.contact_engine is not None:
            Q += self.contact_engine.compute_Q_contact(self.model, state)

        V0 = q0.pack_V()
        R_dyn = M @ (V - V0) / dt + C - Q

        nc = num_joint_constraints(self.model) + num_drive_constraints(self.model)
        if nc > 0:
            from ..dynamics.constraints import eval_constraint_jacobian, eval_drive_constraints
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
        return np.concatenate([R_dyn, Phi])

    def _tangent(self, z, q0, dt):
        eps = 1e-7
        R0 = self._residual(z, q0, dt)
        n = len(z)
        K = np.zeros((n, n))

        for i in range(n):
            zp = z.copy()
            zp[i] += eps
            Rp = self._residual(zp, q0, dt)
            K[:, i] = (Rp - R0) / eps

        return K

    def _line_search(self, z, dz, q0, dt):
        R0_norm = np.linalg.norm(self._residual(z, q0, dt))
        for alpha in [1.0, 0.5, 0.25, 0.125, 0.0625]:
            z_new = z + alpha * dz
            R_norm = np.linalg.norm(self._residual(z_new, q0, dt))
            if R_norm < R0_norm:
                return alpha
        return 1.0

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
        self.history.append(record)
