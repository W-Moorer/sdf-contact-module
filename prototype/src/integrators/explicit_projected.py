import numpy as np
from ..dynamics.assembler import DynamicsAssembler
from ..solvers.kkt import solve_kkt
from .projection import position_projection, velocity_projection


class ExplicitProjectedIntegrator:
    def __init__(self, model, alpha=0.1, beta=0.1, contact_engine=None):
        self.model = model
        self.contact_engine = contact_engine
        self.assembler = DynamicsAssembler(model, contact_engine=contact_engine)
        self.alpha = alpha
        self.beta = beta
        self.history = []

    def step(self, state, dt):
        if self.model.num_movable == 0:
            state.t += dt
            return

        M_sp, C, Q, J, b_c = self.assembler.assemble(state, dt)

        # Cache dense M, J, M_inv for reuse across kkt + projections
        M_dense = M_sp.toarray()
        J_dense = J.toarray() if hasattr(J, 'toarray') else np.asarray(J)
        M_inv = np.linalg.inv(M_dense)

        acc, lam = solve_kkt(M_dense, C, Q, J_dense, b_c, M_inv=M_inv)

        tmp_idx = 0
        for body in self.model.movable_bodies:
            idx = state.body_idx(body.id)
            state.v[idx] += dt * acc[6*tmp_idx:6*tmp_idx+3]
            state.omega[idx] += dt * acc[6*tmp_idx+3:6*tmp_idx+6]
            tmp_idx += 1

        tmp_idx = 0
        for body in self.model.movable_bodies:
            idx = state.body_idx(body.id)
            state.r[idx] += dt * state.v[idx]
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
                state.R[idx] = dR @ state.R[idx]
            tmp_idx += 1

        position_projection(self.model, state,
                            J_precomputed=J_dense, M_precomputed=M_dense, M_inv_precomputed=M_inv)
        velocity_projection(self.model, state,
                            J_precomputed=J_dense, M_precomputed=M_dense, M_inv_precomputed=M_inv)

        state.t += dt

        self._record(state)

    def _record(self, state):
        from ..dynamics.joints_revolute import revolute_joint_coordinate
        record = {'t': state.t}
        for joint in self.model.joints.values():
            if joint.type.name == 'REVOLUTE':
                theta = revolute_joint_coordinate(joint, self.model, state)
                record[f'joint_{joint.name}_theta'] = theta
        for body in self.model.bodies.values():
            idx = state.body_idx(body.id)
            record[f'body_{body.name}_r'] = state.r[idx].copy()
            record[f'body_{body.name}_v'] = state.v[idx].copy()
        self.history.append(record)

    def integrate(self, state, T, dt, callback=None):
        n_steps = int(T / dt)
        for i in range(n_steps):
            self.step(state, dt)
            if callback:
                callback(state, i, n_steps)
        return self.history
