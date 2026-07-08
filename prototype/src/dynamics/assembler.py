import numpy as np
from scipy.sparse import csr_matrix
from .mass_matrix import build_mass_matrix, build_gyroscopic
from .constraints import (eval_constraint_jacobian, compute_bc,
                          compute_bc_from_jacobian,
                          num_joint_constraints, num_drive_constraints)


class DynamicsAssembler:
    def __init__(self, model, contact_engine=None):
        self.model = model
        self.contact_engine = contact_engine

    def assemble(self, state, dt, J_precomputed=None, bc_precomputed=None):
        M = build_mass_matrix(self.model, state)
        C = build_gyroscopic(self.model, state)
        Q_ext = self._external_forces(state)

        if self.contact_engine is not None:
            Q_ext += self.contact_engine.compute_Q_contact(self.model, state)
        Q_ext += self._force_elements(state)

        nc = num_joint_constraints(self.model) + num_drive_constraints(self.model)
        nv = 6 * self.model.num_movable

        if nc > 0 and nv > 0:
            if J_precomputed is not None:
                J = J_precomputed
            else:
                from .drives import motion_drive_jacobian
                J_joint = eval_constraint_jacobian(self.model, state)
                nd = num_drive_constraints(self.model)
                J_drive = np.zeros((nd, nv))
                row = 0
                for drive in self.model.drives.values():
                    Jd = motion_drive_jacobian(drive, self.model, state)
                    J_drive[row:row+1] = Jd
                    row += 1
                J = np.vstack([J_joint, J_drive]) if len(J_joint) > 0 and len(J_drive) > 0 else (J_joint if len(J_joint) > 0 else J_drive)

            J_sparse = csr_matrix(J)

            if bc_precomputed is not None:
                b_c = bc_precomputed
            else:
                b_c = compute_bc_from_jacobian(J, self.model, state, dt)
        else:
            J_sparse = csr_matrix((0, nv))
            b_c = np.zeros(0)

        return csr_matrix(M), C, Q_ext, J_sparse, b_c

    def _external_forces(self, state):
        nb_m = self.model.num_movable
        F = np.zeros(6 * nb_m)
        g = np.asarray(self.model.gravity, dtype=float)
        tmp_idx = 0
        for body in self.model.movable_bodies:
            F[6*tmp_idx:6*tmp_idx+3] = body.mass * g
            tmp_idx += 1
        return F

    def _force_elements(self, state):
        nb_m = self.model.num_movable
        if nb_m == 0:
            return np.zeros(0)
        Q = np.zeros(6 * nb_m)
        for fe in self.model.forces.values():
            F_i, tau_i, F_j, tau_j = fe.evaluate(state, self.model)
            self._add_to_Q(Q, fe.frame_i_id, F_i, tau_i, state)
            self._add_to_Q(Q, fe.frame_j_id, F_j, tau_j, state)
        for te in self.model.torques.values():
            tau = te.evaluate(state, self.model)
            body = self.model.bodies[te.body_id]
            for i, mb in enumerate(self.model.movable_bodies):
                if mb.id == body.id:
                    Q[6*i+3:6*i+6] += tau
                    break
        return Q

    def _add_to_Q(self, Q, frame_id, F, tau, state):
        frame = self.model.frames.get(frame_id)
        if frame is None:
            return
        body = self.model.bodies.get(frame.body_id)
        if body is None:
            return
        for i, mb in enumerate(self.model.movable_bodies):
            if mb.id == body.id:
                Q[6*i:6*i+3] += F
                Q[6*i+3:6*i+6] += tau
                return
