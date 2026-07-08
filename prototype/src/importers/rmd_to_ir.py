import numpy as np
from ..ir import (SystemModel, RigidBody, Frame, FixedJoint, RevoluteJoint,
                  Drive, DriveMode, ConstantFunction, LinearFunction, SinFunction,
                  ContactPair, ForceElement, TorqueElement)
from .rmd_lexer import RMDLexer
from .rmd_parser import RMDParser


JOINT_TYPE_MAP = {
    'fixed': 'FIXED',
    'revolute': 'REVOLUTE',
    'cylindrical': 'REVOLUTE',
    'translational': 'TRANSLATIONAL',
}


def euler_to_rotation(psi, theta, phi):
    """Body-fixed 3-1-3 Euler angles → rotation matrix."""
    c1, s1 = np.cos(psi), np.sin(psi)
    c2, s2 = np.cos(theta), np.sin(theta)
    c3, s3 = np.cos(phi), np.sin(phi)
    R = np.array([
        [c1*c3 - s1*c2*s3, -c1*s3 - s1*c2*c3, s1*s2],
        [s1*c3 + c1*c2*s3, -s1*s3 + c1*c2*c3, -c1*s2],
        [s2*s3, s2*c3, c2],
    ])
    return R


class RMDToIR:
    def __init__(self, gravity=None, length_scale=1.0):
        self.gravity = np.array([0, 0, -9.81]) if gravity is None else np.array(gravity, dtype=float)
        self.length_scale = length_scale
        self.warnings = []
        self._marker_map = {}
        self._body_map = {}
        self._joint_map = {}
        self._surface_map = {}
        self._contact_map = {}
        self._extracted_objs = []

    def import_rmd(self, text):
        lexer = RMDLexer()
        entities = lexer.parse(text)
        parser = RMDParser()
        raw = parser.parse(entities)
        self.warnings = lexer.warnings + parser.warnings
        if raw.units.get('LENGTH', '').upper() in ('MILLIMETER', 'MM'):
            self.length_scale = 0.001
        if raw.gravity is not None:
            self.gravity = np.array(raw.gravity, dtype=float) * self.length_scale
        return self._normalize(raw)

    def load_file(self, path):
        lexer = RMDLexer()
        entities = lexer.load(path)
        parser = RMDParser()
        raw = parser.parse(entities)
        self.warnings = lexer.warnings + parser.warnings
        if raw.units.get('LENGTH', '').upper() in ('MILLIMETER', 'MM'):
            self.length_scale = 0.001
        if raw.gravity is not None:
            self.gravity = np.array(raw.gravity, dtype=float) * self.length_scale
        return self._normalize(raw)

    def _normalize(self, raw):
        model = SystemModel()
        model.gravity = self.gravity.copy()

        self._body_map = {}
        self._marker_map = {}
        self._body_cm_pos = {}
        self._surface_map = {}
        self._contact_map = {}
        self._joint_map = {}

        self._add_bodies(raw, model)
        self._add_markers(raw, model)
        self._add_surfaces(raw, model)
        self._add_contacts(raw, model)
        self._add_joints(raw, model)
        self._add_motions(raw, model)
        self._add_axial_forces(raw, model)

        return model

    def _add_bodies(self, raw, model):
        for bid, rb in raw.bodies.items():
            pos = rb.initial_position
            rot = rb.initial_rotation

            cm_marker = raw.markers.get(rb.cm_marker_id) if rb.cm_marker_id else None
            if cm_marker and cm_marker.position:
                pos = cm_marker.position

            if pos:
                pos = [v * self.length_scale for v in pos]

            initial_R = np.eye(3)
            if rot and len(rot) >= 3:
                initial_R = euler_to_rotation(rot[0], rot[1], rot[2])

            ibody = RigidBody(
                name=rb.name,
                body_id=bid,
                mass=rb.mass,
                inertia_body=np.eye(3) if rb.inertia is None else np.array(rb.inertia, dtype=float) * (self.length_scale**2),
                initial_r=np.zeros(3) if pos is None else np.array(pos, dtype=float),
                initial_R=initial_R,
                fixed=rb.ground,
            )
            model.add_body(ibody)
            self._body_map[rb.name] = bid
            self._body_map[bid] = rb.name

    def _add_markers(self, raw, model):
        self._body_cm_pos = {}
        for bid, rb in raw.bodies.items():
            cm_marker = raw.markers.get(rb.cm_marker_id) if rb.cm_marker_id else None
            if cm_marker and cm_marker.position:
                self._body_cm_pos[bid] = np.array(cm_marker.position, dtype=float)
            else:
                self._body_cm_pos[bid] = np.zeros(3)

        next_fid = 0
        for mid, rm in raw.markers.items():
            body_id = rm.part_id
            if body_id not in model.bodies:
                continue

            pos = rm.position if rm.position else np.zeros(3)
            cm_offset = self._body_cm_pos.get(body_id, np.zeros(3))
            local_pos = (np.array(pos, dtype=float) - cm_offset) * self.length_scale

            if np.linalg.norm(local_pos) < 1e-30 and 'CM' in (rm.name or '').split('.')[-1]:
                continue

            ori = np.eye(3)
            if rm.orientation and len(rm.orientation) >= 3:
                ori = euler_to_rotation(rm.orientation[0], rm.orientation[1], rm.orientation[2])

            iframe = Frame(
                name=rm.name or f"frame_{mid}",
                frame_id=next_fid,
                body_id=body_id,
                local_position=local_pos,
                local_orientation=ori,
            )
            model.add_frame(iframe)
            self._marker_map[mid] = iframe.id
            next_fid += 1

    def _add_surfaces(self, raw, model):
        for sid, surf in raw.surfaces.items():
            rm_marker = self._marker_map.get(surf.rm_marker_id)
            self._surface_map[sid] = {
                'rm_marker': rm_marker,
                'num_patches': surf.num_patches,
                'num_nodes': surf.num_nodes,
                'patch_lines': surf.patch_lines,
                'body_id': None,
            }
            if rm_marker is not None:
                frame = model.frames.get(rm_marker)
                if frame is not None:
                    self._surface_map[sid]['body_id'] = frame.body_id

    def _add_contacts(self, raw, model):
        for cid, rc in raw.contacts.items():
            cp = ContactPair(
                name=rc.name or f"contact_{cid}",
                contact_id=cid,
                body_a_id=self._surface_body(rc.action_ggeom_id),
                body_b_id=self._surface_body(rc.base_ggeom_id),
                normal_stiffness=rc.stiffness / self.length_scale,
                friction_coefficient=rc.dynamic_friction,
                activation_distance=rc.boundary_penetration * self.length_scale,
                contact_mode='A_quad_B_sdf',
                quadrature_settings={'regularizer': 1e-4, 'damping': rc.damping},
            )
            model.add_contact(cp)
            self._contact_map[cid] = cp

    def _surface_body(self, ggeom_id):
        surf = self._surface_map.get(ggeom_id)
        if surf and surf['body_id'] is not None:
            return surf['body_id']
        return None

    def _add_joints(self, raw, model):
        next_jid = 0
        for jid, rj in raw.joints.items():
            jtype = JOINT_TYPE_MAP.get(rj.type)
            if jtype is None:
                continue
            fi = self._marker_map.get(rj.marker_i)
            fj = self._marker_map.get(rj.marker_j)
            if fi is None or fj is None:
                continue
            if jtype == 'FIXED':
                j = FixedJoint(rj.name or f"joint_{next_jid}", next_jid, fi, fj)
            elif jtype == 'REVOLUTE':
                j = RevoluteJoint(rj.name or f"joint_{next_jid}", next_jid, fi, fj, axis='z')
            else:
                continue
            model.add_joint(j)
            self._joint_map[rj.id] = j
            next_jid += 1

    # ============================================================
    # Fixed-body merging (参考意见v2 §11)
    # ============================================================
    def _merge_fixed_bodies(self, model):
        fixed_pairs = []
        for jid, joint in list(model.joints.items()):
            if joint.type.name != 'FIXED':
                continue
            fi = model.frames.get(joint.frame_i_id)
            fj = model.frames.get(joint.frame_j_id)
            if fi is None or fj is None:
                continue
            fixed_pairs.append((fi.body_id, fj.body_id, jid, fi.id, fj.id))

        if not fixed_pairs:
            return

        graph = {}
        for bi, bj, _, _, _ in fixed_pairs:
            graph.setdefault(bi, set()).add(bj)
            graph.setdefault(bj, set()).add(bi)

        visited = set()
        components = []
        for node in graph:
            if node in visited:
                continue
            comp = []
            stack = [node]
            while stack:
                n = stack.pop()
                if n in visited:
                    continue
                visited.add(n)
                comp.append(n)
                for nb in graph.get(n, set()):
                    if nb not in visited:
                        stack.append(nb)
            components.append(comp)

        merged_bodies = {}
        for comp in components:
            if len(comp) < 2:
                continue
            bodies = {bid: model.bodies[bid] for bid in comp if bid in model.bodies}
            if len(bodies) < 2:
                continue
            keep_id = min(bodies.keys())
            merge_ids = [bid for bid in bodies if bid != keep_id]

            merged_bodies[keep_id] = merge_ids
            self._do_merge(model, keep_id, merge_ids, fixed_pairs)

        if merged_bodies:
            self.warnings.append(
                f"Fixed-body merging: {len(merged_bodies)} compounds, "
                f"removed {sum(len(v) for v in merged_bodies.values())} bodies, "
                f"{sum(1 for _, _, j, _, _ in fixed_pairs)} fixed joints")

    def _do_merge(self, model, keep_id, merge_ids, fixed_pairs):
        keep_body = model.bodies[keep_id]

        # Collect world positions and inertias
        all_bodies = [keep_body] + [model.bodies[mid] for mid in merge_ids]

        r_world = np.array([b.initial_r for b in all_bodies])
        m = np.array([b.mass for b in all_bodies])
        m_total = np.sum(m)
        r_com = np.sum(m[:, None] * r_world, axis=0) / m_total

        keep_R = keep_body.initial_R
        I_world_total = np.zeros((3, 3))
        for i, b in enumerate(all_bodies):
            Iw = b.initial_R @ b.inertia_body @ b.initial_R.T
            d = r_world[i] - r_com
            d_outer = np.outer(d, d)
            Iw_steiner = Iw + m[i] * (np.dot(d, d) * np.eye(3) - d_outer)
            I_world_total += Iw_steiner

        I_local = keep_R.T @ I_world_total @ keep_R

        keep_body.mass = m_total
        keep_body.initial_r = r_com
        keep_body.inertia_body = I_local

        # Transfer frames from merged bodies to keep body
        for mid in merge_ids:
            for fid, frame in list(model.frames.items()):
                if frame.body_id == mid:
                    pw = keep_body.initial_r + keep_R @ frame.local_position
                    frame.local_position = keep_R.T @ (pw - r_com)
                    frame.local_orientation = keep_R.T @ (keep_body.initial_R @ frame.local_orientation)
                    frame.body_id = keep_id
            # Remove the merged body
            if mid in model.bodies:
                del model.bodies[mid]

        # Update contact body references
        for cp in model.contacts.values():
            if cp.body_a_id in merge_ids:
                cp.body_a_id = keep_id
            if cp.body_b_id in merge_ids:
                cp.body_b_id = keep_id

        # Update force body references
        for fe in model.forces.values():
            frame_i = model.frames.get(fe.frame_i_id)
            frame_j = model.frames.get(fe.frame_j_id)
            if frame_i and frame_i.body_id in merge_ids:
                frame_i.body_id = keep_id
            if frame_j and frame_j.body_id in merge_ids:
                frame_j.body_id = keep_id

        # Remove fixed joints that connected these bodies
        joints_to_remove = []
        for jid, joint in list(model.joints.items()):
            if joint.type.name != 'FIXED':
                continue
            fi = model.frames.get(joint.frame_i_id)
            fj = model.frames.get(joint.frame_j_id)
            if fi is None or fj is None:
                continue
            if (fi.body_id == keep_id and fj.body_id in merge_ids) or \
               (fj.body_id == keep_id and fi.body_id in merge_ids) or \
               (fi.body_id in merge_ids and fj.body_id in merge_ids):
                joints_to_remove.append(jid)

        for jid in joints_to_remove:
            del model.joints[jid]

        # Update remaining joint references
        for joint in model.joints.values():
            fi = model.frames.get(joint.frame_i_id)
            fj = model.frames.get(joint.frame_j_id)
            if fi and fi.body_id in merge_ids:
                pass
            if fj and fj.body_id in merge_ids:
                pass
        self._joint_map = {k: v for k, v in self._joint_map.items()
                           if v.id not in joints_to_remove}

        # Update surface body references
        for sid, surf in self._surface_map.items():
            if surf['body_id'] in merge_ids:
                surf['body_id'] = keep_id

    def _add_motions(self, raw, model):
        raw_to_ir = {}
        raw_to_ir.update(self._joint_map)

        next_did = 0
        for rm in raw.motions.values():
            ir_jid = raw_to_ir.get(rm.joint_id)
            if ir_jid is None:
                continue

            if rm.motion_type == 'velocity':
                func = self._parse_velocity_expression(rm.function_expr)
            else:
                func = self._parse_expression(rm.function_expr)

            if func is not None:
                fid = f"motion_fn_{rm.id}"
                model.add_function(fid, func)
                fn_name = fid
            else:
                fn_name = rm.function_expr

            drive = Drive(
                name=rm.name or f"drive_{next_did}",
                drive_id=next_did,
                target_joint_id=ir_jid.id,
                mode=DriveMode.MOTION,
                function_id=fn_name,
            )
            model.add_drive(drive)
            next_did += 1

    def _add_axial_forces(self, raw, model):
        next_fid = 0
        for af in raw.axial_forces:
            fi = self._marker_map.get(af['i'])
            fj = self._marker_map.get(af['j'])
            if fi is None or fj is None:
                continue

            func = self._parse_expression(af['function'])
            if func is None:
                func = ConstantFunction(0.0)

            force = ForceElement(
                name=af.get('name', f'force_{af["id"]}'),
                force_id=next_fid,
                frame_i_id=fi,
                frame_j_id=fj,
                magnitude_function=func,
                force_type='axial',
            )
            model.add_force(force)
            next_fid += 1

    def _parse_velocity_expression(self, expr):
        if not expr:
            return None
        e = expr.strip()
        if e == '2*PI':
            return LinearFunction(0.0, 2 * np.pi)
        try:
            val = float(e)
            return LinearFunction(0.0, val)
        except ValueError:
            pass
        if '*' in e and 'time' in e and 'PI' in e and '/' in e:
            parts = e.split('*')
            if len(parts) >= 2:
                try:
                    amp = float(parts[0].strip())
                    return LinearFunction(0.0, amp * np.pi / 180.0)
                except ValueError:
                    pass
        return None

    def _parse_expression(self, expr):
        if not expr:
            return None
        e = expr.strip()
        if e == '2*PI':
            return ConstantFunction(2 * np.pi)
        try:
            val = float(e)
            return ConstantFunction(val)
        except ValueError:
            pass
        if '*' in e and 'time' in e and 'PI' in e and '/' in e:
            parts = e.split('*')
            if len(parts) >= 2:
                try:
                    amp = float(parts[0].strip())
                    return LinearFunction(0.0, amp * np.pi / 180.0)
                except ValueError:
                    pass
        return None
