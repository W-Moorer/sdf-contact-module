from .rmd_entities import (RawRecurDynModel, RawBody, RawMarker,
                           RawJoint, RawMotion, RawExpression,
                           RawSurface, RawContact)
from .rmd_lexer import RMDEntity


class RMDParser:
    def __init__(self):
        self.warnings = []

    def parse(self, entities):
        model = RawRecurDynModel()

        for ent in entities:
            kind = ent.kind
            try:
                if kind == 'PART':
                    self._handle_part(ent, model)
                elif kind == 'MARKER':
                    self._handle_marker(ent, model)
                elif kind == 'JOINT':
                    self._handle_joint(ent, model)
                elif kind == 'MOTION':
                    self._handle_motion(ent, model)
                elif kind == 'EXPRESSION':
                    self._handle_expression(ent, model)
                elif kind == 'GGEOM':
                    self._handle_ggeom(ent, model)
                elif kind == 'GGEOMCONTACT':
                    self._handle_contact(ent, model)
                elif kind == 'AXIAL_FORCE':
                    self._handle_axial_force(ent, model)
                elif kind == 'ACCGRAV':
                    self._handle_gravity(ent, model)
                elif kind == 'UNITS':
                    self._handle_units(ent, model)
                elif kind in ('OUTPUT', 'INTPAR', 'EQUILIBRIUM',
                              'SOLVEROPTION', 'STOPBYCONDITION', 'INFO'):
                    for k, v in ent.attrs.items():
                        model.solver_settings[k] = v
                else:
                    model.unknown_blocks.append(ent)
            except Exception as ex:
                self.warnings.append(f"{kind} / {ent.id}: {ex}")

        self._resolve_marker_bodies(model)
        return model

    def _handle_part(self, ent, model):
        body = RawBody(
            name=ent.get_str('NAME'),
            body_id=ent.id,
        )
        body.ground = 'GROUND' in [l.strip() for l in ent.raw_lines[0].split(',')]

        mass_str = ent.get('MASS')
        if mass_str is not None:
            body.mass = float(mass_str)

        cm_str = ent.get('CM')
        if cm_str is not None:
            body.cm_marker_id = int(float(cm_str))

        ip_vals = ent.get_floats('IP')
        if len(ip_vals) >= 6:
            body.inertia = [
                [ip_vals[0], ip_vals[3], ip_vals[4]],
                [ip_vals[3], ip_vals[1], ip_vals[5]],
                [ip_vals[4], ip_vals[5], ip_vals[2]],
            ]

        qg = ent.get_floats('QG')
        if len(qg) >= 3:
            body.initial_position = qg[:3]

        reuler = ent.get_floats('REULER')
        if len(reuler) >= 3:
            body.initial_rotation = reuler[:3]

        model.add_body(body)

    def _handle_marker(self, ent, model):
        m = RawMarker(
            name=ent.get_str('NAME'),
            marker_id=ent.id,
        )
        part_str = ent.get('PART')
        if part_str is not None:
            m.part_id = int(float(part_str))

        qp = ent.get_floats('QP')
        if len(qp) >= 3:
            m.position = qp[:3]

        reuler = ent.get_floats('REULER')
        if len(reuler) >= 3:
            m.orientation = reuler[:3]

        model.add_marker(m)

    def _handle_joint(self, ent, model):
        j = RawJoint(
            name=ent.get_str('NAME'),
            joint_id=ent.id,
        )

        i_str = ent.get('I')
        if i_str is not None:
            j.marker_i = int(float(i_str))
        j_str = ent.get('J')
        if j_str is not None:
            j.marker_j = int(float(j_str))

        for line in ent.raw_lines:
            s = line.strip().lstrip(', ')
            if s.upper() in ('FIXED', 'REVOLUTE', 'CYLINDRICAL', 'TRANSLATIONAL'):
                j.type = s.lower()
                break

        model.add_joint(j)

    def _handle_motion(self, ent, model):
        m = RawMotion(
            name=ent.get_str('NAME'),
            motion_id=ent.id,
        )
        joint_str = ent.get('JOINT')
        if joint_str is not None:
            m.joint_id = int(float(joint_str))

        m.is_rotation = ent.get('ROTATION') is not None
        if ent.get('VELOCITY') is not None:
            m.motion_type = 'velocity'
        elif ent.get('DISPLACEMENT') is not None:
            m.motion_type = 'displacement'
        else:
            m.motion_type = 'displacement'

        fn = ent.get('FUNCTION', '')
        m.function_expr = fn.rstrip('\\').strip()

        model.add_motion(m)

    def _handle_expression(self, ent, model):
        e = RawExpression(
            name=ent.get_str('NAME'),
            expr_id=ent.id,
        )
        fn = ent.get('FUNCTION', '')
        e.function_expr = fn.rstrip('\\').strip()
        model.add_expression(e)

    def _handle_ggeom(self, ent, model):
        surf = RawSurface(surface_id=ent.id)
        surf.name = ent.get_str('NAME')
        rm_str = ent.get('RM')
        if rm_str is not None:
            surf.rm_marker_id = int(float(rm_str))
        surf.num_patches = ent.get_int('NO_PATCH')
        surf.num_nodes = ent.get_int('NO_NODE')

        patching = False
        for line in ent.raw_lines:
            if line.strip().startswith('PATCHES'):
                patching = True
                continue
            if patching:
                surf.patch_lines.append(line)

        model.add_surface(surf)

    def _handle_contact(self, ent, model):
        c = RawContact(contact_id=ent.id)
        c.name = ent.get_str('NAME')
        c.action_ggeom_id = ent.get_int('IGGEOMID')
        c.base_ggeom_id = ent.get_int('JGGEOMID')
        c.action_marker_id = ent.get_int('IFLOAT')
        c.base_marker_id = ent.get_int('JFLOAT')
        c.stiffness = ent.get_float('K', 100000.0)
        c.damping = ent.get_float('C', 10.0)
        c.dynamic_friction = ent.get_float('D_F_C', 0.0)
        c.static_friction = ent.get_float('S_F_C', 0.0)
        c.static_transition_vel = ent.get_float('S_T_V', 0.1)
        c.boundary_penetration = ent.get_float('BPEN', 0.01)
        model.add_contact(c)

    def _handle_axial_force(self, ent, model):
        model.axial_forces.append({
            'id': ent.id,
            'name': ent.get_str('NAME'),
            'i': ent.get_int('I'),
            'j': ent.get_int('J'),
            'function': ent.get('FUNCTION', '').rstrip('\\').strip(),
            'is_translation': ent.get('TRANSLATION') is not None,
        })

    def _handle_gravity(self, ent, model):
        vals = ent.get_floats('KGRAV')
        if not vals:
            vals = [0.0]
        model.gravity = [0.0, 0.0, vals[0]]

    def _handle_units(self, ent, model):
        for k, v in ent.attrs.items():
            model.units[k] = v.strip("' ")

    def _resolve_marker_bodies(self, model):
        for m in model.markers.values():
            if m.part_id is not None and m.part_id in model.bodies:
                pass
