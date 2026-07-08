class RawBody:
    def __init__(self, name='', body_id=0):
        self.name = name
        self.id = body_id
        self.mass = 1.0
        self.inertia = None
        self.cm_marker_id = None
        self.initial_position = None
        self.initial_rotation = None
        self.initial_velocity = None
        self.initial_angular_velocity = None
        self.ground = False

    def __repr__(self):
        return f"RawBody(name='{self.name}', id={self.id}, mass={self.mass})"


class RawMarker:
    def __init__(self, name='', marker_id=0):
        self.name = name
        self.id = marker_id
        self.part_id = None
        self.position = None
        self.orientation = None

    def __repr__(self):
        return f"RawMarker(name='{self.name}', id={self.id}, part={self.part_id})"


class RawJoint:
    def __init__(self, name='', joint_id=0):
        self.name = name
        self.id = joint_id
        self.type = ''
        self.marker_i = None
        self.marker_j = None

    def __repr__(self):
        return f"RawJoint(name='{self.name}', type={self.type})"


class RawMotion:
    def __init__(self, name='', motion_id=0):
        self.name = name
        self.id = motion_id
        self.joint_id = None
        self.motion_type = 'displacement'
        self.is_rotation = True
        self.function_expr = ''

    def __repr__(self):
        return f"RawMotion(name='{self.name}', joint={self.joint_id})"


class RawExpression:
    def __init__(self, name='', expr_id=0):
        self.name = name
        self.id = expr_id
        self.function_expr = ''

    def __repr__(self):
        return f"RawExpression(name='{self.name}')"


class RawSurface:
    def __init__(self, surface_id=0):
        self.id = surface_id
        self.name = ''
        self.rm_marker_id = None
        self.body_id = None
        self.num_patches = 0
        self.num_nodes = 0
        self.patch_lines = []

    def __repr__(self):
        return f"RawSurface(id={self.id}, patches={self.num_patches})"


class RawContact:
    def __init__(self, contact_id=0):
        self.id = contact_id
        self.name = ''
        self.action_ggeom_id = None
        self.base_ggeom_id = None
        self.action_marker_id = None
        self.base_marker_id = None
        self.stiffness = 100000.0
        self.damping = 10.0
        self.dynamic_friction = 0.0
        self.static_friction = 0.0
        self.static_transition_vel = 0.1
        self.boundary_penetration = 0.01

    def __repr__(self):
        return f"RawContact(id={self.id}, G={self.action_ggeom_id}->{self.base_ggeom_id})"


class RawRecurDynModel:
    def __init__(self):
        self.units = {}
        self.gravity = None
        self.bodies = {}
        self.markers = {}
        self.joints = {}
        self.motions = {}
        self.expressions = {}
        self.surfaces = {}
        self.contacts = {}
        self.axial_forces = []
        self.solver_settings = {}
        self.unknown_blocks = []
        self.part_by_name = {}
        self.marker_by_name = {}

    def add_body(self, body):
        self.bodies[body.id] = body
        self.part_by_name[body.name] = body

    def add_marker(self, marker):
        self.markers[marker.id] = marker
        self.marker_by_name[marker.name] = marker

    def add_joint(self, joint):
        self.joints[joint.id] = joint

    def add_motion(self, motion):
        self.motions[motion.id] = motion

    def add_expression(self, expr):
        self.expressions[expr.id] = expr

    def add_surface(self, surf):
        self.surfaces[surf.id] = surf

    def add_contact(self, contact):
        self.contacts[contact.id] = contact
