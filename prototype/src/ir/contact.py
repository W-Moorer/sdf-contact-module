class ContactPair:
    def __init__(self, name, contact_id, body_a_id, body_b_id,
                 normal_stiffness=1e6, friction_coefficient=0.5,
                 activation_distance=0.0, contact_mode='A_quad_B_sdf',
                 quadrature_settings=None, k_order=1,
                 action_marker_id=None, base_marker_id=None):
        self.name = name
        self.id = contact_id
        self.body_a_id = body_a_id
        self.body_b_id = body_b_id
        self.normal_stiffness = normal_stiffness
        self.friction_coefficient = friction_coefficient
        self.activation_distance = activation_distance
        self.contact_mode = contact_mode
        self.quadrature_settings = quadrature_settings or {}
        self.k_order = k_order
        self.action_marker_id = action_marker_id
        self.base_marker_id = base_marker_id
