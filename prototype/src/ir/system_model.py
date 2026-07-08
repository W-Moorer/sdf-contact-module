import numpy as np


class SystemModel:
    def __init__(self):
        self.gravity = np.array([0.0, 0.0, -9.81], dtype=float)
        self.bodies = {}
        self.frames = {}
        self.joints = {}
        self.drives = {}
        self.forces = {}
        self.torques = {}
        self.contacts = {}
        self.functions = {}
        self.solver_settings = {}

    def add_body(self, body):
        self.bodies[body.id] = body

    def add_frame(self, frame):
        self.frames[frame.id] = frame

    def add_joint(self, joint):
        self.joints[joint.id] = joint

    def add_drive(self, drive):
        self.drives[drive.id] = drive

    def add_force(self, force):
        self.forces[force.id] = force

    def add_torque(self, torque):
        self.torques[torque.id] = torque

    def add_contact(self, contact):
        self.contacts[contact.id] = contact

    def add_function(self, func_id, func):
        self.functions[func_id] = func

    @property
    def num_bodies(self):
        return len(self.bodies)

    @property
    def movable_bodies(self):
        return [b for b in self.bodies.values() if not b.fixed]

    @property
    def num_movable(self):
        return len(self.movable_bodies)
