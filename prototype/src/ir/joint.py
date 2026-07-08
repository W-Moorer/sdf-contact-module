from enum import Enum


class JointType(Enum):
    FIXED = 'fixed'
    REVOLUTE = 'revolute'
    TRANSLATIONAL = 'translational'
    SPHERICAL = 'spherical'


class Joint:
    def __init__(self, name, joint_id, joint_type, frame_i_id, frame_j_id, parameters=None):
        self.name = name
        self.id = joint_id
        self.type = joint_type
        self.frame_i_id = frame_i_id
        self.frame_j_id = frame_j_id
        self.parameters = parameters or {}


class FixedJoint(Joint):
    def __init__(self, name, joint_id, frame_i_id, frame_j_id):
        super().__init__(name, joint_id, JointType.FIXED, frame_i_id, frame_j_id)


class RevoluteJoint(Joint):
    def __init__(self, name, joint_id, frame_i_id, frame_j_id, axis='z'):
        super().__init__(name, joint_id, JointType.REVOLUTE, frame_i_id, frame_j_id)
        self.axis = axis
        self.parameters['axis'] = axis
