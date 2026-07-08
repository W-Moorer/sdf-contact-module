from enum import Enum


class DriveMode(Enum):
    MOTION = 'motion'
    VELOCITY = 'velocity'
    FORCE = 'force'
    TORQUE = 'torque'


class Drive:
    def __init__(self, name, drive_id, target_joint_id, mode=DriveMode.MOTION,
                 function_id=None, coordinate='angle'):
        self.name = name
        self.id = drive_id
        self.target_joint_id = target_joint_id
        self.mode = mode
        self.function_id = function_id
        self.coordinate = coordinate
