import numpy as np


class Frame:
    def __init__(self, name, frame_id, body_id,
                 local_position=None, local_orientation=None):
        self.name = name
        self.id = frame_id
        self.body_id = body_id
        self.local_position = np.zeros(3) if local_position is None else np.asarray(local_position, dtype=float)
        self.local_orientation = np.eye(3) if local_orientation is None else np.asarray(local_orientation, dtype=float)
