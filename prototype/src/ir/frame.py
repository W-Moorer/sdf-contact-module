import numpy as np


class Frame:
    def __init__(self, name, frame_id, body_id,
                 local_position=None, local_orientation=None):
        self.name = name
        self.id = frame_id
        self.body_id = body_id
        self.local_position = np.zeros(3) if local_position is None else np.asarray(local_position, dtype=float)
        self.local_orientation = np.eye(3) if local_orientation is None else np.asarray(local_orientation, dtype=float)

    def world_pose(self, state):
        """Compute this frame's origin position and orientation in world coordinates."""
        idx = state.body_idx(self.body_id)
        r_body = state.r[idx]
        R_body = state.R[idx]
        return r_body + R_body @ self.local_position, R_body @ self.local_orientation

    def relative_to(self, other, state):
        """Return (R, t) such that p_other = R @ p_self + t.
        
        Transforms a point from this frame's local coordinates to
        the other frame's local coordinates, using current body states.
        """
        r_self, R_self = self.world_pose(state)
        r_other, R_other = other.world_pose(state)
        R = R_other.T @ R_self
        t = R_other.T @ (r_self - r_other)
        return R, t
