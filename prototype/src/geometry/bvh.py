"""Simple AABB-BVH for accelerating contact broad-phase culling."""
import numpy as np


class BVHNode:
    __slots__ = ('aabb_min', 'aabb_max', 'left', 'right', 'indices')

    def __init__(self, aabb_min, aabb_max, left=None, right=None, indices=None):
        self.aabb_min = aabb_min
        self.aabb_max = aabb_max
        self.left = left
        self.right = right
        self.indices = indices


class AABB_BVH:
    """Axis-Aligned Bounding Box hierarchy over a set of 3D points."""

    def __init__(self, points, leaf_size=16):
        """
        points: (N, 3) array of point coordinates in body-local frame.
        leaf_size: max points per leaf node.
        """
        self.points = np.asarray(points, dtype=np.float64)
        self.leaf_size = leaf_size
        self.root = self._build(np.arange(len(self.points)), 0)

    def _build(self, indices, depth):
        n = len(indices)
        pts = self.points[indices]
        aabb_min = pts.min(axis=0)
        aabb_max = pts.max(axis=0)

        if n <= self.leaf_size:
            return BVHNode(aabb_min, aabb_max, indices=indices.copy())

        # Split along the longest axis
        extents = aabb_max - aabb_min
        axis = int(np.argmax(extents))
        mid = (aabb_min[axis] + aabb_max[axis]) * 0.5

        left_idx = indices[pts[:, axis] <= mid]
        right_idx = indices[pts[:, axis] > mid]

        # Avoid degenerate splits (all points on one side)
        if len(left_idx) == 0 or len(right_idx) == 0:
            # Use median split instead
            sorted_idx = np.argsort(pts[:, axis])
            half = n // 2
            left_idx = indices[sorted_idx[:half]]
            right_idx = indices[sorted_idx[half:]]

        left = self._build(left_idx, depth + 1)
        right = self._build(right_idx, depth + 1)
        return BVHNode(aabb_min, aabb_max, left=left, right=right)

    def query(self, query_min, query_max, margin=0.0):
        """Return indices of points whose AABB overlaps with query AABB + margin."""
        result = []
        self._query_node(self.root, query_min, query_max, margin, result)
        return np.array(result, dtype=np.int32)

    def _query_node(self, node, qmin, qmax, margin, result):
        # Test overlap: node AABB vs query AABB
        if np.any(node.aabb_max + margin < qmin) or np.any(node.aabb_min - margin > qmax):
            return

        if node.indices is not None:  # Leaf
            result.extend(node.indices.tolist())
        else:  # Internal
            self._query_node(node.left, qmin, qmax, margin, result)
            self._query_node(node.right, qmin, qmax, margin, result)

    def query_sphere(self, center, radius):
        """Return indices of points within radius of center."""
        result = []
        self._query_sphere_node(self.root, center, radius, result)
        return np.array(result, dtype=np.int32)

    def _query_sphere_node(self, node, center, radius, result):
        # Test sphere-AABB overlap
        closest = np.maximum(node.aabb_min, np.minimum(node.aabb_max, center))
        dist_sq = np.sum((closest - center) ** 2)
        if dist_sq > radius * radius:
            return

        if node.indices is not None:
            result.extend(node.indices.tolist())
        else:
            self._query_sphere_node(node.left, center, radius, result)
            self._query_sphere_node(node.right, center, radius, result)
