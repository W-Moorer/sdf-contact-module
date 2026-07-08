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
        # If query fully encloses root AABB, return all points
        if (query_min <= self.root.aabb_min - margin).all() and (query_max >= self.root.aabb_max + margin).all():
            return np.arange(len(self.points), dtype=np.int32)
        stack = [self.root]
        parts = []
        while stack:
            node = stack.pop()
            if np.any(node.aabb_max + margin < query_min) or np.any(node.aabb_min - margin > query_max):
                continue
            if node.indices is not None:
                parts.append(node.indices)
            else:
                stack.append(node.left)
                stack.append(node.right)
        if not parts:
            return np.zeros(0, dtype=np.int32)
        return np.concatenate(parts)

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
