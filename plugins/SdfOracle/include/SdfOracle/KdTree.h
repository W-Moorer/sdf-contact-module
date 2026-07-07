#pragma once

// ============================================================================
// Triangle BVH for nearest-triangle and ray-parity SDF queries.
// The public type name remains KdTree to keep the existing SdfBuilder API stable.
// ============================================================================

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace sdf {

struct ClosestPointResult
{
    Eigen::Vector3d point{0, 0, 0};
    int triangleIndex{-1};
    double distance{0.0};
    Eigen::Vector3d faceNormal{0, 0, 1};
};

inline Eigen::Vector3d closestPointOnTriangle(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c)
{
    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = p - a;
    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a;

    const Eigen::Vector3d bp = p - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3) return b;

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const double v = d1 / (d1 - d3);
        return a + v * ab;
    }

    const Eigen::Vector3d cp = p - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6) return c;

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const double w = d2 / (d2 - d6);
        return a + w * ac;
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b);
    }

    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    return a + ab * v + ac * w;
}

struct KdNode
{
    Eigen::Vector3d bmin{0, 0, 0};
    Eigen::Vector3d bmax{0, 0, 0};
    int left{-1};
    int right{-1};
    int triStart{-1};
    int triCount{0};
};

class KdTree
{
public:
    KdTree() = default;

    void build(
        const std::vector<Eigen::Vector3d>& vertices,
        const std::vector<std::array<int, 3>>& triangles)
    {
        m_vertices = &vertices;
        m_triangles = &triangles;
        const int n = static_cast<int>(triangles.size());
        m_indices.resize(n);
        std::iota(m_indices.begin(), m_indices.end(), 0);
        m_bounds.resize(n);

        for (int i = 0; i < n; ++i)
        {
            const auto& tri = triangles[i];
            const Eigen::Vector3d& a = vertices[tri[0]];
            const Eigen::Vector3d& b = vertices[tri[1]];
            const Eigen::Vector3d& c = vertices[tri[2]];
            TriangleBounds bounds;
            bounds.bmin = a.cwiseMin(b).cwiseMin(c);
            bounds.bmax = a.cwiseMax(b).cwiseMax(c);
            bounds.centroid = (a + b + c) / 3.0;
            m_bounds[i] = bounds;
        }

        m_nodes.clear();
        if (n > 0)
        {
            m_nodes.reserve(2 * n);
            buildRecursive(0, n);
        }
    }

    ClosestPointResult findNearest(const Eigen::Vector3d& p) const
    {
        ClosestPointResult best;
        best.distance = std::numeric_limits<double>::max();
        double bestSq = std::numeric_limits<double>::max();
        if (!m_nodes.empty())
        {
            searchNearest(0, p, best, bestSq);
        }
        return best;
    }

    struct SignedResult
    {
        double signedDistance{0.0};
        Eigen::Vector3d closestPoint{0, 0, 0};
        Eigen::Vector3d normal{0, 0, 1};
        int triangleIndex{-1};
    };

    SignedResult findSignedNearest(const Eigen::Vector3d& p) const
    {
        const ClosestPointResult cp = findNearest(p);
        SignedResult sr;
        sr.closestPoint = cp.point;
        sr.triangleIndex = cp.triangleIndex;
        sr.normal = cp.faceNormal;

        const int rayHits = countRayIntersections(p);
        if (rayHits > 0)
        {
            const bool inside = (rayHits % 2) == 1;
            sr.signedDistance = inside ? -cp.distance : cp.distance;
        }
        else
        {
            const Eigen::Vector3d diff = p - cp.point;
            sr.signedDistance = diff.dot(cp.faceNormal) >= 0.0
                ? cp.distance
                : -cp.distance;
        }
        return sr;
    }

    int triangleCount() const
    {
        return static_cast<int>(m_indices.size());
    }

private:
    struct TriangleBounds
    {
        Eigen::Vector3d bmin{0, 0, 0};
        Eigen::Vector3d bmax{0, 0, 0};
        Eigen::Vector3d centroid{0, 0, 0};
    };

    const std::vector<Eigen::Vector3d>* m_vertices{nullptr};
    const std::vector<std::array<int, 3>>* m_triangles{nullptr};
    std::vector<TriangleBounds> m_bounds;
    std::vector<int> m_indices;
    std::vector<KdNode> m_nodes;

    static constexpr int kLeafSize = 8;

    int buildRecursive(int lo, int hi)
    {
        const int nodeIdx = static_cast<int>(m_nodes.size());
        m_nodes.emplace_back();
        KdNode& node = m_nodes[nodeIdx];

        Eigen::Vector3d bmin(
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max());
        Eigen::Vector3d bmax(
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max());
        Eigen::Vector3d cmin = bmin;
        Eigen::Vector3d cmax = bmax;

        for (int i = lo; i < hi; ++i)
        {
            const TriangleBounds& bounds = m_bounds[m_indices[i]];
            bmin = bmin.cwiseMin(bounds.bmin);
            bmax = bmax.cwiseMax(bounds.bmax);
            cmin = cmin.cwiseMin(bounds.centroid);
            cmax = cmax.cwiseMax(bounds.centroid);
        }

        node.bmin = bmin;
        node.bmax = bmax;

        const int count = hi - lo;
        if (count <= kLeafSize)
        {
            node.triStart = lo;
            node.triCount = count;
            return nodeIdx;
        }

        const Eigen::Vector3d extent = cmax - cmin;
        int axis = 0;
        if (extent.y() > extent.x()) axis = 1;
        if (extent.z() > extent(axis)) axis = 2;

        const int mid = lo + count / 2;
        std::nth_element(
            m_indices.begin() + lo,
            m_indices.begin() + mid,
            m_indices.begin() + hi,
            [&](int a, int b) {
                return m_bounds[a].centroid(axis) < m_bounds[b].centroid(axis);
            });

        node.left = buildRecursive(lo, mid);
        node.right = buildRecursive(mid, hi);
        return nodeIdx;
    }

    static double squaredDistanceToAabb(
        const Eigen::Vector3d& p,
        const Eigen::Vector3d& bmin,
        const Eigen::Vector3d& bmax)
    {
        double sq = 0.0;
        for (int axis = 0; axis < 3; ++axis)
        {
            if (p(axis) < bmin(axis))
            {
                const double d = bmin(axis) - p(axis);
                sq += d * d;
            }
            else if (p(axis) > bmax(axis))
            {
                const double d = p(axis) - bmax(axis);
                sq += d * d;
            }
        }
        return sq;
    }

    void searchNearest(
        int nodeIdx,
        const Eigen::Vector3d& p,
        ClosestPointResult& best,
        double& bestSq) const
    {
        const KdNode& node = m_nodes[nodeIdx];
        if (squaredDistanceToAabb(p, node.bmin, node.bmax) > bestSq)
        {
            return;
        }

        if (node.triCount > 0)
        {
            for (int i = node.triStart; i < node.triStart + node.triCount; ++i)
            {
                const int triIdx = m_indices[i];
                const auto& tri = (*m_triangles)[triIdx];
                const Eigen::Vector3d& a = (*m_vertices)[tri[0]];
                const Eigen::Vector3d& b = (*m_vertices)[tri[1]];
                const Eigen::Vector3d& c = (*m_vertices)[tri[2]];

                const Eigen::Vector3d cp = closestPointOnTriangle(p, a, b, c);
                const double sq = (p - cp).squaredNorm();
                if (sq < bestSq)
                {
                    bestSq = sq;
                    best.point = cp;
                    best.triangleIndex = triIdx;
                    best.distance = std::sqrt(sq);
                    const Eigen::Vector3d n = (b - a).cross(c - a);
                    const double nLen = n.norm();
                    best.faceNormal =
                        nLen > 1.0e-15 ? n / nLen : Eigen::Vector3d{0.0, 0.0, 1.0};
                }
            }
            return;
        }

        const double leftSq = squaredDistanceToAabb(
            p, m_nodes[node.left].bmin, m_nodes[node.left].bmax);
        const double rightSq = squaredDistanceToAabb(
            p, m_nodes[node.right].bmin, m_nodes[node.right].bmax);

        int first = node.left;
        int second = node.right;
        double firstSq = leftSq;
        double secondSq = rightSq;
        if (rightSq < leftSq)
        {
            std::swap(first, second);
            std::swap(firstSq, secondSq);
        }

        if (firstSq <= bestSq)
        {
            searchNearest(first, p, best, bestSq);
        }
        if (secondSq <= bestSq)
        {
            searchNearest(second, p, best, bestSq);
        }
    }

    static Eigen::Vector3d rayDirection()
    {
        return Eigen::Vector3d{1.0, 0.3713906763541037, 0.52999894000318}.normalized();
    }

    static bool rayIntersectsAabb(
        const Eigen::Vector3d& origin,
        const Eigen::Vector3d& dir,
        const Eigen::Vector3d& bmin,
        const Eigen::Vector3d& bmax)
    {
        double tmin = 0.0;
        double tmax = std::numeric_limits<double>::max();
        for (int axis = 0; axis < 3; ++axis)
        {
            const double invD = 1.0 / dir(axis);
            double t0 = (bmin(axis) - origin(axis)) * invD;
            double t1 = (bmax(axis) - origin(axis)) * invD;
            if (t0 > t1)
            {
                std::swap(t0, t1);
            }
            tmin = std::max(tmin, t0);
            tmax = std::min(tmax, t1);
            if (tmax < tmin)
            {
                return false;
            }
        }
        return true;
    }

    static bool rayIntersectsTriangle(
        const Eigen::Vector3d& origin,
        const Eigen::Vector3d& dir,
        const Eigen::Vector3d& a,
        const Eigen::Vector3d& b,
        const Eigen::Vector3d& c)
    {
        constexpr double eps = 1.0e-12;
        const Eigen::Vector3d edge1 = b - a;
        const Eigen::Vector3d edge2 = c - a;
        const Eigen::Vector3d h = dir.cross(edge2);
        const double det = edge1.dot(h);
        if (std::abs(det) < eps)
        {
            return false;
        }

        const double invDet = 1.0 / det;
        const Eigen::Vector3d s = origin - a;
        const double u = invDet * s.dot(h);
        if (u < -eps || u > 1.0 + eps)
        {
            return false;
        }

        const Eigen::Vector3d q = s.cross(edge1);
        const double v = invDet * dir.dot(q);
        if (v < -eps || u + v > 1.0 + eps)
        {
            return false;
        }

        const double t = invDet * edge2.dot(q);
        return t > eps;
    }

    int countRayIntersections(const Eigen::Vector3d& origin) const
    {
        if (m_nodes.empty())
        {
            return 0;
        }
        const Eigen::Vector3d dir = rayDirection();
        return countRayIntersectionsRecursive(0, origin, dir);
    }

    int countRayIntersectionsRecursive(
        int nodeIdx,
        const Eigen::Vector3d& origin,
        const Eigen::Vector3d& dir) const
    {
        const KdNode& node = m_nodes[nodeIdx];
        if (!rayIntersectsAabb(origin, dir, node.bmin, node.bmax))
        {
            return 0;
        }

        if (node.triCount > 0)
        {
            int hits = 0;
            for (int i = node.triStart; i < node.triStart + node.triCount; ++i)
            {
                const int triIdx = m_indices[i];
                const auto& tri = (*m_triangles)[triIdx];
                if (rayIntersectsTriangle(
                        origin,
                        dir,
                        (*m_vertices)[tri[0]],
                        (*m_vertices)[tri[1]],
                        (*m_vertices)[tri[2]]))
                {
                    ++hits;
                }
            }
            return hits;
        }

        return countRayIntersectionsRecursive(node.left, origin, dir) +
            countRayIntersectionsRecursive(node.right, origin, dir);
    }
};

} // namespace sdf
