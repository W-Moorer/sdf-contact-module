#pragma once

#include <Eigen/Dense>
#include <array>
#include <string>
#include <vector>

namespace sdf {

struct TriangleMesh
{
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<int, 3>> triangles;
};

struct CornerNormalTriangleMesh
{
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::vector<std::array<Eigen::Vector3d, 3>> cornerNormals;
};

TriangleMesh loadObjFile(const std::string& path);
CornerNormalTriangleMesh loadContactAwareSurfaceFile(const std::string& path);

} // namespace sdf
