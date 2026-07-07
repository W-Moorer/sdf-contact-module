// ============================================================================
// SdfSurfaceMesh.cpp -- OBJ and contact-aware surface input loaders
// ============================================================================

#include "SdfOracle/SdfSurfaceMesh.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace sdf {
namespace {

std::string stripCommentAndTrim(std::string line)
{
    const size_t comment = line.find('#');
    if (comment != std::string::npos)
    {
        line.erase(comment);
    }
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
    {
        return {};
    }
    const auto last = line.find_last_not_of(" \t\r\n");
    return line.substr(first, last - first + 1);
}

std::runtime_error parseError(
    const std::string& path,
    int lineNumber,
    const std::string& message)
{
    return std::runtime_error(
        path + ":" + std::to_string(lineNumber) + ": " + message);
}

Eigen::Vector3d normalizedOrThrow(
    const std::string& path,
    int lineNumber,
    const Eigen::Vector3d& value,
    const std::string& label)
{
    const double norm = value.norm();
    if (!value.allFinite() || !std::isfinite(norm) || norm <= 1.0e-15)
    {
        throw parseError(path, lineNumber, label + " must be a finite non-zero vector");
    }
    return value / norm;
}

struct ContactAwareSurfaceParseState
{
    bool active{false};
    std::string name;
    int expectedVertices{-1};
    int expectedTriangles{-1};
    int seenVertices{0};
    int seenTriangles{0};
    Eigen::Matrix4d transform{Eigen::Matrix4d::Identity()};
    Eigen::Matrix3d normalTransform{Eigen::Matrix3d::Identity()};
    std::unordered_map<int, int> vertexIds;
};

void closeSurface(
    const std::string& path,
    int lineNumber,
    ContactAwareSurfaceParseState& state)
{
    if (!state.active)
    {
        return;
    }
    if (state.expectedVertices < 0)
    {
        throw parseError(path, lineNumber, "surface '" + state.name + "' is missing vertices block");
    }
    if (state.expectedTriangles < 0)
    {
        throw parseError(path, lineNumber, "surface '" + state.name + "' is missing triangles block");
    }
    if (state.seenVertices != state.expectedVertices)
    {
        throw parseError(path, lineNumber, "surface '" + state.name + "' vertex count mismatch");
    }
    if (state.seenTriangles != state.expectedTriangles)
    {
        throw parseError(path, lineNumber, "surface '" + state.name + "' triangle count mismatch");
    }
    state = ContactAwareSurfaceParseState{};
}

Eigen::Vector3d transformPoint(
    const Eigen::Matrix4d& transform,
    const Eigen::Vector3d& point)
{
    const Eigen::Vector4d homogeneous(point.x(), point.y(), point.z(), 1.0);
    const Eigen::Vector4d out = transform * homogeneous;
    if (std::abs(out.w()) > 1.0e-15)
    {
        return out.head<3>() / out.w();
    }
    return out.head<3>();
}

} // namespace

TriangleMesh loadObjFile(const std::string& path)
{
    TriangleMesh mesh;
    std::ifstream ifs(path);
    if (!ifs.is_open())
    {
        throw std::runtime_error("Cannot open OBJ file: " + path);
    }

    std::string line;
    while (std::getline(ifs, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "v")
        {
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            iss >> x >> y >> z;
            mesh.vertices.emplace_back(x, y, z);
        }
        else if (token == "f")
        {
            std::array<int, 3> tri{};
            for (int i = 0; i < 3; ++i)
            {
                std::string tok;
                iss >> tok;
                // Handle v, v/vt, v/vt/vn, v//vn.
                const int idx = std::stoi(tok.substr(0, tok.find('/'))) - 1;
                tri[i] = idx;
            }
            mesh.triangles.push_back(tri);
        }
    }
    return mesh;
}

CornerNormalTriangleMesh loadContactAwareSurfaceFile(const std::string& path)
{
    std::ifstream ifs(path);
    if (!ifs.is_open())
    {
        throw std::runtime_error("Cannot open contact-aware surface file: " + path);
    }

    CornerNormalTriangleMesh mesh;
    ContactAwareSurfaceParseState state;
    bool headerSeen = false;
    std::string line;
    int lineNumber = 0;

    while (std::getline(ifs, line))
    {
        ++lineNumber;
        line = stripCommentAndTrim(line);
        if (line.empty())
        {
            continue;
        }

        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (!headerSeen)
        {
            int version = 0;
            if (token != "NEXDYN_CONTACT_AWARE_SURFACE" || !(iss >> version) || version != 1)
            {
                throw parseError(path, lineNumber, "expected header 'NEXDYN_CONTACT_AWARE_SURFACE 1'");
            }
            headerSeen = true;
            continue;
        }

        if (token == "surface")
        {
            closeSurface(path, lineNumber, state);
            state.active = true;
            if (!(iss >> state.name) || state.name.empty())
            {
                throw parseError(path, lineNumber, "surface requires a name");
            }
            continue;
        }

        if (token == "end")
        {
            closeSurface(path, lineNumber, state);
            continue;
        }

        if (!state.active)
        {
            throw parseError(path, lineNumber, "expected surface block");
        }

        if (token == "transform")
        {
            for (int row = 0; row < 4; ++row)
            {
                for (int col = 0; col < 4; ++col)
                {
                    if (!(iss >> state.transform(row, col)))
                    {
                        throw parseError(path, lineNumber, "transform requires 16 row-major values");
                    }
                }
            }
            Eigen::Matrix3d linear = state.transform.block<3, 3>(0, 0);
            if (std::abs(linear.determinant()) <= 1.0e-15)
            {
                throw parseError(path, lineNumber, "transform linear part must be invertible");
            }
            state.normalTransform = linear.inverse().transpose();
            continue;
        }

        if (token == "vertices")
        {
            if (!(iss >> state.expectedVertices) || state.expectedVertices < 0)
            {
                throw parseError(path, lineNumber, "vertices requires a non-negative count");
            }
            continue;
        }

        if (token == "triangles")
        {
            if (!(iss >> state.expectedTriangles) || state.expectedTriangles < 0)
            {
                throw parseError(path, lineNumber, "triangles requires a non-negative count");
            }
            continue;
        }

        if (token == "v")
        {
            if (state.expectedVertices < 0)
            {
                throw parseError(path, lineNumber, "vertex appears before vertices block");
            }
            int id = 0;
            Eigen::Vector3d point;
            if (!(iss >> id >> point.x() >> point.y() >> point.z()))
            {
                throw parseError(path, lineNumber, "vertex requires: v <id> <x> <y> <z>");
            }
            if (state.vertexIds.find(id) != state.vertexIds.end())
            {
                throw parseError(path, lineNumber, "duplicate vertex id in surface");
            }
            const int globalIndex = static_cast<int>(mesh.vertices.size());
            state.vertexIds.emplace(id, globalIndex);
            mesh.vertices.push_back(transformPoint(state.transform, point));
            ++state.seenVertices;
            continue;
        }

        if (token == "tri")
        {
            if (state.expectedTriangles < 0)
            {
                throw parseError(path, lineNumber, "triangle appears before triangles block");
            }

            int triId = 0;
            int vertexId[3]{};
            Eigen::Vector3d normals[3];
            if (!(iss >> triId >> vertexId[0] >> vertexId[1] >> vertexId[2] >>
                  normals[0].x() >> normals[0].y() >> normals[0].z() >>
                  normals[1].x() >> normals[1].y() >> normals[1].z() >>
                  normals[2].x() >> normals[2].y() >> normals[2].z()))
            {
                throw parseError(
                    path,
                    lineNumber,
                    "triangle requires id, 3 vertex ids, and 3 corner normals");
            }

            std::array<int, 3> tri{};
            for (int i = 0; i < 3; ++i)
            {
                const auto it = state.vertexIds.find(vertexId[i]);
                if (it == state.vertexIds.end())
                {
                    throw parseError(path, lineNumber, "triangle references unknown vertex id");
                }
                tri[i] = it->second;
                normals[i] = normalizedOrThrow(
                    path,
                    lineNumber,
                    state.normalTransform * normals[i],
                    "corner normal");
            }

            const Eigen::Vector3d edge0 = mesh.vertices[tri[1]] - mesh.vertices[tri[0]];
            const Eigen::Vector3d edge1 = mesh.vertices[tri[2]] - mesh.vertices[tri[0]];
            if (edge0.cross(edge1).norm() <= 1.0e-15)
            {
                throw parseError(path, lineNumber, "triangle is degenerate");
            }

            mesh.triangles.push_back(tri);
            mesh.cornerNormals.push_back(std::array<Eigen::Vector3d, 3>{
                normals[0],
                normals[1],
                normals[2]});
            ++state.seenTriangles;
            continue;
        }

        throw parseError(path, lineNumber, "unknown token '" + token + "'");
    }

    if (!headerSeen)
    {
        throw std::runtime_error(path + ": missing contact-aware surface header");
    }
    closeSurface(path, lineNumber, state);
    if (mesh.vertices.empty() || mesh.triangles.empty())
    {
        throw std::runtime_error(path + ": contact-aware surface contains no triangles");
    }
    return mesh;
}

} // namespace sdf
