#pragma once

#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <vector>

namespace sdf {

inline constexpr int kMaxSdfQueryCandidates = 8;

struct SdfQueryCandidate
{
    double gap{0.0};
    Eigen::Vector3d gapNormal{1.0, 0.0, 0.0};
    Eigen::Vector3d witnessPoint{0.0, 0.0, 0.0};
    int branchId{-1};
};

struct SdfQueryResult
{
    double gap{0.0};
    Eigen::Vector3d gapNormal{1.0, 0.0, 0.0};
    Eigen::Vector3d witnessPoint{0.0, 0.0, 0.0};
    bool valid{false};
    std::vector<SdfQueryCandidate> candidates;
};

struct SdfFastQueryResult
{
    double gap{0.0};
    Eigen::Vector3d gapNormal{1.0, 0.0, 0.0};
    Eigen::Vector3d witnessPoint{0.0, 0.0, 0.0};
    bool valid{false};
    int candidateCount{0};
    std::array<SdfQueryCandidate, kMaxSdfQueryCandidates> candidates{};
};

} // namespace sdf
