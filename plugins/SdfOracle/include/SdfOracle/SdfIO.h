#pragma once

// ============================================================================
// SdfIO - Binary serialization for SdfData and SparseSdfData
// ============================================================================
// Custom binary format (.sdf):
//   Header (128 bytes):
//     magic[4]       "SDFO"
//     version        uint32 = 100 (semantic version 0.1.0)
//     nx, ny, nz     uint32 each
//     bmin[3]        double each
//     bmax[3]        double each
//     voxel[3]       double each
//     flags          bit 0: hasHessian, bit 1: sparse, bit 2: contact-aware
//     reserved       layout-specific fields, padded to 128 bytes
//   Layouts under version 0.1.0:
//     1 dense scalar
//     2 sparse scalar blocks
//     3 dense contact-aware
//     4 sparse contact-aware blocks
// ============================================================================

#include <string>

#include "SdfOracle/SdfData.h"

namespace sdf {

/// Save dense SdfData to binary file.
bool saveSdf(const SdfData& data, const std::string& path);

/// Load dense SdfData from binary file.
SdfData loadSdf(const std::string& path);

/// Save SparseSdfData to binary file.
bool saveSparseSdf(const SparseSdfData& data, const std::string& path);

/// Load SparseSdfData from binary file.
SparseSdfData loadSparseSdf(const std::string& path);

/// Check if a file is a valid SdfOracle data file with format version 0.1.0.
bool isValidSdf(const std::string& path);

/// Check if a SdfOracle 0.1.0 data file uses a sparse block layout.
bool isSparseSdf(const std::string& path);

} // namespace sdf
