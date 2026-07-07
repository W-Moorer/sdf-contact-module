#pragma once

// ============================================================================
// NexDynAdapter - integration note
// ============================================================================
// The production NexDyn adapter is implemented in:
//   src/Contact/SDF/SdfOracleAdapter.h
//
// This plugin header intentionally stays independent of NexDyn core headers so
// SdfOracle remains buildable as a standalone library.
// ============================================================================

#include "SdfOracle/SdfOracle.h"

namespace sdf {

using StandaloneSdfOracle = SdfOracle;

} // namespace sdf
