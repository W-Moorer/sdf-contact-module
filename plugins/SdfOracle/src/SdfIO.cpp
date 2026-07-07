// ============================================================================
// SdfIO.cpp - Binary serialization (.sdf format)
// ============================================================================

#include "SdfOracle/SdfIO.h"

#include <cstring>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace sdf {

namespace {
constexpr uint32_t encodeSdfVersion(uint32_t major, uint32_t minor, uint32_t patch)
{
    return major * 10000u + minor * 100u + patch;
}

constexpr uint32_t kSdfVersionMajor = 0;
constexpr uint32_t kSdfVersionMinor = 1;
constexpr uint32_t kSdfVersionPatch = 0;
constexpr uint32_t kSdfVersion =
    encodeSdfVersion(kSdfVersionMajor, kSdfVersionMinor, kSdfVersionPatch);
constexpr const char* kSdfVersionString = "0.1.0";
constexpr uint32_t kFlagHasHessian = 1u << 0;
constexpr uint32_t kFlagSparseBlocks = 1u << 1;
constexpr uint32_t kFlagContactAwareCandidates = 1u << 2;
constexpr uint32_t kDenseLayout = 1u;
constexpr uint32_t kSparseBlocksLayout = 2u;
constexpr uint32_t kContactAwareDenseLayout = 3u;
constexpr uint32_t kContactAwareSparseLayout = 4u;

void putU32(char* bytes, int offset, uint32_t value)
{
    std::memcpy(bytes + offset, &value, sizeof(value));
}

uint32_t getU32(const char* bytes, int offset)
{
    uint32_t value = 0;
    std::memcpy(&value, bytes + offset, sizeof(value));
    return value;
}

void putDouble(char* bytes, int offset, double value)
{
    std::memcpy(bytes + offset, &value, sizeof(value));
}

double getDouble(const char* bytes, int offset)
{
    double value = 0.0;
    std::memcpy(&value, bytes + offset, sizeof(value));
    return value;
}
} // namespace

// ---------- Binary format header (128 bytes) ----------

#pragma pack(push, 1)
struct SdfHeader
{
    char magic[4]{'S', 'D', 'F', 'O'};
    uint32_t version{kSdfVersion};
    uint32_t nx{0}, ny{0}, nz{0};
    double bmin[3]{0, 0, 0};
    double bmax[3]{0, 0, 0};
    double voxel[3]{0, 0, 0};
    uint32_t flags{1}; // bit 0: hasHessian
    char reserved[32]{}; // pad to 128 bytes total (96 + 32 = 128)

    SdfHeader()
    {
        static_assert(sizeof(SdfHeader) == 128,
                      "SdfHeader must be exactly 128 bytes");
        std::memset(reserved, 0, sizeof(reserved));
    }
};

struct SparseBlockHeader
{
    int32_t bi{0};
    int32_t bj{0};
    int32_t bk{0};
    uint32_t reserved{0};
};
#pragma pack(pop)

// ---------- Save ----------

bool saveSdf(const SdfData& data, const std::string& path)
{
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.is_open())
        return false;

    const auto& spec = data.spec;
    int N = spec.totalVoxels();
    const bool hasContactCandidates = data.hasContactCandidates();

    // Build header
    SdfHeader hdr;
    hdr.version = kSdfVersion;
    hdr.nx = static_cast<uint32_t>(spec.nx());
    hdr.ny = static_cast<uint32_t>(spec.ny());
    hdr.nz = static_cast<uint32_t>(spec.nz());
    hdr.bmin[0] = spec.bmin.x(); hdr.bmin[1] = spec.bmin.y(); hdr.bmin[2] = spec.bmin.z();
    hdr.bmax[0] = spec.bmax.x(); hdr.bmax[1] = spec.bmax.y(); hdr.bmax[2] = spec.bmax.z();
    hdr.voxel[0] = spec.voxel.x(); hdr.voxel[1] = spec.voxel.y(); hdr.voxel[2] = spec.voxel.z();
    hdr.flags = kFlagHasHessian |
        (hasContactCandidates ? kFlagContactAwareCandidates : 0u);
    putU32(
        hdr.reserved,
        0,
        hasContactCandidates ? kContactAwareDenseLayout : kDenseLayout);
    if (hasContactCandidates)
    {
        putU32(
            hdr.reserved,
            4,
            static_cast<uint32_t>(data.contactCandidateStride));
    }

    ofs.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

    // Write phi0[N]
    ofs.write(reinterpret_cast<const char*>(data.phi0.data()),
              N * sizeof(double));

    // Write n0[N] as interleaved doubles (x,y,z)
    for (int idx = 0; idx < N; ++idx)
    {
        const auto& n = data.n0[idx];
        double buf[3] = {n.x(), n.y(), n.z()};
        ofs.write(reinterpret_cast<const char*>(buf), 3 * sizeof(double));
    }

    // Write H0[N] as 6 doubles (Hxx, Hxy, Hxz, Hyy, Hyz, Hzz)
    for (int idx = 0; idx < N; ++idx)
    {
        const auto& H = data.H0[idx];
        double buf[6] = {H(0,0), H(0,1), H(0,2), H(1,1), H(1,2), H(2,2)};
        ofs.write(reinterpret_cast<const char*>(buf), 6 * sizeof(double));
    }

    // Write witness0[N] as interleaved doubles
    for (int idx = 0; idx < N; ++idx)
    {
        const auto& w = data.witness0[idx];
        double buf[3] = {w.x(), w.y(), w.z()};
        ofs.write(reinterpret_cast<const char*>(buf), 3 * sizeof(double));
    }

    if (hasContactCandidates)
    {
        for (int idx = 0; idx < N; ++idx)
        {
            const uint32_t count =
                static_cast<uint32_t>(data.contactCandidateCount[static_cast<size_t>(idx)]);
            ofs.write(reinterpret_cast<const char*>(&count), sizeof(count));
        }

        ofs.write(
            reinterpret_cast<const char*>(data.contactCandidatePhi0.data()),
            static_cast<std::streamsize>(
                data.contactCandidatePhi0.size() * sizeof(double)));

        for (const Eigen::Vector3d& normal : data.contactCandidateNormal0)
        {
            double buf[3] = {normal.x(), normal.y(), normal.z()};
            ofs.write(reinterpret_cast<const char*>(buf), 3 * sizeof(double));
        }

        ofs.write(
            reinterpret_cast<const char*>(data.contactCandidateBranchId.data()),
            static_cast<std::streamsize>(
                data.contactCandidateBranchId.size() * sizeof(int32_t)));
    }

    return ofs.good();
}

bool saveSparseSdf(const SparseSdfData& data, const std::string& path)
{
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.is_open())
        return false;

    const auto& spec = data.spec;
    const int blockVoxelCount = data.blockVoxelCount();
    if (data.blockDim <= 0 || blockVoxelCount <= 0)
        return false;

    const bool hasContactCandidates = data.hasContactCandidates();

    SdfHeader hdr;
    hdr.version = kSdfVersion;
    hdr.nx = static_cast<uint32_t>(spec.nx());
    hdr.ny = static_cast<uint32_t>(spec.ny());
    hdr.nz = static_cast<uint32_t>(spec.nz());
    hdr.bmin[0] = spec.bmin.x(); hdr.bmin[1] = spec.bmin.y(); hdr.bmin[2] = spec.bmin.z();
    hdr.bmax[0] = spec.bmax.x(); hdr.bmax[1] = spec.bmax.y(); hdr.bmax[2] = spec.bmax.z();
    hdr.voxel[0] = spec.voxel.x(); hdr.voxel[1] = spec.voxel.y(); hdr.voxel[2] = spec.voxel.z();
    hdr.flags = kFlagSparseBlocks | (data.hasHessian ? kFlagHasHessian : 0u) |
        (hasContactCandidates ? kFlagContactAwareCandidates : 0u);
    putU32(
        hdr.reserved,
        0,
        hasContactCandidates ? kContactAwareSparseLayout : kSparseBlocksLayout);
    putU32(hdr.reserved, 4, static_cast<uint32_t>(data.blockDim));
    putU32(hdr.reserved, 8, static_cast<uint32_t>(data.haloVoxels));
    putU32(hdr.reserved, 12, static_cast<uint32_t>(data.blocks.size()));
    putDouble(hdr.reserved, 16, data.bandDistance);
    if (hasContactCandidates)
    {
        putU32(
            hdr.reserved,
            24,
            static_cast<uint32_t>(data.contactCandidateStride));
    }

    ofs.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

    for (const SparseSdfBlock& block : data.blocks)
    {
        if (block.phi0.size() < static_cast<size_t>(blockVoxelCount) ||
            block.normal0.size() < static_cast<size_t>(3 * blockVoxelCount) ||
            (data.hasHessian &&
                block.hessian0.size() < static_cast<size_t>(6 * blockVoxelCount)) ||
            (hasContactCandidates &&
                !block.hasContactCandidates(
                    blockVoxelCount,
                    data.contactCandidateStride)))
        {
            return false;
        }

        SparseBlockHeader blockHeader;
        blockHeader.bi = static_cast<int32_t>(block.index[0]);
        blockHeader.bj = static_cast<int32_t>(block.index[1]);
        blockHeader.bk = static_cast<int32_t>(block.index[2]);
        ofs.write(reinterpret_cast<const char*>(&blockHeader), sizeof(blockHeader));
        ofs.write(
            reinterpret_cast<const char*>(block.phi0.data()),
            blockVoxelCount * sizeof(float));
        ofs.write(
            reinterpret_cast<const char*>(block.normal0.data()),
            3 * blockVoxelCount * sizeof(float));
        if (data.hasHessian)
        {
            ofs.write(
                reinterpret_cast<const char*>(block.hessian0.data()),
                6 * blockVoxelCount * sizeof(float));
        }
        if (hasContactCandidates)
        {
            ofs.write(
                reinterpret_cast<const char*>(block.contactCandidateCount.data()),
                blockVoxelCount * sizeof(uint8_t));
            ofs.write(
                reinterpret_cast<const char*>(block.contactCandidatePhi0.data()),
                blockVoxelCount * data.contactCandidateStride * sizeof(float));
            ofs.write(
                reinterpret_cast<const char*>(block.contactCandidateNormal0.data()),
                3 * blockVoxelCount * data.contactCandidateStride * sizeof(float));
            ofs.write(
                reinterpret_cast<const char*>(block.contactCandidateBranchId.data()),
                blockVoxelCount * data.contactCandidateStride * sizeof(int32_t));
        }
    }

    return ofs.good();
}

// ---------- Load ----------

SdfData loadSdf(const std::string& path)
{
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.is_open())
        throw std::runtime_error("Cannot open SdfOracle data file: " + path);

    SdfHeader hdr;
    ifs.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
    if (ifs.gcount() != sizeof(hdr))
        throw std::runtime_error("Failed to read SdfOracle data header");

    if (std::memcmp(hdr.magic, "SDFO", 4) != 0)
        throw std::runtime_error("Invalid SdfOracle data magic number");
    if (hdr.version != kSdfVersion)
        throw std::runtime_error("Unsupported SdfOracle data version: " +
                                 std::to_string(hdr.version) +
                                 " (only " + std::string(kSdfVersionString) +
                                 " is supported)");
    if ((hdr.flags & kFlagSparseBlocks) != 0)
        throw std::runtime_error("SdfOracle data file uses sparse layout; call loadSparseSdf");
    const uint32_t layout = getU32(hdr.reserved, 0);
    if (layout != kDenseLayout && layout != kContactAwareDenseLayout)
        throw std::runtime_error("SdfOracle data file is not dense 0.1.0 layout");
    const bool contactFlag = (hdr.flags & kFlagContactAwareCandidates) != 0;
    if ((layout == kContactAwareDenseLayout) != contactFlag)
        throw std::runtime_error("Dense SDF 0.1.0 contact-aware flag/layout mismatch");
    const bool hasContactCandidates =
        contactFlag && layout == kContactAwareDenseLayout;

    GridSpec spec;
    spec.shape = {static_cast<int>(hdr.nx),
                  static_cast<int>(hdr.ny),
                  static_cast<int>(hdr.nz)};
    spec.bmin = Eigen::Vector3d(hdr.bmin[0], hdr.bmin[1], hdr.bmin[2]);
    spec.bmax = Eigen::Vector3d(hdr.bmax[0], hdr.bmax[1], hdr.bmax[2]);
    spec.voxel = Eigen::Vector3d(hdr.voxel[0], hdr.voxel[1], hdr.voxel[2]);

    SdfData data;
    data.resize(spec);
    int N = spec.totalVoxels();

    // Read phi0
    ifs.read(reinterpret_cast<char*>(data.phi0.data()), N * sizeof(double));

    // Read n0
    for (int idx = 0; idx < N; ++idx)
    {
        double buf[3];
        ifs.read(reinterpret_cast<char*>(buf), 3 * sizeof(double));
        data.n0[idx] = Eigen::Vector3d(buf[0], buf[1], buf[2]);
    }

    // Read H0
    bool hasHessian = (hdr.flags & 1) != 0;
    if (hasHessian)
    {
        for (int idx = 0; idx < N; ++idx)
        {
            double buf[6];
            ifs.read(reinterpret_cast<char*>(buf), 6 * sizeof(double));
            Eigen::Matrix3d H;
            H(0,0) = buf[0]; H(0,1) = buf[1]; H(0,2) = buf[2];
            H(1,0) = buf[1]; H(1,1) = buf[3]; H(1,2) = buf[4];
            H(2,0) = buf[2]; H(2,1) = buf[4]; H(2,2) = buf[5];
            data.H0[idx] = H;
        }
    }

    // Read witness0
    for (int idx = 0; idx < N; ++idx)
    {
        double buf[3];
        ifs.read(reinterpret_cast<char*>(buf), 3 * sizeof(double));
        data.witness0[idx] = Eigen::Vector3d(buf[0], buf[1], buf[2]);
    }

    if (hasContactCandidates)
    {
        const int stride = static_cast<int>(getU32(hdr.reserved, 4));
        if (stride <= 0 || stride > SdfData::kMaxContactCandidates)
            throw std::runtime_error("Invalid contact-aware candidate stride");

        data.resizeContactCandidates(stride);
        for (int idx = 0; idx < N; ++idx)
        {
            uint32_t count = 0;
            ifs.read(reinterpret_cast<char*>(&count), sizeof(count));
            if (count > static_cast<uint32_t>(stride))
                throw std::runtime_error("Invalid contact-aware candidate count");
            data.contactCandidateCount[static_cast<size_t>(idx)] =
                static_cast<uint8_t>(count);
        }

        ifs.read(
            reinterpret_cast<char*>(data.contactCandidatePhi0.data()),
            static_cast<std::streamsize>(
                data.contactCandidatePhi0.size() * sizeof(double)));

        for (Eigen::Vector3d& normal : data.contactCandidateNormal0)
        {
            double buf[3];
            ifs.read(reinterpret_cast<char*>(buf), 3 * sizeof(double));
            normal = Eigen::Vector3d(buf[0], buf[1], buf[2]);
        }

        ifs.read(
            reinterpret_cast<char*>(data.contactCandidateBranchId.data()),
            static_cast<std::streamsize>(
                data.contactCandidateBranchId.size() * sizeof(int32_t)));
    }

    if (!ifs)
        throw std::runtime_error("Failed to read dense SDF payload");

    return data;
}

SparseSdfData loadSparseSdf(const std::string& path)
{
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.is_open())
        throw std::runtime_error("Cannot open SdfOracle data file: " + path);

    SdfHeader hdr;
    ifs.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
    if (ifs.gcount() != sizeof(hdr))
        throw std::runtime_error("Failed to read SdfOracle data header");

    if (std::memcmp(hdr.magic, "SDFO", 4) != 0)
        throw std::runtime_error("Invalid SdfOracle data magic number");
    if (hdr.version != kSdfVersion)
        throw std::runtime_error("Unsupported SdfOracle data version: " +
                                 std::to_string(hdr.version) +
                                 " (only " + std::string(kSdfVersionString) +
                                 " is supported)");
    const uint32_t layout = getU32(hdr.reserved, 0);
    const bool contactFlag = (hdr.flags & kFlagContactAwareCandidates) != 0;
    if ((layout == kContactAwareSparseLayout) != contactFlag)
        throw std::runtime_error("Sparse SDF 0.1.0 contact-aware flag/layout mismatch");
    const bool hasContactCandidates =
        contactFlag && layout == kContactAwareSparseLayout;

    if ((hdr.flags & kFlagSparseBlocks) == 0 ||
        !(layout == kSparseBlocksLayout || hasContactCandidates))
    {
        throw std::runtime_error("SdfOracle data file is not sparse block 0.1.0");
    }

    SparseSdfData data;
    data.spec.shape = {static_cast<int>(hdr.nx),
                       static_cast<int>(hdr.ny),
                       static_cast<int>(hdr.nz)};
    data.spec.bmin = Eigen::Vector3d(hdr.bmin[0], hdr.bmin[1], hdr.bmin[2]);
    data.spec.bmax = Eigen::Vector3d(hdr.bmax[0], hdr.bmax[1], hdr.bmax[2]);
    data.spec.voxel = Eigen::Vector3d(hdr.voxel[0], hdr.voxel[1], hdr.voxel[2]);
    data.blockDim = static_cast<int>(getU32(hdr.reserved, 4));
    data.haloVoxels = static_cast<int>(getU32(hdr.reserved, 8));
    const uint32_t blockCount = getU32(hdr.reserved, 12);
    data.bandDistance = getDouble(hdr.reserved, 16);
    data.hasHessian = (hdr.flags & kFlagHasHessian) != 0;
    data.contactCandidateStride =
        hasContactCandidates ? static_cast<int>(getU32(hdr.reserved, 24)) : 0;

    if (data.blockDim <= 0)
        throw std::runtime_error("Sparse SDF blockDim must be positive");
    if (hasContactCandidates &&
        (data.contactCandidateStride <= 0 ||
         data.contactCandidateStride > SdfData::kMaxContactCandidates))
    {
        throw std::runtime_error("Invalid sparse contact-aware candidate stride");
    }

    const int blockVoxelCount = data.blockVoxelCount();
    data.blocks.reserve(blockCount);
    for (uint32_t b = 0; b < blockCount; ++b)
    {
        SparseBlockHeader blockHeader;
        ifs.read(reinterpret_cast<char*>(&blockHeader), sizeof(blockHeader));
        if (!ifs)
            throw std::runtime_error("Failed to read sparse SDF block header");

        SparseSdfBlock block;
        block.index = {blockHeader.bi, blockHeader.bj, blockHeader.bk};
        block.phi0.resize(blockVoxelCount);
        block.normal0.resize(3 * blockVoxelCount);
        if (data.hasHessian)
        {
            block.hessian0.resize(6 * blockVoxelCount);
        }
        if (hasContactCandidates)
        {
            block.contactCandidateCount.resize(blockVoxelCount);
            block.contactCandidatePhi0.resize(
                static_cast<size_t>(blockVoxelCount * data.contactCandidateStride));
            block.contactCandidateNormal0.resize(
                static_cast<size_t>(3 * blockVoxelCount * data.contactCandidateStride));
            block.contactCandidateBranchId.resize(
                static_cast<size_t>(blockVoxelCount * data.contactCandidateStride));
        }

        ifs.read(
            reinterpret_cast<char*>(block.phi0.data()),
            blockVoxelCount * sizeof(float));
        ifs.read(
            reinterpret_cast<char*>(block.normal0.data()),
            3 * blockVoxelCount * sizeof(float));
        if (data.hasHessian)
        {
            ifs.read(
                reinterpret_cast<char*>(block.hessian0.data()),
                6 * blockVoxelCount * sizeof(float));
        }
        if (hasContactCandidates)
        {
            ifs.read(
                reinterpret_cast<char*>(block.contactCandidateCount.data()),
                blockVoxelCount * sizeof(uint8_t));
            ifs.read(
                reinterpret_cast<char*>(block.contactCandidatePhi0.data()),
                blockVoxelCount * data.contactCandidateStride * sizeof(float));
            ifs.read(
                reinterpret_cast<char*>(block.contactCandidateNormal0.data()),
                3 * blockVoxelCount * data.contactCandidateStride * sizeof(float));
            ifs.read(
                reinterpret_cast<char*>(block.contactCandidateBranchId.data()),
                blockVoxelCount * data.contactCandidateStride * sizeof(int32_t));

            for (uint8_t count : block.contactCandidateCount)
            {
                if (count > data.contactCandidateStride)
                {
                    throw std::runtime_error(
                        "Invalid sparse contact-aware candidate count");
                }
            }
        }
        if (!ifs)
            throw std::runtime_error("Failed to read sparse SDF block payload");

        data.blocks.emplace_back(std::move(block));
    }

    data.rebuildIndex();
    return data;
}

// ---------- Validation ----------

bool isValidSdf(const std::string& path)
{
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.is_open()) return false;

    SdfHeader hdr;
    ifs.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
    if (ifs.gcount() != sizeof(hdr)) return false;
    return std::memcmp(hdr.magic, "SDFO", 4) == 0 &&
        hdr.version == kSdfVersion;
}

bool isSparseSdf(const std::string& path)
{
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.is_open()) return false;

    SdfHeader hdr;
    ifs.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
    if (ifs.gcount() != sizeof(hdr)) return false;
    if (std::memcmp(hdr.magic, "SDFO", 4) != 0) return false;
    if (hdr.version != kSdfVersion) return false;
    const uint32_t layout = getU32(hdr.reserved, 0);
    return (hdr.flags & kFlagSparseBlocks) != 0 &&
        (layout == kSparseBlocksLayout || layout == kContactAwareSparseLayout);
}

} // namespace sdf
