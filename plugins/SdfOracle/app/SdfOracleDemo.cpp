// ============================================================================
// SdfOracleDemo.cpp - Standalone validation executable
// ============================================================================
// Tests the complete SdfOracle pipeline:
//   1. Load OBJ mesh
//   2. Build SDF grid (φ₀, n₀, H₀)
//   3. Save/load .sdf file
//   4. Query points in all 3 modes (1st order, 2nd order, trilinear)
//   5. Report results and timing
//
// Usage:
//   SdfOracleDemo --obj <mesh.obj> [--resolution 128] [--out <prefix>]
//   SdfOracleDemo --load <file.sdf> [--query <x> <y> <z>]
//   SdfOracleDemo --auto-test [--data-dir <obj_library_path>]
// ============================================================================

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "SdfOracle/SdfBuilder.h"
#include "SdfOracle/SdfIO.h"
#include "SdfOracle/SdfOracle.h"

using namespace sdf;

// ---------- Predefined test meshes & resolutions ----------

static const std::vector<std::string> kTestMeshes = {
    "cube.obj",       // 344v/684t
    "sphere.obj",     // 762v/1520t
    "hemisphere.obj", // 2051v/4000t
    "RoughCube.obj",  // 3694v/7384t
};

static const std::vector<int> kResolutions = {64, 128, 256};

// ---------- Helper: parse CLI ----------

struct DemoConfig
{
    std::string objPath;
    std::string contactAwareModelPath;
    std::string loadPath;
    std::string outPrefix{"sdf_output"};
    std::string dataDir;          // for --auto-test
    int resolution{128};
    double padding{0.1};
    bool sparse{false};
    int blockDim{8};
    int tileDim{8};
    int haloVoxels{2};
    int bandVoxels{8};
    double bandDistance{0.0};
    bool storeHessian{true};
    int maxCandidates{4};
    double featureTolerance{0.0};
    double normalMergeAngleDegrees{5.0};
    std::vector<Eigen::Vector3d> queryPoints;
    bool runValidation{true};
    bool autoTest{false};
    bool benchmark{false};
    bool compareDenseSparse{false};
    int queryCount{1000000};
    std::string csvPath;
};

DemoConfig parseArgs(int argc, char* argv[])
{
    DemoConfig cfg;
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--obj" && i + 1 < argc)
            cfg.objPath = argv[++i];
        else if (arg == "--contact-aware-model" && i + 1 < argc)
            cfg.contactAwareModelPath = argv[++i];
        else if (arg == "--load" && i + 1 < argc)
            cfg.loadPath = argv[++i];
        else if (arg == "--out" && i + 1 < argc)
            cfg.outPrefix = argv[++i];
        else if (arg == "--resolution" && i + 1 < argc)
            cfg.resolution = std::atoi(argv[++i]);
        else if (arg == "--padding" && i + 1 < argc)
            cfg.padding = std::atof(argv[++i]);
        else if (arg == "--sparse")
            cfg.sparse = true;
        else if (arg == "--contact-aware-sparse")
            cfg.sparse = true;
        else if (arg == "--block-dim" && i + 1 < argc)
            cfg.blockDim = std::atoi(argv[++i]);
        else if (arg == "--tile-dim" && i + 1 < argc)
            cfg.tileDim = std::atoi(argv[++i]);
        else if (arg == "--halo-voxels" && i + 1 < argc)
            cfg.haloVoxels = std::atoi(argv[++i]);
        else if (arg == "--band-voxels" && i + 1 < argc)
            cfg.bandVoxels = std::atoi(argv[++i]);
        else if (arg == "--band-distance" && i + 1 < argc)
            cfg.bandDistance = std::atof(argv[++i]);
        else if (arg == "--max-candidates" && i + 1 < argc)
            cfg.maxCandidates = std::atoi(argv[++i]);
        else if (arg == "--feature-tolerance" && i + 1 < argc)
            cfg.featureTolerance = std::atof(argv[++i]);
        else if (arg == "--normal-merge-angle" && i + 1 < argc)
            cfg.normalMergeAngleDegrees = std::atof(argv[++i]);
        else if (arg == "--no-hessian")
            cfg.storeHessian = false;
        else if (arg == "--data-dir" && i + 1 < argc)
            cfg.dataDir = argv[++i];
        else if (arg == "--query" && i + 3 < argc)
        {
            double x = std::atof(argv[++i]);
            double y = std::atof(argv[++i]);
            double z = std::atof(argv[++i]);
            cfg.queryPoints.emplace_back(x, y, z);
        }
        else if (arg == "--no-validation")
            cfg.runValidation = false;
        else if (arg == "--auto-test")
            cfg.autoTest = true;
        else if (arg == "--benchmark")
            cfg.benchmark = true;
        else if (arg == "--compare-dense-sparse")
            cfg.compareDenseSparse = true;
        else if (arg == "--query-count" && i + 1 < argc)
            cfg.queryCount = std::atoi(argv[++i]);
        else if (arg == "--csv" && i + 1 < argc)
            cfg.csvPath = argv[++i];
    }
    return cfg;
}

// ---------- Helper: print query result ----------

void printResult(const std::string& label, const SdfQueryResult& r)
{
    std::cout << "  [" << label << "]"
              << " gap=" << std::setw(10) << std::fixed << std::setprecision(6) << r.gap
              << "  normal=(" << std::setprecision(4)
              << r.gapNormal.x() << ", " << r.gapNormal.y() << ", " << r.gapNormal.z() << ")"
              << "  valid=" << (r.valid ? "yes" : "no")
              << "\n";
}

// ---------- Helper: generate test query points ----------

std::vector<Eigen::Vector3d> generateTestPoints(const GridSpec& spec)
{
    std::vector<Eigen::Vector3d> pts;
    Eigen::Vector3d center = (spec.bmin + spec.bmax) * 0.5;
    Eigen::Vector3d range = (spec.bmax - spec.bmin) * 0.3;

    // Grid of test points around center
    for (double fx = -1.0; fx <= 1.0; fx += 0.5)
    for (double fy = -1.0; fy <= 1.0; fy += 0.5)
    for (double fz = -1.0; fz <= 1.0; fz += 0.5)
    {
        pts.emplace_back(
            center.x() + fx * range.x(),
            center.y() + fy * range.y(),
            center.z() + fz * range.z());
    }
    return pts;
}

std::string queryModeName(QueryMode mode)
{
    switch (mode)
    {
    case QueryMode::FirstOrder:
        return "FirstOrder";
    case QueryMode::SecondOrder:
        return "SecondOrder";
    case QueryMode::Trilinear:
        return "Trilinear";
    case QueryMode::Tricubic:
        return "Tricubic";
    case QueryMode::ContactAware:
        return "ContactAware";
    default:
        return "Unknown";
    }
}

uint64_t fileSizeOrZero(const std::string& path)
{
    if (path.empty())
        return 0;
    std::error_code error;
    const auto size = std::filesystem::file_size(path, error);
    return error ? 0 : static_cast<uint64_t>(size);
}

double nextUnit(uint64_t& state)
{
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return static_cast<double>((state >> 11) & ((1ull << 53) - 1)) /
        static_cast<double>(1ull << 53);
}

std::vector<Eigen::Vector3d> generateBenchmarkPoints(
    const SdfOracle& oracle,
    int count)
{
    count = std::max(1, count);
    std::vector<Eigen::Vector3d> points;
    points.reserve(static_cast<size_t>(count));
    const GridSpec& spec = oracle.gridSpec();

    if (oracle.isSparse() && !oracle.sparseData().blocks.empty())
    {
        const SparseSdfData& sparse = oracle.sparseData();
        const int blockDim = sparse.blockDim;
        const int blockVoxelCount = sparse.blockVoxelCount();
        for (int n = 0; n < count; ++n)
        {
            const SparseSdfBlock& block =
                sparse.blocks[static_cast<size_t>(n) % sparse.blocks.size()];
            int local = static_cast<int>(
                (static_cast<uint64_t>(n) * 2654435761ull) %
                static_cast<uint64_t>(blockVoxelCount));
            int lx = local % blockDim;
            int ly = (local / blockDim) % blockDim;
            int lz = local / (blockDim * blockDim);
            if (blockDim > 2)
            {
                lx = 1 + (lx % (blockDim - 2));
                ly = 1 + (ly % (blockDim - 2));
                lz = 1 + (lz % (blockDim - 2));
            }
            const int i = block.index[0] * blockDim + lx;
            const int j = block.index[1] * blockDim + ly;
            const int k = block.index[2] * blockDim + lz;
            if (!spec.inBounds(i, j, k))
            {
                points.push_back((spec.bmin + spec.bmax) * 0.5);
                continue;
            }
            const double ox = (static_cast<double>((n * 13) % 9) - 4.0) / 40.0;
            const double oy = (static_cast<double>((n * 17) % 9) - 4.0) / 40.0;
            const double oz = (static_cast<double>((n * 19) % 9) - 4.0) / 40.0;
            points.push_back(
                spec.centerFromIndex(i, j, k) +
                spec.voxel.cwiseProduct(Eigen::Vector3d{ox, oy, oz}));
        }
        return points;
    }

    uint64_t state = 0x9e3779b97f4a7c15ull;
    const double fxMin = spec.nx() > 3 ? 0.05 : 0.0;
    const double fyMin = spec.ny() > 3 ? 0.05 : 0.0;
    const double fzMin = spec.nz() > 3 ? 0.05 : 0.0;
    const double fxMax = spec.nx() > 3 ? static_cast<double>(spec.nx()) - 1.05 : 0.0;
    const double fyMax = spec.ny() > 3 ? static_cast<double>(spec.ny()) - 1.05 : 0.0;
    const double fzMax = spec.nz() > 3 ? static_cast<double>(spec.nz()) - 1.05 : 0.0;
    for (int n = 0; n < count; ++n)
    {
        const double fx = fxMin + (fxMax - fxMin) * nextUnit(state);
        const double fy = fyMin + (fyMax - fyMin) * nextUnit(state);
        const double fz = fzMin + (fzMax - fzMin) * nextUnit(state);
        points.push_back(
            spec.bmin +
            spec.voxel.cwiseProduct(Eigen::Vector3d{fx + 0.5, fy + 0.5, fz + 0.5}));
    }
    return points;
}

struct QueryBenchmarkRow
{
    std::string sourcePath;
    std::string layout;
    std::string mode;
    int queryCount{0};
    int validCount{0};
    double avgAbsGap{0.0};
    double queryMs{0.0};
    double queriesPerSecond{0.0};
    int nx{0};
    int ny{0};
    int nz{0};
    size_t activeBlocks{0};
    size_t storedVoxels{0};
    uint64_t fileBytes{0};
    double buildMs{0.0};
    double loadMs{0.0};
};

void appendBenchmarkCsv(
    const std::string& path,
    const std::vector<QueryBenchmarkRow>& rows)
{
    if (path.empty() || rows.empty())
    {
        return;
    }
    const std::filesystem::path csvPath(path);
    if (csvPath.has_parent_path())
    {
        std::error_code error;
        std::filesystem::create_directories(csvPath.parent_path(), error);
    }
    const bool writeHeader = !std::filesystem::exists(path) ||
        std::filesystem::file_size(path) == 0;
    std::ofstream out(path, std::ios::app);
    if (!out.is_open())
    {
        std::cerr << "[Benchmark] Failed to open CSV: " << path << "\n";
        return;
    }
    if (writeHeader)
    {
        out << "source,layout,mode,query_count,valid_count,valid_rate,"
            << "avg_abs_gap,query_ms,queries_per_second,nx,ny,nz,"
            << "active_blocks,stored_voxels,file_bytes,build_ms,load_ms\n";
    }
    for (const QueryBenchmarkRow& row : rows)
    {
        const double validRate = row.queryCount > 0
            ? static_cast<double>(row.validCount) / row.queryCount
            : 0.0;
        out << row.sourcePath << ','
            << row.layout << ','
            << row.mode << ','
            << row.queryCount << ','
            << row.validCount << ','
            << validRate << ','
            << row.avgAbsGap << ','
            << row.queryMs << ','
            << row.queriesPerSecond << ','
            << row.nx << ','
            << row.ny << ','
            << row.nz << ','
            << row.activeBlocks << ','
            << row.storedVoxels << ','
            << row.fileBytes << ','
            << row.buildMs << ','
            << row.loadMs << '\n';
    }
}

std::vector<QueryBenchmarkRow> runQueryBenchmark(
    const SdfOracle& oracle,
    int queryCount,
    const std::string& sourcePath,
    const std::string& sdfPath,
    double buildMs,
    double loadMs)
{
    std::cout << "\n=== Benchmark: Query Throughput ===\n";
    std::cout << "  Query points: " << queryCount << "\n";
    const std::vector<Eigen::Vector3d> points =
        generateBenchmarkPoints(oracle, queryCount);
    const GridSpec& spec = oracle.gridSpec();
    const std::vector<QueryMode> modes{
        QueryMode::FirstOrder,
        QueryMode::SecondOrder,
        QueryMode::Trilinear,
        QueryMode::Tricubic,
        QueryMode::ContactAware};

    std::vector<QueryBenchmarkRow> rows;
    rows.reserve(modes.size());
    for (QueryMode mode : modes)
    {
        const auto t0 = std::chrono::high_resolution_clock::now();
        int valid = 0;
        double sumAbsGap = 0.0;
#pragma omp parallel for reduction(+:valid,sumAbsGap) schedule(static)
        for (int i = 0; i < static_cast<int>(points.size()); ++i)
        {
            const SdfFastQueryResult result =
                oracle.queryFast(points[static_cast<size_t>(i)], mode);
            if (result.valid)
            {
                ++valid;
                sumAbsGap += std::abs(result.gap);
            }
        }
        const auto t1 = std::chrono::high_resolution_clock::now();
        const double ms =
            std::chrono::duration<double, std::milli>(t1 - t0).count();

        QueryBenchmarkRow row;
        row.sourcePath = sourcePath;
        row.layout = oracle.isSparse() ? "sparse_blocks" : "dense";
        row.mode = queryModeName(mode);
        row.queryCount = static_cast<int>(points.size());
        row.validCount = valid;
        row.avgAbsGap = valid > 0 ? sumAbsGap / valid : 0.0;
        row.queryMs = ms;
        row.queriesPerSecond = ms > 0.0
            ? 1000.0 * static_cast<double>(points.size()) / ms
            : 0.0;
        row.nx = spec.nx();
        row.ny = spec.ny();
        row.nz = spec.nz();
        if (oracle.isSparse())
        {
            row.activeBlocks = oracle.sparseData().blocks.size();
            row.storedVoxels = oracle.sparseData().totalStoredVoxels();
        }
        else
        {
            row.storedVoxels = static_cast<size_t>(spec.totalVoxels());
        }
        row.fileBytes = fileSizeOrZero(sdfPath);
        row.buildMs = buildMs;
        row.loadMs = loadMs;
        rows.push_back(row);

        std::cout << "  " << std::setw(14) << std::left << row.mode
                  << " | valid=" << row.validCount << "/" << row.queryCount
                  << " | avg|gap|=" << std::setprecision(6) << row.avgAbsGap
                  << " | " << std::setprecision(3) << row.queryMs << " ms"
                  << " | " << std::setprecision(3)
                  << (row.queriesPerSecond / 1.0e6) << " Mq/s\n";
    }
    return rows;
}

void compareDenseSparseIfRequested(
    const DemoConfig& cfg,
    const TriangleMesh& mesh,
    const SdfOracle& sparseOracle)
{
    if (!cfg.compareDenseSparse || !sparseOracle.isSparse())
    {
        return;
    }
    const int totalVoxels = sparseOracle.gridSpec().totalVoxels();
    constexpr int kMaxDenseCompareVoxels = 5000000;
    if (totalVoxels > kMaxDenseCompareVoxels)
    {
        std::cout << "\n[Compare] Skipping dense/sparse comparison: "
                  << totalVoxels << " dense voxels exceeds safe limit "
                  << kMaxDenseCompareVoxels << "\n";
        return;
    }

    BuilderParams denseParams;
    denseParams.targetResolution = cfg.resolution;
    denseParams.maxResolution = cfg.resolution;
    denseParams.padding = cfg.padding;
    std::cout << "\n[Compare] Building dense reference for sparse comparison...\n";
    SdfData denseData = SdfBuilder::build(mesh, denseParams);
    SdfOracle denseOracle(std::move(denseData));
    const std::vector<Eigen::Vector3d> points =
        generateBenchmarkPoints(sparseOracle, std::min(cfg.queryCount, 20000));

    for (QueryMode mode : {
        QueryMode::FirstOrder,
        QueryMode::SecondOrder,
        QueryMode::Trilinear,
        QueryMode::Tricubic})
    {
        double sum = 0.0;
        double maxError = 0.0;
        int validBoth = 0;
        for (const Eigen::Vector3d& point : points)
        {
            const SdfQueryResult dense = denseOracle.query(point, mode);
            const SdfQueryResult sparse = sparseOracle.query(point, mode);
            if (dense.valid && sparse.valid)
            {
                const double error = std::abs(dense.gap - sparse.gap);
                sum += error;
                maxError = std::max(maxError, error);
                ++validBoth;
            }
        }
        std::cout << "  " << std::setw(14) << std::left << queryModeName(mode)
                  << " | valid_both=" << validBoth << "/" << points.size()
                  << " | mean_abs=" << std::setprecision(6)
                  << (validBoth > 0 ? sum / validBoth : 0.0)
                  << " | max_abs=" << maxError << "\n";
    }
}

// ---------- Validation: compare query modes ----------

void runValidation(const SdfOracle& oracle)
{
    std::cout << "\n=== Validation: Query Mode Comparison ===\n";

    const auto& spec = oracle.gridSpec();
    auto testPts = generateTestPoints(spec);

    std::cout << "  Testing " << testPts.size() << " points across 4 query modes\n\n";

    // Time each mode
    auto timeQuery = [&](QueryMode mode, const std::string& name)
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto results = oracle.queryBatch(testPts, mode);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        int valid = 0;
        double sumGap = 0.0;
        for (const auto& r : results)
        {
            if (r.valid) { ++valid; sumGap += std::abs(r.gap); }
        }
        std::cout << "  " << std::setw(14) << std::left << name
                  << " | " << valid << "/" << results.size() << " valid"
                  << " | avg|gap|=" << std::setprecision(6)
                  << (valid > 0 ? sumGap / valid : 0.0)
                  << " | " << std::setprecision(2) << ms << " ms\n";
    };

    timeQuery(QueryMode::FirstOrder, "1st-order");
    timeQuery(QueryMode::SecondOrder, "2nd-order");
    timeQuery(QueryMode::Trilinear, "Trilinear");
    timeQuery(QueryMode::Tricubic, "Tricubic");
    if ((oracle.isSparse() && oracle.sparseData().hasContactCandidates()) ||
        (!oracle.isSparse() && oracle.data().hasContactCandidates()))
    {
        timeQuery(QueryMode::ContactAware, "ContactAware");
    }

    // Compare modes for a few points
    std::cout << "\n  Sample comparisons (first 5 valid points):\n";
    int shown = 0;
    for (size_t i = 0; i < testPts.size() && shown < 5; ++i)
    {
        auto r1 = oracle.query(testPts[i], QueryMode::FirstOrder);
        auto r2 = oracle.query(testPts[i], QueryMode::SecondOrder);
        auto r3 = oracle.query(testPts[i], QueryMode::Trilinear);
        auto r4 = oracle.query(testPts[i], QueryMode::Tricubic);
        if (r1.valid && r2.valid)
        {
            std::cout << "    p=(" << std::setprecision(3)
                      << testPts[i].x() << ", "
                      << testPts[i].y() << ", "
                      << testPts[i].z() << ")"
                      << " | 2nd=" << std::setprecision(6) << r2.gap
                      << " | tri=" << r3.gap
                      << " | cubic=" << r4.gap
                      << "\n";
            ++shown;
        }
    }
}

// ---------- Auto-test: multi-mesh, multi-resolution ----------

struct AutoTestRow
{
    std::string mesh;
    int resolution;
    int nx, ny, nz;
    int totalVoxels;
    double buildMs;
    int validCount;
    int totalCount;
    double avgGap;
    double queryMs;
};

void runAutoTest(const std::string& dataDir)
{
    namespace fs = std::filesystem;

    std::cout << "\n============================================================\n";
    std::cout << "  Auto-Test: Multi-Mesh × Multi-Resolution\n";
    std::cout << "============================================================\n";
    std::cout << "  Data directory: " << dataDir << "\n";
    std::cout << "  Meshes:         " << kTestMeshes.size() << "\n";
    std::cout << "  Resolutions:    ";
    for (int r : kResolutions) std::cout << r << " ";
    std::cout << "\n\n";

    std::vector<AutoTestRow> rows;

    for (const auto& meshName : kTestMeshes)
    {
        std::string meshPath = (fs::path(dataDir) / meshName).string();
        if (!fs::exists(meshPath))
        {
            std::cout << "  [SKIP] " << meshName << " (not found)\n";
            continue;
        }

        TriangleMesh mesh = loadObjFile(meshPath);
        std::cout << "--- " << meshName << " (" << mesh.vertices.size()
                  << "v / " << mesh.triangles.size() << "t) ---\n";

        for (int res : kResolutions)
        {
            BuilderParams params;
            params.targetResolution = res;
            params.maxResolution = res;  // cap at target
            params.padding = 0.1;

            // Build
            auto t0 = std::chrono::high_resolution_clock::now();
            SdfData data = SdfBuilder::build(mesh, params);
            auto t1 = std::chrono::high_resolution_clock::now();
            double buildMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

            const auto& spec = data.spec;

            // Save + Load (I/O round-trip)
            std::string sdfFile = "autotest_" + meshName.substr(0, meshName.size() - 4)
                                   + "_r" + std::to_string(res) + ".sdf";
            saveSdf(data, sdfFile);
            auto loadedData = loadSdf(sdfFile);

            SdfOracle oracle(std::move(loadedData));

            // Query validation
            auto testPts = generateTestPoints(spec);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto results = oracle.queryBatch(testPts, QueryMode::SecondOrder);
            auto t3 = std::chrono::high_resolution_clock::now();
            double queryMs = std::chrono::duration<double, std::milli>(t3 - t2).count();

            int valid = 0;
            double sumGap = 0.0;
            for (const auto& r : results)
            {
                if (r.valid) { ++valid; sumGap += std::abs(r.gap); }
            }

            AutoTestRow row;
            row.mesh = meshName;
            row.resolution = res;
            row.nx = spec.nx();
            row.ny = spec.ny();
            row.nz = spec.nz();
            row.totalVoxels = spec.totalVoxels();
            row.buildMs = buildMs;
            row.validCount = valid;
            row.totalCount = static_cast<int>(testPts.size());
            row.avgGap = (valid > 0) ? sumGap / valid : 0.0;
            row.queryMs = queryMs;
            rows.push_back(row);

            std::cout << "  R=" << std::setw(4) << res
                      << "  grid=" << spec.nx() << "x" << spec.ny() << "x" << spec.nz()
                      << " (" << spec.totalVoxels() << " vox)"
                      << "  build=" << std::setprecision(1) << buildMs / 1000.0 << "s"
                      << "  valid=" << valid << "/" << testPts.size()
                      << "  query=" << std::setprecision(2) << queryMs << "ms"
                      << "\n";

            // Clean up temp file
            fs::remove(sdfFile);
        }
    }

    // Summary table
    std::cout << "\n============================================================\n";
    std::cout << "  Summary Table\n";
    std::cout << "============================================================\n";
    std::cout << std::left
              << std::setw(18) << "Mesh"
              << std::setw(6)  << "Res"
              << std::setw(16) << "Grid"
              << std::setw(12) << "Voxels"
              << std::setw(10) << "Build(s)"
              << std::setw(10) << "Valid"
              << std::setw(12) << "AvgGap"
              << std::setw(10) << "Query(ms)"
              << "\n";
    std::cout << std::string(94, '-') << "\n";

    for (const auto& r : rows)
    {
        std::string gridStr = std::to_string(r.nx) + "x" +
                              std::to_string(r.ny) + "x" +
                              std::to_string(r.nz);
        std::cout << std::left
                  << std::setw(18) << r.mesh
                  << std::setw(6)  << r.resolution
                  << std::setw(16) << gridStr
                  << std::setw(12) << r.totalVoxels
                  << std::setprecision(2) << std::setw(10) << r.buildMs / 1000.0
                  << std::setw(10) << (std::to_string(r.validCount) + "/" + std::to_string(r.totalCount))
                  << std::setprecision(6) << std::setw(12) << r.avgGap
                  << std::setprecision(2) << std::setw(10) << r.queryMs
                  << "\n";
    }

    // Verify all dimensions are standard resolutions
    std::cout << "\n--- Resolution Validation ---\n";
    bool allOk = true;
    auto isStandard = [](int n) {
        static constexpr int levels[] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};
        for (int v : levels) if (v == n) return true;
        return false;
    };
    for (const auto& r : rows)
    {
        bool okX = isStandard(r.nx);
        bool okY = isStandard(r.ny);
        bool okZ = isStandard(r.nz);
        if (!okX || !okY || !okZ)
        {
            std::cout << "  FAIL: " << r.mesh << " R=" << r.resolution
                      << " grid=" << r.nx << "x" << r.ny << "x" << r.nz << "\n";
            allOk = false;
        }
    }
    if (allOk)
        std::cout << "  All grid dimensions are standard power-of-2 resolutions. PASS\n";
}

// ---------- Main ----------

int main(int argc, char* argv[])
{
    std::cout << "========================================\n";
    std::cout << "  SdfOracle - Demo/Validator\n";
    std::cout << "========================================\n\n";

    DemoConfig cfg = parseArgs(argc, argv);

    // --- Auto-test mode ---
    if (cfg.autoTest)
    {
        std::string dataDir = cfg.dataDir;
        if (dataDir.empty())
        {
            // Default: relative to build/Release/
            dataDir = "../../../../data/input/obj_library";
        }
        runAutoTest(dataDir);
        std::cout << "\n=== Done ===\n";
        return 0;
    }

    SdfOracle oracle;
    std::string benchmarkSourcePath;
    std::string benchmarkSdfPath;
    double benchmarkBuildMs = 0.0;
    double benchmarkLoadMs = 0.0;
    TriangleMesh loadedMeshForCompare;

    if (!cfg.loadPath.empty())
    {
        // Load from file
        std::cout << "[Load] Loading from: " << cfg.loadPath << "\n";
        auto t0 = std::chrono::high_resolution_clock::now();
        oracle = SdfOracle::loadFromFile(cfg.loadPath);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        benchmarkLoadMs = ms;
        benchmarkSourcePath = cfg.loadPath;
        benchmarkSdfPath = cfg.loadPath;

        const auto& spec = oracle.gridSpec();
        std::cout << "[Load] Layout: " << (oracle.isSparse() ? "sparse_blocks" : "dense") << "\n";
        std::cout << "[Load] Grid: " << spec.nx() << " x " << spec.ny()
                  << " x " << spec.nz() << " = " << spec.totalVoxels()
                  << " voxels (" << std::setprecision(2) << ms << " ms)\n";
        if (oracle.isSparse())
        {
            const SparseSdfData& sparse = oracle.sparseData();
            std::cout << "[Load] Sparse blocks: " << sparse.blocks.size()
                      << " stored voxels=" << sparse.totalStoredVoxels()
                      << " blockDim=" << sparse.blockDim
                      << " band=" << sparse.bandDistance
                      << " contactStride=" << sparse.contactCandidateStride
                      << "\n";
        }
    }
    else if (!cfg.contactAwareModelPath.empty())
    {
        std::cout << "[Build] Loading contact-aware model: "
                  << cfg.contactAwareModelPath << "\n";
        CornerNormalTriangleMesh mesh =
            loadContactAwareSurfaceFile(cfg.contactAwareModelPath);
        benchmarkSourcePath = cfg.contactAwareModelPath;
        std::cout << "[Build] Contact-aware mesh: " << mesh.vertices.size()
                  << " vertices, " << mesh.triangles.size()
                  << " triangles\n";

        ContactAwareBuilderParams params;
        params.targetResolution = cfg.resolution;
        params.maxResolution = cfg.resolution;
        params.padding = cfg.padding;
        params.maxCandidates = cfg.maxCandidates;
        params.featureDistanceTolerance = cfg.featureTolerance;
        params.normalMergeAngleDegrees = cfg.normalMergeAngleDegrees;
        params.tileDim = cfg.tileDim;
        params.storeHessian = cfg.storeHessian;
        params.blockDim = cfg.blockDim;
        params.haloVoxels = cfg.haloVoxels;
        params.bandVoxels = cfg.bandVoxels;
        params.bandDistance = cfg.bandDistance;

        auto t0 = std::chrono::high_resolution_clock::now();
        if (cfg.sparse)
        {
            SparseSdfData data = SdfBuilder::buildSparseContactAware(mesh, params);
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            benchmarkBuildMs = ms;

            const auto& spec = data.spec;
            std::cout << "[Build] Layout: sparse contact-aware"
                      << " blockDim=" << data.blockDim
                      << " halo=" << data.haloVoxels
                      << " band=" << data.bandDistance
                      << " tileDim=" << cfg.tileDim
                      << "\n";
            std::cout << "[Build] Grid: " << spec.nx() << " x " << spec.ny()
                      << " x " << spec.nz() << " = " << spec.totalVoxels()
                      << " dense voxels\n";
            std::cout << "[Build] Sparse blocks: " << data.blocks.size()
                      << " stored voxels=" << data.totalStoredVoxels()
                      << " candidates stride=" << data.contactCandidateStride
                      << "\n";
            std::cout << "[Build] Completed in " << std::setprecision(1)
                      << ms / 1000.0 << " s\n";

            std::string sdfPath = cfg.outPrefix + ".sdf";
            benchmarkSdfPath = sdfPath;
            std::cout << "[Save] Writing: " << sdfPath << "\n";
            if (saveSparseSdf(data, sdfPath))
                std::cout << "[Save] Done.\n";
            else
                std::cerr << "[Save] FAILED!\n";

            oracle = SdfOracle(std::move(data));
        }
        else
        {
            SdfData data = SdfBuilder::buildContactAware(mesh, params);
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            benchmarkBuildMs = ms;

            const auto& spec = data.spec;
            std::cout << "[Build] Layout: dense contact-aware\n";
            std::cout << "[Build] Grid: " << spec.nx() << " x " << spec.ny()
                      << " x " << spec.nz() << " = " << spec.totalVoxels()
                      << " voxels\n";
            std::cout << "[Build] Candidates: stride="
                      << data.contactCandidateStride
                      << " max=" << cfg.maxCandidates << "\n";
            std::cout << "[Build] Completed in " << std::setprecision(1)
                      << ms / 1000.0 << " s\n";

            std::string sdfPath = cfg.outPrefix + ".sdf";
            benchmarkSdfPath = sdfPath;
            std::cout << "[Save] Writing: " << sdfPath << "\n";
            if (saveSdf(data, sdfPath))
                std::cout << "[Save] Done.\n";
            else
                std::cerr << "[Save] FAILED!\n";

            oracle = SdfOracle(std::move(data));
        }
    }
    else if (!cfg.objPath.empty())
    {
        // Build from OBJ
        std::cout << "[Build] Loading mesh: " << cfg.objPath << "\n";
        TriangleMesh mesh = loadObjFile(cfg.objPath);
        loadedMeshForCompare = mesh;
        benchmarkSourcePath = cfg.objPath;
        std::cout << "[Build] Mesh: " << mesh.vertices.size() << " vertices, "
                  << mesh.triangles.size() << " triangles\n";

        if (cfg.sparse)
        {
            SparseBuilderParams params;
            params.targetResolution = cfg.resolution;
            params.maxResolution = cfg.resolution;  // cap at target
            params.padding = cfg.padding;
            params.blockDim = cfg.blockDim;
            params.haloVoxels = cfg.haloVoxels;
            params.bandVoxels = cfg.bandVoxels;
            params.bandDistance = cfg.bandDistance;
            params.storeHessian = cfg.storeHessian;

            auto t0 = std::chrono::high_resolution_clock::now();
            SparseSdfData data = SdfBuilder::buildSparse(mesh, params);
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            benchmarkBuildMs = ms;

            const auto& spec = data.spec;
            std::cout << "[Build] Layout: sparse_blocks"
                      << " blockDim=" << data.blockDim
                      << " halo=" << data.haloVoxels
                      << " band=" << data.bandDistance << "\n";
            std::cout << "[Build] Grid: " << spec.nx() << " x " << spec.ny()
                      << " x " << spec.nz() << " = " << spec.totalVoxels()
                      << " dense voxels\n";
            std::cout << "[Build] Sparse blocks: " << data.blocks.size()
                      << " stored voxels=" << data.totalStoredVoxels()
                      << "\n";
            std::cout << "[Build] Completed in " << std::setprecision(1)
                      << ms / 1000.0 << " s\n";

            std::string sdfPath = cfg.outPrefix + ".sdf";
            benchmarkSdfPath = sdfPath;
            std::cout << "[Save] Writing: " << sdfPath << "\n";
            if (saveSparseSdf(data, sdfPath))
                std::cout << "[Save] Done.\n";
            else
                std::cerr << "[Save] FAILED!\n";

            oracle = SdfOracle(std::move(data));
            compareDenseSparseIfRequested(cfg, loadedMeshForCompare, oracle);
        }
        else
        {
            BuilderParams params;
            params.targetResolution = cfg.resolution;
            params.maxResolution = cfg.resolution;  // cap at target
            params.padding = cfg.padding;

            auto t0 = std::chrono::high_resolution_clock::now();
            SdfData data = SdfBuilder::build(mesh, params);
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            benchmarkBuildMs = ms;

            const auto& spec = data.spec;
            std::cout << "[Build] Layout: dense\n";
            std::cout << "[Build] Grid: " << spec.nx() << " x " << spec.ny()
                      << " x " << spec.nz() << " = " << spec.totalVoxels()
                      << " voxels\n";
            std::cout << "[Build] Completed in " << std::setprecision(1)
                      << ms / 1000.0 << " s\n";

            // Save
            std::string sdfPath = cfg.outPrefix + ".sdf";
            benchmarkSdfPath = sdfPath;
            std::cout << "[Save] Writing: " << sdfPath << "\n";
            if (saveSdf(data, sdfPath))
                std::cout << "[Save] Done.\n";
            else
                std::cerr << "[Save] FAILED!\n";

            oracle = SdfOracle(std::move(data));
        }
    }
    else
    {
        std::cerr << "Usage:\n"
                  << "  SdfOracleDemo --obj <mesh.obj> [--resolution {64,128,256,512,1024}] [--out prefix]\n"
                  << "      [--sparse] [--block-dim 8] [--band-voxels 8] [--halo-voxels 2]\n"
                  << "      [--benchmark] [--query-count 1000000] [--csv results.csv]\n"
                  << "      [--compare-dense-sparse]\n"
                  << "  SdfOracleDemo --contact-aware-model <model.ncas> [--out prefix]\n"
                  << "      [--resolution 128] [--padding 0.1] [--max-candidates 4]\n"
                  << "      [--feature-tolerance d] [--normal-merge-angle deg]\n"
                  << "      [--tile-dim 8] [--contact-aware-sparse] [--block-dim 8]\n"
                  << "      [--band-voxels 8] [--halo-voxels 2] [--no-hessian]\n"
                  << "  SdfOracleDemo --load <file.sdf> [--query x y z] [--benchmark]\n"
                  << "  SdfOracleDemo --auto-test [--data-dir <obj_library_path>]\n";
        return 1;
    }

    // Custom query points
    if (!cfg.queryPoints.empty())
    {
        std::cout << "\n=== Custom Query Points ===\n";
        for (const auto& p : cfg.queryPoints)
        {
            std::cout << "\n  Point: (" << p.x() << ", " << p.y() << ", " << p.z() << ")\n";
            printResult("1st-order", oracle.query(p, QueryMode::FirstOrder));
            printResult("2nd-order", oracle.query(p, QueryMode::SecondOrder));
            printResult("Trilinear", oracle.query(p, QueryMode::Trilinear));
            printResult("Tricubic",  oracle.query(p, QueryMode::Tricubic));
            printResult("ContactAware", oracle.query(p, QueryMode::ContactAware));
        }
    }

    // Validation
    if (cfg.runValidation)
    {
        runValidation(oracle);
    }

    if (cfg.benchmark)
    {
        std::vector<QueryBenchmarkRow> rows = runQueryBenchmark(
            oracle,
            cfg.queryCount,
            benchmarkSourcePath,
            benchmarkSdfPath,
            benchmarkBuildMs,
            benchmarkLoadMs);
        appendBenchmarkCsv(cfg.csvPath, rows);
        if (!cfg.csvPath.empty())
        {
            std::cout << "[Benchmark] CSV appended: " << cfg.csvPath << "\n";
        }
    }

    std::cout << "\n=== Done ===\n";
    return 0;
}
