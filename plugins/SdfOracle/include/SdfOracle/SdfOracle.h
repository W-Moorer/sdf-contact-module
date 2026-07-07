#pragma once

#include <string>
#include <vector>

#include "SdfOracle/SdfCommon.h"
#include "SdfOracle/SdfData.h"

namespace sdf {

enum class QueryMode
{
    FirstOrder,
    SecondOrder,
    Trilinear,
    Tricubic,
    ContactAware
};

class SdfOracle
{
public:
    SdfOracle() = default;

    explicit SdfOracle(SdfData data);
    explicit SdfOracle(SparseSdfData data);

    static SdfOracle loadFromFile(const std::string& path);

    SdfQueryResult query(
        const Eigen::Vector3d& point,
        QueryMode mode = QueryMode::SecondOrder) const;

    SdfFastQueryResult queryFast(
        const Eigen::Vector3d& point,
        QueryMode mode = QueryMode::SecondOrder) const;

    std::vector<SdfQueryResult> queryBatch(
        const std::vector<Eigen::Vector3d>& points,
        QueryMode mode = QueryMode::SecondOrder) const;

    const SdfData& data() const { return m_data; }
    const SparseSdfData& sparseData() const { return m_sparseData; }
    bool isSparse() const { return m_sparse; }
    const GridSpec& gridSpec() const { return m_sparse ? m_sparseData.spec : m_data.spec; }

    bool isValid() const { return m_sparse ? !m_sparseData.empty() : !m_data.phi0.empty(); }

private:
    SdfQueryResult secondOrderQuery(const Eigen::Vector3d& p) const;
    SdfQueryResult firstOrderQuery(const Eigen::Vector3d& p) const;
    SdfQueryResult trilinearQuery(const Eigen::Vector3d& p) const;
    SdfQueryResult tricubicQuery(const Eigen::Vector3d& p) const;
    SdfQueryResult contactAwareQuery(const Eigen::Vector3d& p) const;

    SdfQueryResult sparseFirstOrderQuery(const Eigen::Vector3d& p) const;
    SdfQueryResult sparseSecondOrderQuery(const Eigen::Vector3d& p) const;
    SdfQueryResult sparseTrilinearQuery(const Eigen::Vector3d& p) const;
    SdfQueryResult sparseTricubicQuery(const Eigen::Vector3d& p) const;
    SdfQueryResult sparseContactAwareQuery(const Eigen::Vector3d& p) const;

    Eigen::Vector3d estimateWitness(int i, int j, int k) const;

    SdfData m_data;
    SparseSdfData m_sparseData;
    bool m_sparse{false};
};

} // namespace sdf
