#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>

namespace rtfem{

template<class T>
class SparseMatrixCSR;

}

class GPUSparseLinearSolverTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    Eigen::Matrix<double, 5, 5> dense_matrix;
    Eigen::Vector<double, 5> b;
    Eigen::Vector<double, 5> x;

    std::unique_ptr<rtfem::SparseMatrixCSR<double>> sparse_matrix;
};