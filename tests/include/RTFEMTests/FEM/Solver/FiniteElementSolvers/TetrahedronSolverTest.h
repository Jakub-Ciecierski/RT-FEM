#include "gtest/gtest.h"

#include <memory>
#include <RTFEM/DataTypes.h>

namespace rtfem{
template<class T>
class TetrahedronFiniteElement;

template<class T>
class TetrahedronSolver;
}

class TetrahedronSolverTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    std::shared_ptr<rtfem::TetrahedronFiniteElement<double>> finite_element_;
    std::unique_ptr<rtfem::TetrahedronSolver<double>> solver_;
};