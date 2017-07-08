#include "gtest/gtest.h"

#include <memory>

namespace rtfem{
class TetrahedronFiniteElement;
class TetrahedronSolver;
}

class TetrahedronSolverTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    std::shared_ptr<rtfem::TetrahedronFiniteElement> finite_element_;
    std::unique_ptr<rtfem::TetrahedronSolver> solver_;
};