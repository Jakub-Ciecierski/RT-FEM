#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>

#include <memory>

namespace rtfem{
template<class T>
class TetrahedronFiniteElement;

template<class T>
class TetrahedronSolver;

template<class T>
class Vertex;
}

class TetrahedronSolverTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    std::shared_ptr<rtfem::TetrahedronFiniteElement<double>> finite_element_;
    std::vector<std::shared_ptr<rtfem::Vertex<double>>> vertices_;

    std::unique_ptr<rtfem::TetrahedronSolver<double>> solver_;
};