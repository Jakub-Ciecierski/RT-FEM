#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>

#include <memory>

namespace rtfem{
template<class T>
class TetrahedronFiniteElement;
}

class TetrahedronFiniteElementTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    std::unique_ptr<rtfem::TetrahedronFiniteElement<double>> tetrahedron_;
};