#include "gtest/gtest.h"

#include <memory>
#include <RTFEM/DataTypes.h>

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