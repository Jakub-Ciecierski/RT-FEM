#include "gtest/gtest.h"

#include <memory>

namespace rtfem{
class TetrahedronFiniteElement;
}

class TetrahedronFiniteElementTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    std::unique_ptr<rtfem::TetrahedronFiniteElement> tetrahedron_;
};