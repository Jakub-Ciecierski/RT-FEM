#include "gtest/gtest.h"

#include <memory>

namespace rtfem {
template<class T>
struct TriangleMesh;
}

class IntersectionsTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;

    std::unique_ptr<rtfem::TriangleMesh<float>> triangle_mesh_;
};