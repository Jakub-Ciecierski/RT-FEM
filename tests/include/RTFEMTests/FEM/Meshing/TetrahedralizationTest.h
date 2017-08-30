#include "gtest/gtest.h"

#include <memory>

namespace rtfem {
template<class T>
class Tetrahedralization;

template<class T>
struct TriangleMeshIndexed;

}

class TetrahedralizationTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;

    std::unique_ptr<rtfem::TriangleMeshIndexed<float>> triangle_mesh_cube_;
};