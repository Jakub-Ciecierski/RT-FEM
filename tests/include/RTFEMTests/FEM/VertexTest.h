#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/Vector3.h>

namespace rtfem{
class Vertex;
}

struct VertexPack{
    rtfem::UInt id;
    rtfem::Vector3 coordinates;
    std::unique_ptr<rtfem::Vertex> vertex;
};

class VertexTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;

    VertexPack default_vertex_;
};