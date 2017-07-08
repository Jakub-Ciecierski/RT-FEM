#include "RTFEMTests/FEM/VertexTest.h"

#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/DataStructure/Vector3.h>
#include <RTFEM/Memory/UniquePointer.h>

void VertexTest::SetUp() {
    default_vertex_.id = 2;
    default_vertex_.coordinates = rtfem::Vector3();
    default_vertex_.vertex = rtfem::make_unique<rtfem::Vertex>(default_vertex_.id,
                                                               std::move(default_vertex_.coordinates));
}

void VertexTest::TearDown() {
}

TEST_F(VertexTest, CreatedVertex_ProperIDValue){
    EXPECT_EQ(default_vertex_.id, default_vertex_.vertex->id());
}