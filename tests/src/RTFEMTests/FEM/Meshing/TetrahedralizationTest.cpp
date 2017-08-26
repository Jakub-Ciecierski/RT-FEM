#include <RTFEMTests/FEM/Meshing/TetrahedralizationTest.h>

#include <RTFEMTests/Builder/TriangleMeshBuilder.h>

#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/FEM/Meshing/Tetrahedralization.h>
#include <RTFEM/FEM/Meshing/TriangleMesh.h>

#include <RTFEM/FEM/FEMGeometry.h>

void TetrahedralizationTest::SetUp() {
    tetrahedralization_ = rtfem::make_unique<rtfem::Tetrahedralization<float>>();
    triangle_mesh_cube_ = TriangleMeshBuilder().BuildCube();
}

void TetrahedralizationTest::TearDown() {
}

TEST_F(TetrahedralizationTest, Compute_NumberOfElements){
    auto fem_geometry = tetrahedralization_->Compute(*triangle_mesh_cube_);

    constexpr unsigned int expected_vertex_count = 8;
    constexpr unsigned int expected_tetra_count = 6;
    constexpr unsigned int expected_face_count = 12;

    EXPECT_EQ(expected_vertex_count, fem_geometry.vertices.size());
    EXPECT_EQ(expected_tetra_count, fem_geometry.finite_elements.size());
    EXPECT_EQ(expected_face_count, fem_geometry.finite_element_faces.size());
}