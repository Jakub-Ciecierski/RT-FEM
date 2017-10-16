#include <RTFEMTests/FEM/Meshing/TetrahedralizationTest.h>

#include <RTFEMTests/Builder/TriangleMeshBuilder.h>

#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/FEM/Meshing/Tetrahedralization.h>
#include <RTFEM/FEM/Meshing/TriangleMesh.h>

#include <RTFEM/FEM/FEMGeometry.h>

void TetrahedralizationTest::SetUp() {
    triangle_mesh_cube_ = TriangleMeshBuilder().BuildCube();
}

void TetrahedralizationTest::TearDown() {
}

TEST_F(TetrahedralizationTest, Compute_NumberOfElements) {
    rtfem::Tetrahedralization<float> tetrahedralization;
    auto fem_geometry = tetrahedralization.Compute(*triangle_mesh_cube_);

    constexpr unsigned int expected_vertex_count = 8;
    constexpr unsigned int expected_tetra_count = 6;

    EXPECT_EQ(expected_vertex_count, fem_geometry->vertices.size());
    EXPECT_EQ(expected_tetra_count, fem_geometry->finite_elements.size());
}

TEST_F(TetrahedralizationTest, Compute_SetVolumeConstriant) {
    rtfem::Tetrahedralization<float> tetrahedralization;

    rtfem::TetrahedralizationOptions options;
    options.maximum_volume = 0.1;
    tetrahedralization.SetOptions(options);

    auto fem_geometry = tetrahedralization.Compute(*triangle_mesh_cube_);

    constexpr unsigned int expected_vertex_count = 50;
    EXPECT_EQ(expected_vertex_count, fem_geometry->vertices.size());
}