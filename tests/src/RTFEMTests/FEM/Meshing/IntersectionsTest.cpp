#include "RTFEMTests/FEM/Meshing/IntersectionsTest.h"

#include <RTFEM/FEM/Meshing/TriangleMesh.h>
#include <RTFEM/DataTypes.h>
#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/FEM/Meshing/Intersections.h>

void IntersectionsTest::SetUp() {
    triangle_mesh_ = rtfem::make_unique<rtfem::TriangleMesh<float>>();
    triangle_mesh_->triangles = {
        // Bottom
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{0, 0, 0},
                               Eigen::Vector3<float>{0, 0, 1},
                               Eigen::Vector3<float>{1, 0, 0}},
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{0, 0, 1},
                               Eigen::Vector3<float>{1, 0, 0},
                               Eigen::Vector3<float>{1, 0, 1}},
        // Top
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{0, 1, 0},
                               Eigen::Vector3<float>{0, 1, 1},
                               Eigen::Vector3<float>{1, 1, 0}},
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{0, 1, 1},
                               Eigen::Vector3<float>{1, 1, 0},
                               Eigen::Vector3<float>{1, 1, 1}},
        // Left
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{0, 1, 1},
                               Eigen::Vector3<float>{0, 1, 0},
                               Eigen::Vector3<float>{0, 0, 0}},
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{0, 1, 1},
                               Eigen::Vector3<float>{0, 0, 1},
                               Eigen::Vector3<float>{0, 0, 0}},
        // Right
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{1, 1, 0},
                               Eigen::Vector3<float>{1, 1, 1},
                               Eigen::Vector3<float>{1, 0, 0}},
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{1, 1, 1},
                               Eigen::Vector3<float>{1, 0, 1},
                               Eigen::Vector3<float>{1, 0, 0}},
        // Back
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{0, 1, 0},
                               Eigen::Vector3<float>{1, 1, 0},
                               Eigen::Vector3<float>{0, 0, 0}},
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{1, 1, 0},
                               Eigen::Vector3<float>{1, 0, 0},
                               Eigen::Vector3<float>{0, 0, 0}},
        // Front
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{0, 1, 1},
                               Eigen::Vector3<float>{0, 0, 1},
                               Eigen::Vector3<float>{1, 0, 1}},
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{0, 1, 1},
                               Eigen::Vector3<float>{1, 1, 1},
                               Eigen::Vector3<float>{1, 0, 1}},

    };
}

void IntersectionsTest::TearDown() {
}

TEST_F(IntersectionsTest, Contains_PointInsideTheCube_ReturnsTrue) {
    Eigen::Vector3<float> point1(0.5f, 0.5f, 0.5f);
    Eigen::Vector3<float> point2(0.001f, 0.0f, 0.0f);
    Eigen::Vector3<float> point3(0.999f, 0.999f, 0.999f);
    Eigen::Vector3<float> point4(0.1f, 0.0f, 0.0f);
    Eigen::Vector3<float> point5(0.1f, 0.0f, 0.1f);
    Eigen::Vector3<float> point6(1.0f, 0.0f, 0.9f);

    EXPECT_EQ(true, rtfem::Contains(point1, *triangle_mesh_));
    EXPECT_EQ(true, rtfem::Contains(point2, *triangle_mesh_));
    EXPECT_EQ(true, rtfem::Contains(point3, *triangle_mesh_));
    EXPECT_EQ(true, rtfem::Contains(point4, *triangle_mesh_));
    EXPECT_EQ(true, rtfem::Contains(point5, *triangle_mesh_));
    EXPECT_EQ(true, rtfem::Contains(point6, *triangle_mesh_));
}

TEST_F(IntersectionsTest, Contains_PointInsideTheCube_ReturnsFalse) {
    Eigen::Vector3<float> point1(-0.1f, 0.0f, 0.0f);
    Eigen::Vector3<float> point2(-0.1f, 0.0f, 1.0f);
    Eigen::Vector3<float> point3(0.0f, 0.0f, 1.1f);
    Eigen::Vector3<float> point4(0.1f, 0.0f, 1.1f);
    Eigen::Vector3<float> point5(0.0f, 1.1f, 1.0f);

    Eigen::Vector3<float> point6(1.0f, 1.0f, 1.0f);
    Eigen::Vector3<float> point7(0.0f, 0.0f, 0.0f);

    EXPECT_EQ(false, rtfem::Contains(point1, *triangle_mesh_));
    EXPECT_EQ(false, rtfem::Contains(point2, *triangle_mesh_));
    EXPECT_EQ(false, rtfem::Contains(point3, *triangle_mesh_));
    EXPECT_EQ(false, rtfem::Contains(point4, *triangle_mesh_));
    EXPECT_EQ(false, rtfem::Contains(point5, *triangle_mesh_));
    EXPECT_EQ(false, rtfem::Contains(point6, *triangle_mesh_));
    EXPECT_EQ(false, rtfem::Contains(point7, *triangle_mesh_));
}