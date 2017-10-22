#include "RTFEMTests/FEM/Meshing/TriangleFaceTest.h"

#include <RTFEM/FEM/Meshing/TriangleMesh.h>

void TriangleFaceTest::SetUp() {
}

void TriangleFaceTest::TearDown() {
}

TEST_F(TriangleFaceTest, TwoTriangleFaces_Equal1){
    rtfem::TriangleFace face1{1, 2, 3};
    rtfem::TriangleFace face2{2, 1, 3};

    const bool expected_value = true;

    EXPECT_EQ(expected_value, face1 == face2);
}

TEST_F(TriangleFaceTest, TwoTriangleFaces_Equal2){
    rtfem::TriangleFace face1{1, 2, 3};
    rtfem::TriangleFace face2{3, 1, 2};

    const bool expected_value = true;

    EXPECT_EQ(expected_value, face1 == face2);
}

TEST_F(TriangleFaceTest, TwoTriangleFaces_Equal3){
    rtfem::TriangleFace face1{1, 2, 3};
    rtfem::TriangleFace face2{3, 2, 1};

    const bool expected_value = true;

    EXPECT_EQ(expected_value, face1 == face2);
}

TEST_F(TriangleFaceTest, TwoTriangleFaces_Equal4){
    rtfem::TriangleFace face1{2, 1, 3};
    rtfem::TriangleFace face2{2, 1, 3};

    const bool expected_value = true;

    EXPECT_EQ(expected_value, face1 == face2);
}

TEST_F(TriangleFaceTest, TwoTriangleFaces_Equal5){
    rtfem::TriangleFace face1{1, 2, 3};
    rtfem::TriangleFace face2{2, 3, 1};

    const bool expected_value = true;

    EXPECT_EQ(expected_value, face1 == face2);
}

TEST_F(TriangleFaceTest, TwoTriangleFaces_NotEqual1){
    rtfem::TriangleFace face1{1, 2, 3};
    rtfem::TriangleFace face2{2, 2, 3};

    const bool expected_value = false;

    EXPECT_EQ(expected_value, face1 == face2);
}