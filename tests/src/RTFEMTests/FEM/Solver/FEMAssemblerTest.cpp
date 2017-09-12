#include "RTFEMTests/FEM/Solver/FEMAssemblerTest.h"

#include <RTFEMTests/Builder/FEMModelSampleBuilder.h>

#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/FEM/Solver/FEMAssembler.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/FEM/Material.h>
#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/FEMGeometry.h>

#include <cmath>

void FEMAssemblerTest::SetUp() {
    FEMModelSampleBuilder builder;
    fem_model_ = builder.CreateRandomFEMModel();

    fem_assembler_ = rtfem::make_unique<rtfem::FEMAssembler<double>>();
}

void FEMAssemblerTest::TearDown() {}

TEST_F(FEMAssemblerTest,
       SingleFiniteElement_CorrectGlobalStiffnessMatrix_Mathematica){
    auto fem_geometry = rtfem::make_unique<rtfem::FEMGeometry<double>>();
    fem_geometry->finite_elements =
        std::vector<std::shared_ptr<rtfem::FiniteElement<double>>>(1);
    fem_geometry->vertices = std::vector<std::shared_ptr<rtfem::Vertex<double>>>(4);

    fem_geometry->vertices[0] = std::make_shared<rtfem::Vertex<double>>(0, Eigen::Vector3<double>(2,3,4));
    fem_geometry->vertices[1] = std::make_shared<rtfem::Vertex<double>>(1, Eigen::Vector3<double>(6,3,2));
    fem_geometry->vertices[2] = std::make_shared<rtfem::Vertex<double>>(2, Eigen::Vector3<double>(2,5,1));
    fem_geometry->vertices[3] = std::make_shared<rtfem::Vertex<double>>(3, Eigen::Vector3<double>(4,3,6));

    fem_geometry->finite_elements[0] = std::make_shared<rtfem::TetrahedronFiniteElement<double>>(
            fem_geometry->vertices[0],
            fem_geometry->vertices[1],
            fem_geometry->vertices[2],
            fem_geometry->vertices[3]);

    auto fem_model = std::make_shared<rtfem::FEMModel<double>>(
        std::move(fem_geometry),
            rtfem::Material<double>{480, 1.0 / 3.0});

    auto fem_assembler_data = fem_assembler_->Compute(fem_model);

    for(unsigned int i = 0; i < 12; i++){
        for(unsigned int j = 0; j < 12; j++){
            fem_assembler_data.global_stiffness(i, j) = (int)(fem_assembler_data.global_stiffness(i,j));
        }
    }
    Eigen::Matrix<double, 12, 12> expected_stiffness;
    expected_stiffness <<
                 745,   540,    120,    -5,     30,     60,     -270,   -240,   0,      -470,   -330,   -180 ,
                 540,   1720,   270,    -120,   520,    210,    -120,   -1080,  -60,    -300,   -1160,  -420 ,
                 120,   270,    565,    0,      150,    175,    0,      -120,    -270,  -120,   -300,   -470 ,
                 -5,    -120,   0,      145,    -90,    -60,    -90,    120,    0,      -50,    90,     60 ,
                 30,    520,    150,    -90,    220,    90,     60,     -360,   -60,    0,      -380,   -180 ,
                 60,    210,    175,    -60,    90,     145,    0,      -120,   -90,    0,      -180,   -230 ,
                 -270,  -120,   0,      -90,    60,     0,      180,    0,      0,      180,    60,     0 ,
                 -240,  -1080,  -120,   120,    -360,   -120,   0,      720,    0,      120,    720,    240 ,
                 0,     -60,    -270,   0,      -60,    -90,    0,      0,      180,    0,      120,    180 ,
                 -470,  -300,   -120,   -50,    0,      0,      180,    120,    0,      340,    180,    120 ,
                 -330,  -1160,  -300,   90,     -380,   -180,   60,     720,    120,    180,    820,    360 ,
                 -180,  -420,   -470,   60,     -180,   -230,   0,      240,    180,    120,    360,    520;

    EXPECT_EQ(expected_stiffness, fem_assembler_data.global_stiffness);
};

TEST_F(FEMAssemblerTest, FEMAssembler_Compute_ProperForceVectorDimensions){
    auto fem_model = FEMModelSampleBuilder().CreateRandomFEMModel();
    auto fem_assembler_data = fem_assembler_->Compute(fem_model);

    unsigned int expected_row_count = fem_model_->fem_geometry().vertices.size() * 3;
    EXPECT_EQ(expected_row_count, fem_assembler_data.global_force.rows());
    EXPECT_EQ((unsigned int)1, fem_assembler_data.global_force.cols());
}

TEST_F(FEMAssemblerTest, FEMAssembler_Compute_ProperStiffnessDimensions){
    auto fem_model = FEMModelSampleBuilder().CreateRandomFEMModel();
    auto fem_assembler_data = fem_assembler_->Compute(fem_model);

    unsigned int expected_row_count = fem_model_->fem_geometry().vertices.size() * 3;
    unsigned int expected_column_count = expected_row_count;
    EXPECT_EQ(expected_row_count, fem_assembler_data.global_stiffness.rows());
    EXPECT_EQ(expected_column_count, fem_assembler_data.global_stiffness.cols());
}