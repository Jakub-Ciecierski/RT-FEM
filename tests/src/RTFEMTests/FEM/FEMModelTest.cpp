#include "RTFEMTests/FEM/FEMModelTest.h"

#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>

#include <RTFEMTests/Builder/FEMModelSampleBuilder.h>

void FEMModelTest::SetUp() {
    FEMModelSampleBuilder builder;

    fem_model_pack_.finite_element_count = builder.finite_element_count();
    fem_model_pack_.vertex_count = builder.vertex_count();

    fem_model_pack_.fem_model = builder.CreateRandomFEMModel();
}

void FEMModelTest::TearDown() {
}

TEST_F(FEMModelTest, CreatedFEMModel_ProperVertexCount){
    EXPECT_EQ(fem_model_pack_.vertex_count, fem_model_pack_.fem_model->VertexCount());
}

TEST_F(FEMModelTest, CreatedFEMModel_ProperFiniteElementCount){
    EXPECT_EQ(fem_model_pack_.finite_element_count, fem_model_pack_.fem_model->FiniteElementCount());
}