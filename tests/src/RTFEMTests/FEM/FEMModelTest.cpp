#include "RTFEMTests/FEM/FEMModelTest.h"

#include <RTFEMTests/Builder/FEMModelSampleBuilder.h>

#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/FEMGeometry.h>

void FEMModelTest::SetUp() {
    FEMModelSampleBuilder builder;

    fem_model_pack_.finite_element_count = builder.finite_element_count();
    fem_model_pack_.vertex_count = builder.vertex_count();

    fem_model_pack_.fem_model = builder.CreateRandomFEMModel();
}

void FEMModelTest::TearDown() {
}

TEST_F(FEMModelTest, CreatedFEMModel_ProperVertexCount) {
    EXPECT_EQ(fem_model_pack_.vertex_count,
              fem_model_pack_.fem_model->fem_geometry().vertices.size());
}

TEST_F(FEMModelTest, CreatedFEMModel_ProperFiniteElementCount) {
    EXPECT_EQ(fem_model_pack_.finite_element_count,
              fem_model_pack_.fem_model->fem_geometry().finite_elements.size());
}