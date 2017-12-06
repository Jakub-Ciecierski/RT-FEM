#include "RTFEMTests/FEM/Solver/FEMGlobalAssemblers/FEMGlobalDynamicAssemblerTest.h"

#include <RTFEMTests/Builder/FEMModelSampleBuilder.h>

#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMGlobalDynamicAssembler.h>
#include <RTFEM/Memory/UniquePointer.h>

void FEMGlobalDynamicAssemblerTest::SetUp() {
    FEMModelSampleBuilder builder;
    fem_model_ = builder.CreateRandomFEMModel();

    fem_assembler_
        = rtfem::make_unique<rtfem::FEMGlobalDynamicAssembler<double>>();
}

void FEMGlobalDynamicAssemblerTest::TearDown() {
}

TEST_F(FEMGlobalDynamicAssemblerTest, MassMatrix_IsSymetric){
    auto output = fem_assembler_->Compute(*fem_model_);
    EXPECT_EQ(output.global_mass.transpose(), output.global_mass);
}

TEST_F(FEMGlobalDynamicAssemblerTest, DampingMatrix_IsSymetric) {
    auto material = fem_model_->material();
    material.damping_mass = 1;
    material.damping_stiffness = 1;
    fem_model_->material(material);
    auto output = fem_assembler_->Compute(*fem_model_);
    EXPECT_EQ(output.global_damping.transpose(), output.global_damping);
}