#include "RTFEMTests/FEM/Solver/GlobalStiffnessAssemblerTest.h"

#include <RTFEMTests/Builder/FEMModelSampleBuilder.h>
#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/FEM/Solver/GlobalStiffnessAssembler.h>
#include <RTFEM/DataStructure/Matrix.h>

void GlobalStiffnessAssemblerTest::SetUp() {
    FEMModelSampleBuilder builder;
    fem_model = builder.CreateRandomFEMModel();

    global_stiffness_assembler_ = rtfem::make_unique<rtfem::GlobalStiffnessAssembler>();
}

void GlobalStiffnessAssemblerTest::TearDown() {
}

TEST_F(GlobalStiffnessAssemblerTest, Test1){
    auto global_stiffness_matrix = global_stiffness_assembler_->Compute(fem_model);
}