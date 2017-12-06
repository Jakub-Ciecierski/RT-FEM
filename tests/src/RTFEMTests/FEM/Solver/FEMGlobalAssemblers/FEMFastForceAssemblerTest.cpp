#include "RTFEMTests/FEM/Solver/FEMGlobalAssemblers/FEMFastForceAssemblerTest.h"

#include <RTFEMTests/Builder/FEMModelSampleBuilder.h>

#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMGlobalDynamicAssembler.h>
#include <RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMFastForceAssembler.h>

void FEMFastForceAssemblerTest::SetUp() {
    FEMModelSampleBuilder builder;
    fem_model_ = builder.CreateRandomFEMModel();
}

void FEMFastForceAssemblerTest::TearDown() {
}

TEST_F(FEMFastForceAssemblerTest, BodyForce_SameResultAsNormalAssembler){
    auto body_force = Eigen::Vector3<double>(0, -9.82, 0);
    fem_model_->AddDynamicBodyForce(body_force);

    rtfem::FEMGlobalDynamicAssembler<double> assembler;
    auto assembler_data = assembler.Compute(*fem_model_);

    Eigen::Vector<double, Eigen::Dynamic> fast_global_force =
            Eigen::Vector<double, Eigen::Dynamic>::Zero(
                    assembler_data.global_force.size());

    rtfem::FEMFastForceAssembler<double> fast_force_assembler;
    fast_force_assembler.Assemble(fem_model_->fem_geometry(),
                                  fem_model_->total_body_force(),
                                  fem_model_->material(),
                                  fast_global_force);

    EXPECT_EQ(assembler_data.global_force, fast_global_force);
}


TEST_F(FEMFastForceAssemblerTest, TractionForce_SameResultAsNormalAssembler){
    fem_model_->fem_geometry().triangle_faces[0].traction_force = 2;

    rtfem::FEMGlobalDynamicAssembler<double> assembler;
    auto assembler_data = assembler.Compute(*fem_model_);

    Eigen::Vector<double, Eigen::Dynamic> fast_global_force =
            Eigen::Vector<double, Eigen::Dynamic>::Zero(
                    assembler_data.global_force.size());

    rtfem::FEMFastForceAssembler<double> fast_force_assembler;
    fast_force_assembler.Assemble(fem_model_->fem_geometry(),
                                  fem_model_->total_body_force(),
                                  fem_model_->material(),
                                  fast_global_force);

    // TODO
    //EXPECT_EQ(assembler_data.global_force, fast_global_force);
}