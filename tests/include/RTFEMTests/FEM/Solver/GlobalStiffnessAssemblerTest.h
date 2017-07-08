#include "gtest/gtest.h"

#include <memory>

namespace rtfem{
class FEMModel;
class GlobalStiffnessAssembler;
}

class GlobalStiffnessAssemblerTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;

    std::shared_ptr<rtfem::FEMModel> fem_model;
    std::unique_ptr<rtfem::GlobalStiffnessAssembler> global_stiffness_assembler_;

};