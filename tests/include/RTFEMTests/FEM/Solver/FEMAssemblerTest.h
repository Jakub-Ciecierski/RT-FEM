#include "gtest/gtest.h"

#include <memory>

namespace rtfem{
class FEMModel;
class FEMAssembler;
}

class FEMAssemblerTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;

    std::shared_ptr<rtfem::FEMModel> fem_model;
    std::unique_ptr<rtfem::FEMAssembler> fem_assembler_;

};