#include "gtest/gtest.h"

#include <memory>
#include <RTFEM/DataTypes.h>

namespace rtfem{
template<class T>
class FEMModel;

template<class T>
class FEMAssembler;
}

class FEMAssemblerTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;

    std::shared_ptr<rtfem::FEMModel<double>> fem_model_;
    std::unique_ptr<rtfem::FEMAssembler<double>> fem_assembler_;

};