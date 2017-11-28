#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>

#include <memory>

namespace rtfem {

template<class T>
class FEMModel;

template<class T>
class FEMGlobalDynamicAssembler;
}

class FEMGlobalDynamicAssemblerTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    std::shared_ptr<rtfem::FEMModel<double>> fem_model_;
    std::unique_ptr<rtfem::FEMGlobalDynamicAssembler<double>> fem_assembler_;
};