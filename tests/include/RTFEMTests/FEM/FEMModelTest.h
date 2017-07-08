#include "gtest/gtest.h"

#include <memory>
#include <RTFEM/DataTypes.h>

namespace rtfem{
class FEMModel;
}

struct FEMModelPack{
    rtfem::UInt vertex_count;
    rtfem::UInt finite_element_count;
    std::shared_ptr<rtfem::FEMModel> fem_model;
};

class FEMModelTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    FEMModelPack fem_model_pack_;
};