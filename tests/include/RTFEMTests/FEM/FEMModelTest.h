#include "gtest/gtest.h"

#include <memory>
#include <RTFEM/DataTypes.h>

namespace rtfem{
template<class T>
class FEMModel;
}

struct FEMModelPack{
    unsigned int vertex_count;
    unsigned int finite_element_count;
    std::shared_ptr<rtfem::FEMModel<double>> fem_model;
};

class FEMModelTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    FEMModelPack fem_model_pack_;
};