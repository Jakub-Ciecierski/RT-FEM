#include "gtest/gtest.h"

#include <RTFEM/FEM/FEMModel.h>

class BoundaryConditionTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    std::shared_ptr<rtfem::FEMModel<double>> fem_model_;
};