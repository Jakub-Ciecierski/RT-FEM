#include "gtest/gtest.h"

class FEMSolverTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;
};