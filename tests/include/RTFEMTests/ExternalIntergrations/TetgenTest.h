#include "gtest/gtest.h"

class TetgenTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;
};