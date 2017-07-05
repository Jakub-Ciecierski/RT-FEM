#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>

TEST(DataTypes, DataTypes_ProperDoubleValue) {
    rtfem::Float x = 5.0;

    EXPECT_EQ(5.0, x);
}