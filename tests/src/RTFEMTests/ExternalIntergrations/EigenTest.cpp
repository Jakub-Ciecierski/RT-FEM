#include <RTFEMTests/ExternalIntergrations/EigenTest.h>

#include <Eigen/Core>

void EigenTest::SetUp() {
}

void EigenTest::TearDown() {
}

TEST_F(EigenTest, IntergartionTest) {
    Eigen::Matrix<double, 3, 3> matrix;

    matrix(0, 0) = 0;
    EXPECT_EQ(0, matrix(0, 0));
}