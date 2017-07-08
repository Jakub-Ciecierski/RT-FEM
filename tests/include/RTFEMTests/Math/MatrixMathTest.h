#include "gtest/gtest.h"

namespace rtfem{
class Matrix;
class MatrixMath;
}

class MatrixMathTest : public ::testing::Test {
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    std::unique_ptr<rtfem::Matrix> matrix_;
    std::unique_ptr<rtfem::MatrixMath> matrix_math_;
};