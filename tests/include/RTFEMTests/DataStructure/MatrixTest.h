#ifndef PROJECT_MATRIXTEST_H
#define PROJECT_MATRIXTEST_H

#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>

#include <memory>

namespace rtfem{
class Matrix;
}

struct MatrixTestPack{
    unsigned int row_count;
    unsigned int column_count;
    std::unique_ptr<rtfem::Matrix> matrix;
};

class MatrixTest : public ::testing::Test{
protected:
    virtual void SetUp() override;
    virtual void TearDown() override;

    MatrixTestPack& small_matrix(){return small_matrix_;};

private:
    std::unique_ptr<rtfem::Matrix> CreateMatrix(unsigned int r, unsigned int c);

    MatrixTestPack small_matrix_;
};



#endif //PROJECT_MATRIXTEST_H
