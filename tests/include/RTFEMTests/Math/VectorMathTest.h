#include "gtest/gtest.h"

#include <memory>

namespace rtfem {
class VectorMath;
}

class VectorMathTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;

    std::unique_ptr<rtfem::VectorMath> vector_math_;
};