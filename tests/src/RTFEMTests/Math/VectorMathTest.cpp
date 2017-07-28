#include <RTFEMTests/Math/VectorMathTest.h>

#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/DataStructure/Vector3.h>
#include <RTFEM/Math/VectorMath.h>

void VectorMathTest::SetUp() {
    vector_math_ = rtfem::make_unique<rtfem::VectorMath>();
}

void VectorMathTest::TearDown() {
}

TEST_F(VectorMathTest, Cross_OrthonormalBasis){
    rtfem::Vector3 v1(1, 0, 0);
    rtfem::Vector3 v2(0, 1, 0);

    rtfem::Vector3 expected_vector(0, 0, 1);

    EXPECT_EQ(expected_vector, vector_math_->Cross(v1, v2));
}

TEST_F(VectorMathTest, Magnitude_UnitLength){
    rtfem::Vector3 v(1, 0, 0);
    double expected_magnitude = 1.0;
    EXPECT_EQ(expected_magnitude, vector_math_->Magnitude(v));
}