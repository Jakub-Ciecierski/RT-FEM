#include "RTFEMTests/DataStructure/Vector3Test.h"

#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/DataStructure/Vector3.h>

void Vector3Test::SetUp(){
    vector_ = rtfem::make_unique<rtfem::Vector3>();
}

void Vector3Test::TearDown(){

}

TEST_F(Vector3Test, CreatedVector_InitiedWithZeros) {
    EXPECT_EQ(0, vector_->x);
    EXPECT_EQ(0, vector_->y);
    EXPECT_EQ(0, vector_->z);
}

TEST_F(Vector3Test, CreatedVector_TrippleConstructor_ProperValues) {
    rtfem::Float expexted_x = 1.0;
    rtfem::Float expexted_y = 2.0;
    rtfem::Float expexted_z = 3.0;

    rtfem::Vector3 vec3(expexted_x, expexted_y, expexted_z);

    EXPECT_EQ(expexted_x, vec3.x);
    EXPECT_EQ(expexted_y, vec3.y);
    EXPECT_EQ(expexted_z, vec3.z);
}

TEST_F(Vector3Test, Vector_ChangedXValue_ValueRegisteredProperly) {
    rtfem::Float expected_value = 5.5f;
    vector_->x = expected_value;

    EXPECT_EQ(expected_value, vector_->x);
}

TEST_F(Vector3Test, Vector_ChangedYValue_ValueRegisteredProperly) {
    rtfem::Float expected_value = 5.5f;
    vector_->y = expected_value;

    EXPECT_EQ(expected_value, vector_->y);
}

TEST_F(Vector3Test, Vector_ChangedZValue_ValueRegisteredProperly) {
    rtfem::Float expected_value = 5.5f;
    vector_->z = expected_value;

    EXPECT_EQ(expected_value, vector_->z);
}