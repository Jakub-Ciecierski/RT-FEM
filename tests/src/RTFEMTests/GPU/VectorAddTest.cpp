#include "RTFEMTests/GPU/VectorAddTest.h"

#include <RTFEM/GPU/VectorAdd/VectorAdd.cuh>

void VectorAddTest::SetUp() {
}

void VectorAddTest::TearDown() {
}

TEST_F(VectorAddTest, CUDA_IntegrationTest){
    rtfem::VectorAddition vector_addition;
    vector_addition.RunVectorAddTest();
}