#include "RTFEMTests/FEM/FiniteElements/TetrahedronFiniteElementTest.h"

#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/Memory/UniquePointer.h>

void TetrahedronFiniteElementTest::SetUp() {
    tetrahedron_ = rtfem::make_unique<rtfem::TetrahedronFiniteElement<double>>(
        0,0,0,0
    );
}

void TetrahedronFiniteElementTest::TearDown() {
}

TEST_F(TetrahedronFiniteElementTest, CreatedTetrahedron_ProperType){
    EXPECT_EQ(rtfem::FiniteElementType::Tetrahedron, tetrahedron_->type());
}

TEST_F(TetrahedronFiniteElementTest, CreatedTetrahedron_ProperVertexCount){
    constexpr unsigned int expected_count = 4;
    EXPECT_EQ(expected_count, tetrahedron_->GetVertexCount());
}