#include "RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver.h"

#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/Math/MatrixMath.h>

namespace rtfem {

TetrahedronSolver::TetrahedronSolver() {}

TetrahedronSolver::~TetrahedronSolver() {}

FiniteElementSolverData TetrahedronSolver::Solve(std::shared_ptr<FiniteElement> finite_element){
    auto a0 = ComputeA(0, finite_element);
    auto a1 = ComputeA(1, finite_element);
    auto a2 = ComputeA(2, finite_element);
    auto a3 = ComputeA(3, finite_element);

    auto b0 = ComputeB(0, finite_element);
    auto b1 = ComputeB(1, finite_element);
    auto b2 = ComputeB(2, finite_element);
    auto b3 = ComputeB(3, finite_element);

    auto c0 = ComputeC(0, finite_element);
    auto c1 = ComputeC(1, finite_element);
    auto c2 = ComputeC(2, finite_element);
    auto c3 = ComputeC(3, finite_element);

    auto d0 = ComputeD(0, finite_element);
    auto d1 = ComputeD(1, finite_element);
    auto d2 = ComputeD(2, finite_element);
    auto d3 = ComputeD(3, finite_element);

    auto volume = ComputeVolume(finite_element);
    if(volume == 0){
        throw std::invalid_argument("TetrahedronSolver::Solve: Element with 0 volume");
    }

    FiniteElementSolverData data;
    data.volume = volume;

    data.geometry_matrix = Matrix(6, 12);
    AssemblyGeometryMatrix(data.geometry_matrix, 0, b0, c0, d0);
    AssemblyGeometryMatrix(data.geometry_matrix, 3, b1, c1, d1);
    AssemblyGeometryMatrix(data.geometry_matrix, 6, b2, c2, d2);
    AssemblyGeometryMatrix(data.geometry_matrix, 9, b3, c3, d3);

    data.geometry_matrix = data.geometry_matrix / (6*volume);

    return data;
}

Float TetrahedronSolver::ComputeA(UInt vertex_local_index,
                                  const std::shared_ptr<FiniteElement> finite_element){
    const UInt n = 3;
    Matrix m(n, n);
    UInt current_column = 0;
    for(UInt i = 0; i < finite_element->GetVertexCount(); i++){
        if(i == vertex_local_index)
            continue;
        const auto& vertex = finite_element->vertices()[i];
        m[0][current_column] = vertex->coordinates().x;
        m[1][current_column] = vertex->coordinates().y;
        m[2][current_column] = vertex->coordinates().z;
        current_column++;
    }
    MatrixMath matrix_math;

    return matrix_math.ComputeDeterminant(m);
}

Float TetrahedronSolver::ComputeB(UInt vertex_local_index,
                                  const std::shared_ptr<FiniteElement> finite_element){
    const UInt n = 3;
    Matrix m(n, n);
    UInt current_column = 0;
    for(UInt i = 0; i < finite_element->GetVertexCount(); i++){
        if(i == vertex_local_index)
            continue;

        const auto& vertex = finite_element->vertices()[i];
        m[0][current_column] = 1;
        m[1][current_column] = vertex->coordinates().y;
        m[2][current_column] = vertex->coordinates().z;
        current_column++;
    }
    MatrixMath matrix_math;

    return -matrix_math.ComputeDeterminant(m);
}

Float TetrahedronSolver::ComputeC(UInt vertex_local_index,
                                  const std::shared_ptr<FiniteElement> finite_element){
    const UInt n = 3;
    Matrix m(n, n);
    UInt current_column = 0;
    for(UInt i = 0; i < finite_element->GetVertexCount(); i++){
        if(i == vertex_local_index)
            continue;

        const auto& vertex = finite_element->vertices()[i];
        m[0][current_column] = vertex->coordinates().x;
        m[1][current_column] = 1;
        m[2][current_column] = vertex->coordinates().z;
        current_column++;
    }
    MatrixMath matrix_math;

    return matrix_math.ComputeDeterminant(m);
}

Float TetrahedronSolver::ComputeD(UInt vertex_local_index,
                                  const std::shared_ptr<FiniteElement> finite_element){
    const UInt n = 3;
    Matrix m(n, n);
    UInt current_column = 0;
    for(UInt i = 0; i < finite_element->GetVertexCount(); i++){
        if(i == vertex_local_index)
            continue;

        const auto& vertex = finite_element->vertices()[i];
        m[0][current_column] = vertex->coordinates().x;
        m[1][current_column] = vertex->coordinates().y;
        m[2][current_column] = 1;
        current_column++;
    }
    MatrixMath matrix_math;

    return matrix_math.ComputeDeterminant(m);
}

Float TetrahedronSolver::ComputeVolume(const std::shared_ptr<FiniteElement> finite_element){
    const UInt n = 4;
    Matrix m(n, n);
    UInt current_column = 0;
    for(UInt i = 0; i < finite_element->GetVertexCount(); i++){
        const auto& vertex = finite_element->vertices()[i];

        m[0][current_column] = 1;
        m[1][current_column] = vertex->coordinates().x;
        m[2][current_column] = vertex->coordinates().y;
        m[3][current_column] = vertex->coordinates().z;
        current_column++;
    }
    MatrixMath matrix_math;

    Float volume = matrix_math.ComputeDeterminant(m) / (Float)6.0;
    return volume;
}

Float TetrahedronSolver::ComputeShapeFunctionValue(Float a, Float b,
                                                   Float c, Float d,
                                                   Float volume,
                                                   const Vector3& x){
    return (a + b*x.x + c*x.y + d*x.z) / 6*volume;
}

void TetrahedronSolver::AssemblyGeometryMatrix(Matrix& B, UInt column_offset,
                                               Float b, Float c, Float d){
    B[0][0 + column_offset] = b;
    B[1][1 + column_offset] = c;
    B[2][2 + column_offset] = d;
    B[3][0 + column_offset] = c;
    B[3][1 + column_offset] = b;
    B[4][1 + column_offset] = d;
    B[4][2 + column_offset] = c;
    B[5][0 + column_offset] = d;
    B[5][2 + column_offset] = b;
}

}