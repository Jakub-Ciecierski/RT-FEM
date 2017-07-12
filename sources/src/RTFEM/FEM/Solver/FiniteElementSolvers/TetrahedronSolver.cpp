#include "RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver.h"

#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/Math/MatrixMath.h>
#include <RTFEM/DataStructure/Vector4.h>

namespace rtfem {

TetrahedronSolver::TetrahedronSolver() {}

TetrahedronSolver::~TetrahedronSolver() {}

FiniteElementSolverData TetrahedronSolver::Solve(std::shared_ptr<FiniteElement> finite_element){
    return Solve(finite_element, Vector3(0,0,0), Vector3(0,0,0));
}

FiniteElementSolverData TetrahedronSolver::Solve(
        std::shared_ptr<FiniteElement> finite_element,
        const Vector3& body_force,
        const Vector3& traction_force){
    auto v1 = finite_element->vertices()[0];
    auto v2 = finite_element->vertices()[1];
    auto v3 = finite_element->vertices()[2];
    auto v4 = finite_element->vertices()[3];

    auto shape_function_coefficients = ComputeShapeFunctionCoefficients(*v1, *v2 ,*v3, *v4);

    FiniteElementSolverData data;
    data.volume = ComputeVolume(shape_function_coefficients);
    data.geometry_matrix = ComputeGeometryMatrix(shape_function_coefficients, data.volume);
    data.force_vector = ComputeForceVector(shape_function_coefficients, body_force, traction_force);

    return data;
}

Matrix TetrahedronSolver::SolveJacobianInverse(std::shared_ptr<FiniteElement> finite_element){
    auto v1 = finite_element->vertices()[0];
    auto v2 = finite_element->vertices()[1];
    auto v3 = finite_element->vertices()[2];
    auto v4 = finite_element->vertices()[3];

    auto shape_function_coefficients = ComputeShapeFunctionCoefficients(*v1, *v2 ,*v3, *v4);
    auto volume = ComputeVolume(shape_function_coefficients);

    return AssemblyJacobianInverse(shape_function_coefficients, volume);
}

TetrahedronShapeFunctionCoefficients TetrahedronSolver::ComputeShapeFunctionCoefficients(const Vertex &v1,
                                                                              const Vertex &v2,
                                                                              const Vertex &v3,
                                                                              const Vertex &v4) {
    auto x32 = v3.x() - v2.x();
    auto x34 = v3.x() - v4.x();
    auto x43 = v4.x() - v3.x();
    auto x14 = v1.x() - v4.x();
    auto x21 = v2.x() - v1.x();
    auto x31 = v3.x() - v1.x();
    auto x24 = v2.x() - v4.x();
    auto x42 = v4.x() - v2.x();
    auto x13 = v1.x() - v3.x();
    auto x12 = v1.x() - v2.x();

    auto z43 = v4.z() - v3.z();
    auto z31 = v3.z() - v1.z();
    auto z32 = v3.z() - v2.z();
    auto z24 = v2.z() - v4.z();
    auto z34 = v3.z() - v4.z();
    auto z13 = v1.z() - v3.z();
    auto z14 = v1.z() - v4.z();
    auto z21 = v2.z() - v1.z();
    auto z42 = v4.z() - v2.z();
    auto z12 = v1.z() - v2.z();

    auto y42 = v4.y() - v2.y();
    auto y31 = v3.y() - v1.y();
    auto y24 = v2.y() - v4.y();
    auto y13 = v1.y() - v3.y();
    auto y32 = v3.y() - v2.y();
    auto y34 = v3.y() - v4.y();
    auto y14 = v1.y() - v4.y();
    auto y12 = v1.y() - v2.y();
    auto y43 = v4.y() - v3.y();
    auto y21 = v2.y() - v1.y();

    TetrahedronShapeFunctionCoefficients coefficients;

    auto V0i = ComputeAi(v1, v2, v3, v4);
    coefficients.A1 = V0i.x;
    coefficients.A2 = V0i.y;
    coefficients.A3 = V0i.z;
    coefficients.A4 = V0i.w;

    coefficients.B1 = (y42 * z32) - (y32 * z42);
    coefficients.B2 = (y31 * z43) - (y34 * z13);
    coefficients.B3 = (y24 * z14) - (y14 * z24);
    coefficients.B4 = (y13 * z21) - (y12 * z31);

    coefficients.C1 = (x32 * z42) - (x42 * z32);
    coefficients.C2 = (x43 * z31) - (x13 * z34);
    coefficients.C3 = (x14 * z24) - (x24 * z14);
    coefficients.C4 = (x21 * z13) - (x31 * z12);

    coefficients.D1 = (x42 * y32) - (x32 * y42);
    coefficients.D2 = (x31 * y43) - (x34 * y13);
    coefficients.D3 = (x24 * y14) - (x14 * y24);
    coefficients.D4 = (x13 * y21) - (x12 * y31);

    return coefficients;
}

Vector4 TetrahedronSolver::ComputeAi(const Vertex &v1,
                                     const Vertex &v2,
                                     const Vertex &v3,
                                     const Vertex &v4) {
    auto V01 =
            v2.x() * ((v3.y() * v4.z()) - (v4.y() * v3.z()))
            + v3.x() * ((v4.y() * v2.z()) - (v2.y() * v4.z()))
            + v4.x() * ((v2.y() * v3.z()) - (v3.y() * v2.z()));
    auto V02 =
            v1.x() * ((v4.y() * v3.z()) - (v3.y() * v4.z()))
            + v3.x() * ((v1.y() * v4.z()) - (v4.y() * v1.z()))
            + v4.x() * ((v3.y() * v1.z()) - (v1.y() * v3.z()));
    auto V03 =
            v1.x() * ((v2.y() * v4.z()) - (v4.y() * v2.z()))
            + v2.x() * ((v4.y() * v1.z()) - (v1.y() * v4.z()))
            + v4.x() * ((v1.y() * v2.z()) - (v2.y() * v1.z()));
    auto V04 =
            v1.x() * ((v3.y() * v2.z()) - (v2.y() * v3.z()))
            + v2.x() * ((v1.y() * v3.z()) - (v3.y() * v1.z()))
            + v3.x() * ((v2.y() * v1.z()) - (v1.y() * v2.z()));

    Vector4 V0i(V01, V02, V03, V04);

    return V0i;
}

Matrix TetrahedronSolver::ComputeGeometryMatrix(const TetrahedronShapeFunctionCoefficients& coefficients,
                                                Float volume){
    Matrix geometry_matrix(6, 12);
    AssemblyGeometryMatrix(geometry_matrix, 0, coefficients.B1, coefficients.C1, coefficients.D1);
    AssemblyGeometryMatrix(geometry_matrix, 3, coefficients.B2, coefficients.C2, coefficients.D2);
    AssemblyGeometryMatrix(geometry_matrix, 6, coefficients.B3, coefficients.C3, coefficients.D3);
    AssemblyGeometryMatrix(geometry_matrix, 9, coefficients.B4, coefficients.C4, coefficients.D4);
    geometry_matrix = geometry_matrix / (6*volume);

    return geometry_matrix;
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

Float TetrahedronSolver::ComputeVolume(const TetrahedronShapeFunctionCoefficients& coefficients){
    auto volume = (coefficients.A1 / 6.0)
                  + (coefficients.A2 / 6.0)
                  + (coefficients.A3 / 6.0)
                  + (coefficients.A4 / 6.0);
    if(volume == 0){
        throw std::invalid_argument("TetrahedronSolver::Solve: Element with 0 volume");
    }

    return volume;
}

Matrix TetrahedronSolver::AssemblyJacobianInverse(const TetrahedronShapeFunctionCoefficients& coefficients,
                                                  Float volume){
    Matrix jacobian_inverse(4,4);
    jacobian_inverse[0][0] = coefficients.A1;
    jacobian_inverse[1][0] = coefficients.A2;
    jacobian_inverse[2][0] = coefficients.A3;
    jacobian_inverse[3][0] = coefficients.A4;

    jacobian_inverse[0][1] = coefficients.B1;
    jacobian_inverse[1][1] = coefficients.B2;
    jacobian_inverse[2][1] = coefficients.B3;
    jacobian_inverse[3][1] = coefficients.B4;

    jacobian_inverse[0][2] = coefficients.C1;
    jacobian_inverse[1][2] = coefficients.C2;
    jacobian_inverse[2][2] = coefficients.C3;
    jacobian_inverse[3][2] = coefficients.C4;

    jacobian_inverse[0][3] = coefficients.D1;
    jacobian_inverse[1][3] = coefficients.D2;
    jacobian_inverse[2][3] = coefficients.D3;
    jacobian_inverse[3][3] = coefficients.D4;

    jacobian_inverse = jacobian_inverse / (6 * volume);

    return jacobian_inverse;
}

Matrix TetrahedronSolver::ComputeForceVector(const TetrahedronShapeFunctionCoefficients& shape_function_coefficients,
                                             const Vector3& body_force,
                                             const Vector3& traction_force){

}

}