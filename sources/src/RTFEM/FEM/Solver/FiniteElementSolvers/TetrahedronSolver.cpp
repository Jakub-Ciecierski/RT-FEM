#include "RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver.h"

#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/Math/MatrixMath.h>
#include <RTFEM/DataStructure/Vector4.h>
#include <cmath>

namespace rtfem {

TetrahedronSolver::TetrahedronSolver() {}

TetrahedronSolver::~TetrahedronSolver() {}

FiniteElementSolverData TetrahedronSolver::Solve(std::shared_ptr<FiniteElement> finite_element){
    return Solve(finite_element, Vector3(0,0,0), TractionForce{});
}

FiniteElementSolverData TetrahedronSolver::Solve(
        std::shared_ptr<FiniteElement> finite_element,
        const Vector3& body_force,
        const TractionForce& traction_force){
    auto v1 = finite_element->vertices()[0];
    auto v2 = finite_element->vertices()[1];
    auto v3 = finite_element->vertices()[2];
    auto v4 = finite_element->vertices()[3];

    ComputeEdgesCache(*v1, *v2, *v3, *v4);

    auto shape_function_coefficients = ComputeShapeFunctionCoefficients(*v1, *v2 ,*v3, *v4);

    FiniteElementSolverData data;
    data.volume = ComputeVolume(shape_function_coefficients);
    data.geometry_matrix = ComputeGeometryMatrix(shape_function_coefficients, data.volume);
    data.force_vector = ComputeForceVector(shape_function_coefficients,
                                           data.volume,
                                           body_force, traction_force);

    return data;
}

Matrix TetrahedronSolver::SolveJacobianInverse(std::shared_ptr<FiniteElement> finite_element){
    auto v1 = finite_element->vertices()[0];
    auto v2 = finite_element->vertices()[1];
    auto v3 = finite_element->vertices()[2];
    auto v4 = finite_element->vertices()[3];

    ComputeEdgesCache(*v1, *v2, *v3, *v4);

    auto shape_function_coefficients = ComputeShapeFunctionCoefficients(*v1, *v2 ,*v3, *v4);
    auto volume = ComputeVolume(shape_function_coefficients);

    return AssemblyJacobianInverse(shape_function_coefficients, volume);
}

void TetrahedronSolver::ComputeEdgesCache(const Vertex &v1, const Vertex &v2,
                                          const Vertex &v3, const Vertex &v4) {
    edges_cache_ = EdgesCache{};

    edges_cache_.x32 = v3.x() - v2.x();
    edges_cache_.x34 = v3.x() - v4.x();
    edges_cache_.x43 = v4.x() - v3.x();
    edges_cache_.x14 = v1.x() - v4.x();
    edges_cache_.x21 = v2.x() - v1.x();
    edges_cache_.x31 = v3.x() - v1.x();
    edges_cache_.x24 = v2.x() - v4.x();
    edges_cache_.x42 = v4.x() - v2.x();
    edges_cache_.x13 = v1.x() - v3.x();
    edges_cache_.x12 = v1.x() - v2.x();
    edges_cache_.z43 = v4.z() - v3.z();
    edges_cache_.z31 = v3.z() - v1.z();
    edges_cache_.z32 = v3.z() - v2.z();
    edges_cache_.z24 = v2.z() - v4.z();
    edges_cache_.z34 = v3.z() - v4.z();
    edges_cache_.z13 = v1.z() - v3.z();
    edges_cache_.z14 = v1.z() - v4.z();
    edges_cache_.z21 = v2.z() - v1.z();
    edges_cache_.z42 = v4.z() - v2.z();
    edges_cache_.z12 = v1.z() - v2.z();
    edges_cache_.y42 = v4.y() - v2.y();
    edges_cache_.y31 = v3.y() - v1.y();
    edges_cache_.y24 = v2.y() - v4.y();
    edges_cache_.y13 = v1.y() - v3.y();
    edges_cache_.y32 = v3.y() - v2.y();
    edges_cache_.y34 = v3.y() - v4.y();
    edges_cache_.y14 = v1.y() - v4.y();
    edges_cache_.y12 = v1.y() - v2.y();
    edges_cache_.y43 = v4.y() - v3.y();
    edges_cache_.y21 = v2.y() - v1.y();
}

TetrahedronShapeFunctionCoefficients TetrahedronSolver::ComputeShapeFunctionCoefficients(const Vertex &v1,
                                                                                         const Vertex &v2,
                                                                                         const Vertex &v3,
                                                                                         const Vertex &v4) {
    TetrahedronShapeFunctionCoefficients coefficients;

    auto V0i = ComputeAi(v1, v2, v3, v4);
    coefficients.A1 = V0i.x;
    coefficients.A2 = V0i.y;
    coefficients.A3 = V0i.z;
    coefficients.A4 = V0i.w;

    coefficients.B1 = (edges_cache_.y42 * edges_cache_.z32) - (edges_cache_.y32 * edges_cache_.z42);
    coefficients.B2 = (edges_cache_.y31 * edges_cache_.z43) - (edges_cache_.y34 * edges_cache_.z13);
    coefficients.B3 = (edges_cache_.y24 * edges_cache_.z14) - (edges_cache_.y14 * edges_cache_.z24);
    coefficients.B4 = (edges_cache_.y13 * edges_cache_.z21) - (edges_cache_.y12 * edges_cache_.z31);

    coefficients.C1 = (edges_cache_.x32 * edges_cache_.z42) - (edges_cache_.x42 * edges_cache_.z32);
    coefficients.C2 = (edges_cache_.x43 * edges_cache_.z31) - (edges_cache_.x13 * edges_cache_.z34);
    coefficients.C3 = (edges_cache_.x14 * edges_cache_.z24) - (edges_cache_.x24 * edges_cache_.z14);
    coefficients.C4 = (edges_cache_.x21 * edges_cache_.z13) - (edges_cache_.x31 * edges_cache_.z12);

    coefficients.D1 = (edges_cache_.x42 * edges_cache_.y32) - (edges_cache_.x32 * edges_cache_.y42);
    coefficients.D2 = (edges_cache_.x31 * edges_cache_.y43) - (edges_cache_.x34 * edges_cache_.y13);
    coefficients.D3 = (edges_cache_.x24 * edges_cache_.y14) - (edges_cache_.x14 * edges_cache_.y24);
    coefficients.D4 = (edges_cache_.x13 * edges_cache_.y21) - (edges_cache_.x12 * edges_cache_.y31);

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
                                             Float volume,
                                             const Vector3& body_force,
                                             const TractionForce& traction_force){
    auto consistent_body_force = ComputeBodyForceVector(volume, body_force);
    auto consistent_traction_force = ComputeTractionForceVector(shape_function_coefficients, traction_force);

    return consistent_body_force + consistent_traction_force;
}

Matrix TetrahedronSolver::ComputeBodyForceVector(
        Float volume,
        const Vector3& body_force){
    Matrix consistent_body_force(12,1);

    consistent_body_force[0][0] = body_force.x;
    consistent_body_force[1][0] = body_force.y;
    consistent_body_force[2][0] = body_force.z;
    consistent_body_force[3][0] = body_force.x;
    consistent_body_force[4][0] = body_force.y;
    consistent_body_force[5][0] = body_force.z;
    consistent_body_force[6][0] = body_force.x;
    consistent_body_force[7][0] = body_force.y;
    consistent_body_force[8][0] = body_force.z;
    consistent_body_force[9][0] = body_force.x;
    consistent_body_force[10][0] = body_force.y;
    consistent_body_force[11][0] = body_force.z;

    consistent_body_force = consistent_body_force * (volume/4.0);

    return consistent_body_force;
}

Matrix TetrahedronSolver::ComputeTractionForceVector(
        const TetrahedronShapeFunctionCoefficients& shape_function_coefficients,
        const TractionForce& traction_force){
    auto& coeff = shape_function_coefficients;
    auto S1 = std::sqrt((coeff.A1 * coeff.A1) + (coeff.B1 * coeff.B1) + (coeff.C1 * coeff.C1));
    auto A1_normal = coeff.A1 / S1;
    auto B1_normal = coeff.B1 / S1;
    auto C1_normal = coeff.C1 / S1;


    return Matrix(12,1);
}

}