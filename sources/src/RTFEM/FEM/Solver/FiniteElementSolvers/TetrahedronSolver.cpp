#include "RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver.h"

#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/Math/MatrixMath.h>
#include <RTFEM/DataStructure/Vector4.h>
#include <cmath>
#include <RTFEM/Math/VectorMath.h>

namespace rtfem {

constexpr UInt DOF_COUNT = 12;

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

    auto edges = ComputeEdgesCache(*v1, *v2, *v3, *v4);
    auto faces_area = ComputeFacesArea(edges);
    auto shape_function_coefficients = ComputeShapeFunctionCoefficients(*v1, *v2 ,*v3, *v4, edges);

    FiniteElementSolverData data;
    data.volume = ComputeVolume(shape_function_coefficients);
    data.geometry_matrix = ComputeGeometryMatrix(shape_function_coefficients, data.volume);
    data.force_vector = ComputeForceVector(shape_function_coefficients,
                                           data.volume, faces_area,
                                           body_force, traction_force);

    return data;
}

Matrix TetrahedronSolver::SolveJacobianInverse(std::shared_ptr<FiniteElement> finite_element){
    auto v1 = finite_element->vertices()[0];
    auto v2 = finite_element->vertices()[1];
    auto v3 = finite_element->vertices()[2];
    auto v4 = finite_element->vertices()[3];

    auto edges = ComputeEdgesCache(*v1, *v2, *v3, *v4);

    auto shape_function_coefficients = ComputeShapeFunctionCoefficients(*v1, *v2 ,*v3, *v4, edges);
    auto volume = ComputeVolume(shape_function_coefficients);

    return AssemblyJacobianInverse(shape_function_coefficients, volume);
}

Edges TetrahedronSolver::ComputeEdgesCache(const Vertex &v1, const Vertex &v2,
                                           const Vertex &v3, const Vertex &v4) {
    Edges edges;

    edges.x32 = v3.x() - v2.x();
    edges.x34 = v3.x() - v4.x();
    edges.x43 = v4.x() - v3.x();
    edges.x14 = v1.x() - v4.x();
    edges.x21 = v2.x() - v1.x();
    edges.x31 = v3.x() - v1.x();
    edges.x24 = v2.x() - v4.x();
    edges.x42 = v4.x() - v2.x();
    edges.x13 = v1.x() - v3.x();
    edges.x12 = v1.x() - v2.x();
    edges.z43 = v4.z() - v3.z();
    edges.z31 = v3.z() - v1.z();
    edges.z32 = v3.z() - v2.z();
    edges.z24 = v2.z() - v4.z();
    edges.z34 = v3.z() - v4.z();
    edges.z13 = v1.z() - v3.z();
    edges.z14 = v1.z() - v4.z();
    edges.z21 = v2.z() - v1.z();
    edges.z42 = v4.z() - v2.z();
    edges.z12 = v1.z() - v2.z();
    edges.y42 = v4.y() - v2.y();
    edges.y31 = v3.y() - v1.y();
    edges.y24 = v2.y() - v4.y();
    edges.y13 = v1.y() - v3.y();
    edges.y32 = v3.y() - v2.y();
    edges.y34 = v3.y() - v4.y();
    edges.y14 = v1.y() - v4.y();
    edges.y12 = v1.y() - v2.y();
    edges.y43 = v4.y() - v3.y();
    edges.y21 = v2.y() - v1.y();

    return edges;
}

TetrahedronShapeFunctionCoefficients TetrahedronSolver::ComputeShapeFunctionCoefficients(
        const Vertex &v1, const Vertex &v2,
        const Vertex &v3, const Vertex &v4,
        const Edges& edges) {
    TetrahedronShapeFunctionCoefficients coefficients;

    auto V0i = ComputeAi(v1, v2, v3, v4);
    coefficients.A1 = V0i.x;
    coefficients.A2 = V0i.y;
    coefficients.A3 = V0i.z;
    coefficients.A4 = V0i.w;

    coefficients.B1 = (edges.y42 * edges.z32) - (edges.y32 * edges.z42);
    coefficients.B2 = (edges.y31 * edges.z43) - (edges.y34 * edges.z13);
    coefficients.B3 = (edges.y24 * edges.z14) - (edges.y14 * edges.z24);
    coefficients.B4 = (edges.y13 * edges.z21) - (edges.y12 * edges.z31);

    coefficients.C1 = (edges.x32 * edges.z42) - (edges.x42 * edges.z32);
    coefficients.C2 = (edges.x43 * edges.z31) - (edges.x13 * edges.z34);
    coefficients.C3 = (edges.x14 * edges.z24) - (edges.x24 * edges.z14);
    coefficients.C4 = (edges.x21 * edges.z13) - (edges.x31 * edges.z12);

    coefficients.D1 = (edges.x42 * edges.y32) - (edges.x32 * edges.y42);
    coefficients.D2 = (edges.x31 * edges.y43) - (edges.x34 * edges.y13);
    coefficients.D3 = (edges.x24 * edges.y14) - (edges.x14 * edges.y24);
    coefficients.D4 = (edges.x13 * edges.y21) - (edges.x12 * edges.y31);

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
    Matrix geometry_matrix(6, DOF_COUNT);
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
                                             Float volume, const FacesArea& faces_area,
                                             const Vector3& body_force,
                                             const TractionForce& traction_force){
    auto consistent_body_force = ComputeBodyForceVector(volume, body_force);
    auto consistent_traction_force = ComputeTractionForceVector(shape_function_coefficients,
                                                                faces_area,
                                                                traction_force);

    return consistent_body_force + consistent_traction_force;
}

Matrix TetrahedronSolver::ComputeBodyForceVector(
        Float volume,
        const Vector3& body_force){
    Matrix consistent_body_force(DOF_COUNT,1);

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
        const FacesArea& faces_area,
        const TractionForce& traction_force){
    auto& coeff = shape_function_coefficients;

    auto traction_force_face1 =
            ComputeTractionForceVectorFace(0, traction_force.force_face1, faces_area.area1,
                                           coeff.B1, coeff.C1, coeff.D1);
    auto traction_force_face2 =
            ComputeTractionForceVectorFace(1, traction_force.force_face2, faces_area.area2,
                                           coeff.B2, coeff.C2, coeff.D2);
    auto traction_force_face3 =
            ComputeTractionForceVectorFace(2, traction_force.force_face3, faces_area.area3,
                                           coeff.B3, coeff.C3, coeff.D3);
    auto traction_force_face4 =
            ComputeTractionForceVectorFace(3, traction_force.force_face4, faces_area.area4,
                                           coeff.B4, coeff.C4, coeff.D4);

    return traction_force_face1 + traction_force_face2 + traction_force_face3 + traction_force_face4;
}

Matrix TetrahedronSolver::ComputeTractionForceVectorFace(UInt face_index, Float traction_force,
                                                         Float area,
                                                         Float B, Float C, Float D){
    auto S = std::sqrt((B * B) + (C * C) + (D * D));
    auto B_normal = B / S;
    auto C_normal = C / S;
    auto D_normal = D / S;

    auto magnitude = (1.0 / 3.0) * traction_force * area;

    Matrix traction_force_vector_face(DOF_COUNT, 1);
    traction_force_vector_face[0][0] = B_normal;
    traction_force_vector_face[1][0] = C_normal;
    traction_force_vector_face[2][0] = D_normal;
    traction_force_vector_face[3][0] = B_normal;
    traction_force_vector_face[4][0] = C_normal;
    traction_force_vector_face[5][0] = D_normal;
    traction_force_vector_face[6][0] = B_normal;
    traction_force_vector_face[7][0] = C_normal;
    traction_force_vector_face[8][0] = D_normal;
    traction_force_vector_face[9][0] = B_normal;
    traction_force_vector_face[10][0] = C_normal;
    traction_force_vector_face[11][0] = D_normal;

    traction_force_vector_face = traction_force_vector_face * magnitude;

    traction_force_vector_face[face_index + 0][0] = 0;
    traction_force_vector_face[face_index + 1][0] = 0;
    traction_force_vector_face[face_index + 2][0] = 0;

    return traction_force_vector_face;
}

FacesArea TetrahedronSolver::ComputeFacesArea(const Edges& edges){
    VectorMath vector_math_;
    FacesArea faces_area;

    faces_area.area1 = vector_math_.Magnitude(
            vector_math_.Cross(Vector3(edges.x32, edges.y32, edges.z32),
                               Vector3(edges.x42, edges.y42, edges.z42)));
    faces_area.area2 = vector_math_.Magnitude(
            vector_math_.Cross(Vector3(edges.x43, edges.y43, edges.z43),
                               Vector3(edges.x13, edges.y13, edges.z13)));
    faces_area.area3 = vector_math_.Magnitude(
            vector_math_.Cross(Vector3(edges.x14, edges.y14, edges.z14),
                               Vector3(edges.x24, edges.y24, edges.z24)));
    faces_area.area4 = vector_math_.Magnitude(
            vector_math_.Cross(Vector3(edges.x21, edges.y21, edges.z21),
                               Vector3(edges.x31, edges.y31, edges.z31)));

    return faces_area;

}

}