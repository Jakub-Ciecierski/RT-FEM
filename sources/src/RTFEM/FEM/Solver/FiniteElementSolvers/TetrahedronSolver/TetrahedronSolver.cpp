#include "RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver/TetrahedronSolver.h"

#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/FEM/Meshing/TriangleMesh.h>
#include <RTFEM/FEM/Material.h>

#include <Eigen/Geometry>
#include <iostream>

namespace rtfem {

template<class T>
FiniteElementSolverData<T> TetrahedronSolver<T>::Solve(
    std::shared_ptr<FiniteElement<T>> finite_element,
    const std::vector<std::shared_ptr<Vertex<T>>> &vertices,
    std::vector<TriangleFace<T>> &triangle_faces,
    const Eigen::Vector3<T> &body_force,
    const Material<T>& material){
    auto vertex_index1 = finite_element->vertices_indices()[0];
    auto vertex_index2 = finite_element->vertices_indices()[1];
    auto vertex_index3 = finite_element->vertices_indices()[2];
    auto vertex_index4 = finite_element->vertices_indices()[3];

    auto vertex1 = vertices[vertex_index1];
    auto vertex2 = vertices[vertex_index2];
    auto vertex3 = vertices[vertex_index3];
    auto vertex4 = vertices[vertex_index4];

    auto& triangle_face1
        = triangle_faces[finite_element->faces_indices()[0]];
    auto& triangle_face2
        = triangle_faces[finite_element->faces_indices()[1]];
    auto& triangle_face3
        = triangle_faces[finite_element->faces_indices()[2]];
    auto& triangle_face4
        = triangle_faces[finite_element->faces_indices()[3]];

    auto ordered_triangle_faces = FetchOrderedTriangleFaces(
        vertex_index1, vertex_index2, vertex_index3, vertex_index4,
        triangle_face1, triangle_face2, triangle_face3, triangle_face4);

    TractionForces<T> traction_forces = FetchTractionForce(
        ordered_triangle_faces);

    auto edges = ComputeEdgesCache(*vertex1, *vertex2, *vertex3, *vertex4);
    auto faces_normal = ComputeFacesNormal(edges);
    auto faces_area = ComputeFacesArea(faces_normal);
    NormalizeFacesNormal(faces_normal);
    SaveTriangleFaceData(ordered_triangle_faces, faces_area, faces_normal);

    auto shape_function_coefficients =
        ComputeShapeFunctionCoefficients(*vertex1, *vertex2,
                                         *vertex3, *vertex4,
                                         edges);

    FiniteElementSolverData<T> data;
    data.volume = ComputeVolume(shape_function_coefficients);
    finite_element->SetVolume(data.volume);
    data.geometry_matrix =
        ComputeGeometryMatrix(shape_function_coefficients, data.volume);
    data.force_vector = ComputeForceVector(shape_function_coefficients,
                                           data.volume, faces_area,
                                           body_force, traction_forces,
                                           material);

    return data;
}

template<class T>
Eigen::Matrix<T, TETRAHEDRON_JACOBIAN_MATRIX_N, TETRAHEDRON_JACOBIAN_MATRIX_N>
TetrahedronSolver<T>::SolveJacobianInverse(
    std::shared_ptr<FiniteElement<T>> finite_element,
    const std::vector<std::shared_ptr<Vertex<T>>> &vertices) {
    auto v1 = vertices[finite_element->vertices_indices()[0]];
    auto v2 = vertices[finite_element->vertices_indices()[1]];
    auto v3 = vertices[finite_element->vertices_indices()[2]];
    auto v4 = vertices[finite_element->vertices_indices()[3]];

    auto edges = ComputeEdgesCache(*v1, *v2, *v3, *v4);

    auto shape_function_coefficients =
        ComputeShapeFunctionCoefficients(*v1, *v2, *v3, *v4, edges);
    auto volume = ComputeVolume(shape_function_coefficients);

    return AssemblyJacobianInverse(shape_function_coefficients, volume);
}

template<class T>
std::vector<TriangleFace<T>*> TetrahedronSolver<T>::FetchOrderedTriangleFaces(
    unsigned int vertex_index1,
    unsigned int vertex_index2,
    unsigned int vertex_index3,
    unsigned int vertex_index4,
    TriangleFace<T>& triangle_face1,
    TriangleFace<T>& triangle_face2,
    TriangleFace<T>& triangle_face3,
    TriangleFace<T>& triangle_face4){
    std::vector<TriangleFace<T>*> ordered_triangle_faces;

    const auto triangle_first = TriangleFace<T>{vertex_index2,
                                                vertex_index3,
                                                vertex_index4};
    const auto triangle_second = TriangleFace<T>{vertex_index1,
                                                 vertex_index3,
                                                 vertex_index4};
    const auto triangle_third = TriangleFace<T>{vertex_index1,
                                                vertex_index2,
                                                vertex_index4};
    const auto triangle_fourth = TriangleFace<T>{vertex_index1,
                                                 vertex_index2,
                                                 vertex_index3};
    // 0
    if(triangle_face1 == triangle_first){
        ordered_triangle_faces.push_back(&triangle_face1);
    }
    else if(triangle_face2 == triangle_first){
        ordered_triangle_faces.push_back(&triangle_face2);
    }
    else if(triangle_face3 == triangle_first){
        ordered_triangle_faces.push_back(&triangle_face3);
    }
    else if(triangle_face4 == triangle_first){
        ordered_triangle_faces.push_back(&triangle_face4);
    }

    // 1
    if(triangle_face1 == triangle_second){
        ordered_triangle_faces.push_back(&triangle_face1);
    }
    else if(triangle_face2 == triangle_second){
        ordered_triangle_faces.push_back(&triangle_face2);
    }
    else if(triangle_face3 == triangle_second){
        ordered_triangle_faces.push_back(&triangle_face3);
    }
    else if(triangle_face4 == triangle_second){
        ordered_triangle_faces.push_back(&triangle_face4);
    }

    // 2
    if(triangle_face1 == triangle_third){
        ordered_triangle_faces.push_back(&triangle_face1);
    }
    else if(triangle_face2 == triangle_third){
        ordered_triangle_faces.push_back(&triangle_face2);
    }
    else if(triangle_face3 == triangle_third){
        ordered_triangle_faces.push_back(&triangle_face3);
    }
    else if(triangle_face4 == triangle_third){
        ordered_triangle_faces.push_back(&triangle_face4);
    }

    // 3
    if(triangle_face1 == triangle_fourth){
        ordered_triangle_faces.push_back(&triangle_face1);
    }
    else if(triangle_face2 == triangle_fourth){
        ordered_triangle_faces.push_back(&triangle_face2);
    }
    else if(triangle_face3 == triangle_fourth){
        ordered_triangle_faces.push_back(&triangle_face3);
    }
    else if(triangle_face4 == triangle_fourth){
        ordered_triangle_faces.push_back(&triangle_face4);
    }

    return ordered_triangle_faces;
}

template<class T>
TractionForces<T> TetrahedronSolver<T>::FetchTractionForce(
    const std::vector<TriangleFace<T>*>& ordered_triangle_faces){

    TractionForces<T> traction_forces;
    traction_forces.force_face1 = ordered_triangle_faces[0]->traction_force;
    traction_forces.force_face2 = ordered_triangle_faces[1]->traction_force;
    traction_forces.force_face3 = ordered_triangle_faces[2]->traction_force;
    traction_forces.force_face4 = ordered_triangle_faces[3]->traction_force;

    return traction_forces;
}

template<class T>
Edges<T> TetrahedronSolver<T>::ComputeEdgesCache(const Vertex<T> &v1,
                                                 const Vertex<T> &v2,
                                                 const Vertex<T> &v3,
                                                 const Vertex<T> &v4) {
    Edges<T> edges;

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

template<class T>
FacesNormal<T> TetrahedronSolver<T>::ComputeFacesNormal(const Edges<T> &edges) {
    FacesNormal<T> faces_normals;

    faces_normals.normal1 =
        Eigen::Vector<T, 3>(edges.x32, edges.y32, edges.z32).cross(
            Eigen::Vector<T, 3>(edges.x42, edges.y42, edges.z42));
    faces_normals.normal2 =
        Eigen::Vector<T, 3>(edges.x43, edges.y43, edges.z43).cross(
            Eigen::Vector<T, 3>(edges.x13, edges.y13, edges.z13));
    faces_normals.normal3 =
        Eigen::Vector<T, 3>(edges.x14, edges.y14, edges.z14).cross(
            Eigen::Vector<T, 3>(edges.x24, edges.y24, edges.z24));
    faces_normals.normal4 =
        Eigen::Vector<T, 3>(edges.x21, edges.y21, edges.z21).cross(
            Eigen::Vector<T, 3>(edges.x31, edges.y31, edges.z31));

    return faces_normals;
}


template<class T>
FacesArea<T> TetrahedronSolver<T>::ComputeFacesArea(
    const FacesNormal<T>& face_normals) {
    FacesArea<T> faces_area;
    T half = 0.5;
    faces_area.area1 = half * face_normals.normal1.norm();
    faces_area.area2 = half * face_normals.normal2.norm();
    faces_area.area3 = half * face_normals.normal3.norm();
    faces_area.area4 = half * face_normals.normal4.norm();

    return faces_area;
}

template<class T>
void TetrahedronSolver<T>::NormalizeFacesNormal(FacesNormal<T> &normals){
    normals.normal1.normalize();
    normals.normal2.normalize();
    normals.normal3.normalize();
    normals.normal4.normalize();
}

template<class T>
void TetrahedronSolver<T>::SaveTriangleFaceData(
    std::vector<TriangleFace<T>*>& ordered_triangle_faces,
    const FacesArea<T>& faces_area,
    const FacesNormal<T>& faces_normal){
    ordered_triangle_faces[0]->area = faces_area.area1;
    ordered_triangle_faces[1]->area = faces_area.area2;
    ordered_triangle_faces[2]->area = faces_area.area3;
    ordered_triangle_faces[3]->area = faces_area.area4;

    ordered_triangle_faces[0]->normal = faces_normal.normal1;
    ordered_triangle_faces[1]->normal = faces_normal.normal2;
    ordered_triangle_faces[2]->normal = faces_normal.normal3;
    ordered_triangle_faces[3]->normal = faces_normal.normal4;
}

template<class T>
TetrahedronShapeFunctionCoefficients<T>
TetrahedronSolver<T>::ComputeShapeFunctionCoefficients(
    const Vertex<T> &v1, const Vertex<T> &v2,
    const Vertex<T> &v3, const Vertex<T> &v4,
    const Edges<T> &edges) {
    TetrahedronShapeFunctionCoefficients<T> coefficients;

    auto V0i = ComputeAi(v1, v2, v3, v4);
    coefficients.A1 = V0i.x();
    coefficients.A2 = V0i.y();
    coefficients.A3 = V0i.z();
    coefficients.A4 = V0i.w();

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

template<class T>
Eigen::Vector4<T> TetrahedronSolver<T>::ComputeAi(const Vertex<T> &v1,
                                                  const Vertex<T> &v2,
                                                  const Vertex<T> &v3,
                                                  const Vertex<T> &v4) {
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

    Eigen::Vector4<T> V0i(V01, V02, V03, V04);

    return V0i;
}

template<class T>
T TetrahedronSolver<T>::ComputeVolume(const TetrahedronShapeFunctionCoefficients<
    T> &coefficients) {
    auto volume = (coefficients.A1 / 6.0)
        + (coefficients.A2 / 6.0)
        + (coefficients.A3 / 6.0)
        + (coefficients.A4 / 6.0);
    if (volume == 0) {
        throw std::invalid_argument(
            "TetrahedronSolver<T>::Solve: Element with 0 volume");
    }
    if (volume < 0) {
        throw std::invalid_argument(
            "TetrahedronSolver<T>::Solve: Element with negative volume");
    }

    return volume;
}

template<class T>
Eigen::Matrix<T, TETRAHEDRON_GEOMETRIC_MATRIX_N, TETRAHEDRON_GEOMETRIC_MATRIX_M>
TetrahedronSolver<T>::ComputeGeometryMatrix(const TetrahedronShapeFunctionCoefficients<
    T> &coefficients, T volume) {
    Eigen::Matrix<T,
                  TETRAHEDRON_GEOMETRIC_MATRIX_N,
                  TETRAHEDRON_GEOMETRIC_MATRIX_M> geometry_matrix =
        Eigen::Matrix<T,
                      TETRAHEDRON_GEOMETRIC_MATRIX_N,
                      TETRAHEDRON_GEOMETRIC_MATRIX_M>::Zero();

    AssemblyGeometryMatrix(geometry_matrix,
                           0,
                           coefficients.B1,
                           coefficients.C1,
                           coefficients.D1);
    AssemblyGeometryMatrix(geometry_matrix,
                           3,
                           coefficients.B2,
                           coefficients.C2,
                           coefficients.D2);
    AssemblyGeometryMatrix(geometry_matrix,
                           6,
                           coefficients.B3,
                           coefficients.C3,
                           coefficients.D3);
    AssemblyGeometryMatrix(geometry_matrix,
                           9,
                           coefficients.B4,
                           coefficients.C4,
                           coefficients.D4);
    geometry_matrix = geometry_matrix / (6 * volume);

    return geometry_matrix;
}

template<class T>
void TetrahedronSolver<T>::AssemblyGeometryMatrix(
    Eigen::Matrix<T,
                  TETRAHEDRON_GEOMETRIC_MATRIX_N,
                  TETRAHEDRON_GEOMETRIC_MATRIX_M> &B,
    unsigned int column_offset,
    T b,
    T c,
    T d) {
    B(0, 0 + column_offset) = b;
    B(1, 1 + column_offset) = c;
    B(2, 2 + column_offset) = d;
    B(3, 0 + column_offset) = c;
    B(3, 1 + column_offset) = b;
    B(4, 1 + column_offset) = d;
    B(4, 2 + column_offset) = c;
    B(5, 0 + column_offset) = d;
    B(5, 2 + column_offset) = b;
}

template<class T>
Eigen::Vector<T, TETRAHEDRON_FORCE_VECTOR_N>
TetrahedronSolver<T>::ComputeForceVector(
        const TetrahedronShapeFunctionCoefficients<T> &shape_function_coefficients,
        T volume,
        const FacesArea<T> &faces_area,
        const Eigen::Vector3<T> &body_force,
        const TractionForces<T>
        &traction_force,
        const Material<T> &material) {

    auto consistent_body_force = ComputeBodyForceVector(volume, body_force,
                                                        material);
    auto consistent_traction_force =
            ComputeTractionForceVector(shape_function_coefficients,
                                       faces_area,
                                       traction_force);

    return consistent_body_force + consistent_traction_force;
}

template<class T>
Eigen::Vector<T, TETRAHEDRON_FORCE_VECTOR_N>
TetrahedronSolver<T>::ComputeBodyForceVector(
    T volume,
    const Eigen::Vector3<T> &body_force,
    const Material<T>& material) {
    Eigen::Vector<T, TETRAHEDRON_FORCE_VECTOR_N> consistent_body_force;

    consistent_body_force(0) = body_force.x();
    consistent_body_force(1) = body_force.y();
    consistent_body_force(2) = body_force.z();
    consistent_body_force(3) = body_force.x();
    consistent_body_force(4) = body_force.y();
    consistent_body_force(5) = body_force.z();
    consistent_body_force(6) = body_force.x();
    consistent_body_force(7) = body_force.y();
    consistent_body_force(8) = body_force.z();
    consistent_body_force(9) = body_force.x();
    consistent_body_force(10) = body_force.y();
    consistent_body_force(11) = body_force.z();

    consistent_body_force *= material.density * (volume / 4.0);

    return consistent_body_force;
}

template<class T>
Eigen::Vector<T, TETRAHEDRON_FORCE_VECTOR_N>
TetrahedronSolver<T>::ComputeTractionForceVector(
    const TetrahedronShapeFunctionCoefficients<T> &shape_function_coefficients,
    const FacesArea<T> &faces_area,
    const TractionForces<T> &traction_force) {
    auto &coeff = shape_function_coefficients;

    auto traction_force_face1 =
        ComputeTractionForceVectorFace(0,
                                       traction_force.force_face1,
                                       faces_area.area1,
                                       coeff.B1,
                                       coeff.C1,
                                       coeff.D1);
    auto traction_force_face2 =
        ComputeTractionForceVectorFace(1,
                                       traction_force.force_face2,
                                       faces_area.area2,
                                       coeff.B2,
                                       coeff.C2,
                                       coeff.D2);
    auto traction_force_face3 =
        ComputeTractionForceVectorFace(2,
                                       traction_force.force_face3,
                                       faces_area.area3,
                                       coeff.B3,
                                       coeff.C3,
                                       coeff.D3);
    auto traction_force_face4 =
        ComputeTractionForceVectorFace(3,
                                       traction_force.force_face4,
                                       faces_area.area4,
                                       coeff.B4,
                                       coeff.C4,
                                       coeff.D4);

    return traction_force_face1 + traction_force_face2 + traction_force_face3
        + traction_force_face4;
}

template<class T>
Eigen::Vector<T, TETRAHEDRON_FORCE_VECTOR_N>
TetrahedronSolver<T>::ComputeTractionForceVectorFace(unsigned int face_index,
                                                     T traction_force,
                                                     T area,
                                                     T B,
                                                     T C,
                                                     T D) {
    auto S = std::sqrt((B * B) + (C * C) + (D * D));
    auto B_normal = -(B / S);
    auto C_normal = -(C / S);
    auto D_normal = -(D / S);

    auto magnitude = (1.0 / 3.0) * traction_force * area;
    Eigen::Vector<T, TETRAHEDRON_FORCE_VECTOR_N> traction_force_vector_face;

    traction_force_vector_face(0) = B_normal;
    traction_force_vector_face(1) = C_normal;
    traction_force_vector_face(2) = D_normal;
    traction_force_vector_face(3) = B_normal;
    traction_force_vector_face(4) = C_normal;
    traction_force_vector_face(5) = D_normal;
    traction_force_vector_face(6) = B_normal;
    traction_force_vector_face(7) = C_normal;
    traction_force_vector_face(8) = D_normal;
    traction_force_vector_face(9) = B_normal;
    traction_force_vector_face(10) = C_normal;
    traction_force_vector_face(11) = D_normal;

    traction_force_vector_face *= magnitude;

    traction_force_vector_face(face_index + 0) = 0;
    traction_force_vector_face(face_index + 1) = 0;
    traction_force_vector_face(face_index + 2) = 0;

    return traction_force_vector_face;
}

template<class T>
Eigen::Matrix<T, TETRAHEDRON_JACOBIAN_MATRIX_N, TETRAHEDRON_JACOBIAN_MATRIX_N>
TetrahedronSolver<T>::AssemblyJacobianInverse(const TetrahedronShapeFunctionCoefficients<
    T> &coefficients,
                                              T volume) {
    Eigen::Matrix<T,
                  TETRAHEDRON_JACOBIAN_MATRIX_N,
                  TETRAHEDRON_JACOBIAN_MATRIX_N> jacobian_inverse;
    jacobian_inverse(0, 0) = coefficients.A1;
    jacobian_inverse(1, 0) = coefficients.A2;
    jacobian_inverse(2, 0) = coefficients.A3;
    jacobian_inverse(3, 0) = coefficients.A4;

    jacobian_inverse(0, 1) = coefficients.B1;
    jacobian_inverse(1, 1) = coefficients.B2;
    jacobian_inverse(2, 1) = coefficients.B3;
    jacobian_inverse(3, 1) = coefficients.B4;

    jacobian_inverse(0, 2) = coefficients.C1;
    jacobian_inverse(1, 2) = coefficients.C2;
    jacobian_inverse(2, 2) = coefficients.C3;
    jacobian_inverse(3, 2) = coefficients.C4;

    jacobian_inverse(0, 3) = coefficients.D1;
    jacobian_inverse(1, 3) = coefficients.D2;
    jacobian_inverse(2, 3) = coefficients.D3;
    jacobian_inverse(3, 3) = coefficients.D4;

    jacobian_inverse = jacobian_inverse / (6 * volume);

    return jacobian_inverse;
}

template
class TetrahedronSolver<double>;
template
class TetrahedronSolver<float>;

}