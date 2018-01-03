#include <RTFEM/FEM/Material.h>
#include <iostream>
#include "RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMFastForceAssembler.h"

#include "RTFEM/FEM/Meshing/TriangleMesh.h"
#include "RTFEM/FEM/FiniteElement.h"
#include "RTFEM/FEM/Vertex.h"

namespace rtfem {

constexpr int DIMENSION_COUNT = 3;

template<class T>
void FEMFastForceAssembler<T>::Assemble(
        FEMGeometry<T> &fem_geometry,
        const Eigen::Vector3<T>& body_force,
        const Material<T>& material,
        Eigen::Vector<T, Eigen::Dynamic> &global_force) {
    int vertex_count = fem_geometry.vertices.size();
    assert(vertex_count * DIMENSION_COUNT == global_force.size());

    ResetGlobalForce(global_force);
    for(unsigned int i = 0; i < fem_geometry.finite_elements.size(); i++){
        auto& finite_element = *fem_geometry.finite_elements[i];

        AddTractionForces(finite_element, fem_geometry, global_force);

        AddBodyForce(finite_element,
                     material, fem_geometry,
                     body_force, global_force);
    }
}

template<class T>
void FEMFastForceAssembler<T>::ResetGlobalForce(
    Eigen::Vector<T, Eigen::Dynamic> &global_force){
    for(unsigned int i = 0; i < global_force.size(); i++){
        global_force(i) = 0;
    }
}

template<class T>
void FEMFastForceAssembler<T>::AddTractionForces(
    FiniteElement<T>& finite_element,
    FEMGeometry<T>& fem_geometry,
    Eigen::Vector<T, Eigen::Dynamic> &global_force){
    AddTractionForce(
        fem_geometry.triangle_faces[finite_element.faces_indices()[0]],
        global_force);
    AddTractionForce(
        fem_geometry.triangle_faces[finite_element.faces_indices()[1]],
        global_force);
    AddTractionForce(
        fem_geometry.triangle_faces[finite_element.faces_indices()[2]],
        global_force);
    AddTractionForce(
        fem_geometry.triangle_faces[finite_element.faces_indices()[3]],
        global_force);
}

template<class T>
void FEMFastForceAssembler<T>::AddTractionForce(
    const TriangleFace<T>& triangle_face,
    Eigen::Vector<T, Eigen::Dynamic> &global_force){
    auto magnitude = (1.0 / 3.0) * triangle_face.traction_force *
        triangle_face.area;
    // Should be cashed!
    auto S = std::sqrt((triangle_face.B * triangle_face.B) +
        (triangle_face.C * triangle_face.C) +
        (triangle_face.D * triangle_face.D));
    if(S == 0){
        return;
    }
    auto B_normal = -(triangle_face.B / S) * magnitude;
    auto C_normal = -(triangle_face.C / S) * magnitude;
    auto D_normal = -(triangle_face.D / S) * magnitude;

    AddTractionForceToVertex(triangle_face.v1 * DIMENSION_COUNT,
                             B_normal, C_normal, D_normal,
                             global_force);
    AddTractionForceToVertex(triangle_face.v2 * DIMENSION_COUNT,
                             B_normal, C_normal, D_normal,
                             global_force);
    AddTractionForceToVertex(triangle_face.v3 * DIMENSION_COUNT,
                             B_normal, C_normal, D_normal,
                             global_force);
}

template<class T>
void FEMFastForceAssembler<T>::AddTractionForceToVertex(
    unsigned int start_index,
    const T& x_value,
    const T& y_value,
    const T& z_value,
    Eigen::Vector<T, Eigen::Dynamic> &global_force){
    global_force(start_index + 0) += x_value;
    global_force(start_index + 1) += y_value;
    global_force(start_index + 2) += z_value;
}

template<class T>
void FEMFastForceAssembler<T>::AddBodyForce(
    FiniteElement<T>& finite_element,
    const Material<T>& material,
    const FEMGeometry<T>& fem_geometry,
    const Eigen::Vector3<T> &body_force,
    Eigen::Vector<T, Eigen::Dynamic> &global_force){
    const auto& vertices_indices = finite_element.vertices_indices();

    auto multiplier = material.density * (finite_element.volume() / 4.0);

    AddBodyForceToVertex(*fem_geometry.vertices[vertices_indices[0]],
                         body_force, multiplier, global_force);
    AddBodyForceToVertex(*fem_geometry.vertices[vertices_indices[1]],
                         body_force, multiplier, global_force);
    AddBodyForceToVertex(*fem_geometry.vertices[vertices_indices[2]],
                         body_force, multiplier, global_force);
    AddBodyForceToVertex(*fem_geometry.vertices[vertices_indices[3]],
                         body_force, multiplier, global_force);
}

template<class T>
void FEMFastForceAssembler<T>::AddBodyForceToVertex(
    const Vertex<T> &vertex,
    const Eigen::Vector3<T> &body_force,
    T multiplier,
    Eigen::Vector<T, Eigen::Dynamic> &global_force){
    auto start_index = vertex.id() * DIMENSION_COUNT;

    global_force(start_index + 0) += body_force.x() * multiplier;
    global_force(start_index + 1) += body_force.y() * multiplier;
    global_force(start_index + 2) += body_force.z() * multiplier;
}


template
class FEMFastForceAssembler<double>;
template
class FEMFastForceAssembler<float>;

}