#include <RTFEM/FEM/Material.h>
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

    for(unsigned int i = 0; i < global_force.size(); i++){
        global_force(i) = 0;
    }

    for(unsigned int i = 0; i < fem_geometry.finite_elements.size(); i++){
        auto& finite_element = fem_geometry.finite_elements[i];
        const auto& vertices_indices = finite_element->vertices_indices();

        auto multiplier = material.density * (finite_element->volume() / 4.0);

        AddBodyForce(*fem_geometry.vertices[vertices_indices[0]],
                     body_force, multiplier, global_force);
        AddBodyForce(*fem_geometry.vertices[vertices_indices[1]],
                     body_force, multiplier, global_force);
        AddBodyForce(*fem_geometry.vertices[vertices_indices[2]],
                     body_force, multiplier, global_force);
        AddBodyForce(*fem_geometry.vertices[vertices_indices[3]],
                     body_force, multiplier, global_force);
    }
}

template<class T>
void FEMFastForceAssembler<T>::AddBodyForce(
        const Vertex<T>& vertex,
        const Eigen::Vector3<T>& body_force,
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