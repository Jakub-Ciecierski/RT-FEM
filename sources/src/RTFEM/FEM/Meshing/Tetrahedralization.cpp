#include "RTFEM/FEM/Meshing/Tetrahedralization.h"

#include <RTFEM/FEM/Meshing/TriangleMesh.h>
#include <RTFEM/FEM/FEMGeometry.h>

namespace rtfem {

template<class T>
FEMGeometry<T> Tetrahedralization<T>::Compute(const TriangleMesh &triangle_mesh,
                                              unsigned int vertex_count) {
   //auto

    return FEMGeometry<T>();
}

template<class T>
std::vector<std::shared_ptr<Eigen::Vector3<T>>>
Tetrahedralization<T>::GenerateRandomPointsInsideTriangleMesh(const TriangleMesh &triangle_mesh){

}

}