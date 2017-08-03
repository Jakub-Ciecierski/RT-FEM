#include "RTFEM/FEM/Meshing/Tetrahedralization.h"

#include <RTFEM/FEM/Meshing/TriangleMesh.h>
#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/Meshing/AABBAlgorithms.h>
#include <RTFEM/FEM/Meshing/AABB.h>
#include <RTFEM/FEM/Meshing/RandomDistributions.h>

namespace rtfem {

template<class T>
FEMGeometry<T> Tetrahedralization<T>::Compute(const TriangleMesh<T> &triangle_mesh,
                                              unsigned int vertex_count) {
   //auto

    return FEMGeometry<T>();
}

template<class T>
std::vector<Eigen::Vector3<T>>
Tetrahedralization<T>::GenerateRandomPointsInsideTriangleMesh(const TriangleMesh<T> &triangle_mesh,
                                                              unsigned int vertex_count){
    std::vector<Eigen::Vector3<T>> vertices(vertex_count);
    auto aabb = CreateAABB<T>(triangle_mesh);

    auto random_distributions = CreateRandomDistributions(aabb);
    for(unsigned int i = 0; i < vertex_count; i++){
        auto vertex = GenereRandomValidPoint(triangle_mesh, random_distributions);
    }

}

template<class T>
Eigen::Vector3<T>
Tetrahedralization<T>::GenereRandomValidPoint(const TriangleMesh<T> triangle_mesh,
                                              RandomDistributions<T>& random_distributions) {
    Eigen::Vector3<T> random_point{
            random_distributions.GenerateX(),
            random_distributions.GenerateY(),
            random_distributions.GenerateZ()
    };

    //IsInside(triangle_mesh, random_point);
}

template<class T>
RandomDistributions<T> Tetrahedralization<T>::CreateRandomDistributions(const AABB<T>& aabb){
    RandomDistributions<T> random_distributions(aabb);

    return random_distributions;
}

template class Tetrahedralization<float>;
template class Tetrahedralization<double>;

}