#include "RTFEM/FEM/Meshing/AABBAlgorithms.h"

#include <RTFEM/FEM/Meshing/AABB.h>
#include <RTFEM/FEM/Meshing/TriangleMesh.h>

#include <algorithm>
#include <vector>

namespace rtfem {

template<class T>
AABB<T> CreateAABB(const TriangleMesh<T> &triangle_mesh) {
    Eigen::Vector3<T> min = Eigen::Vector3<T>::Zero();
    Eigen::Vector3<T> max = Eigen::Vector3<T>::Zero();

    if (triangle_mesh.triangles.size() == 0)
        return AABB<T>{min, max};

    min = triangle_mesh.triangles[0].v1;
    max = triangle_mesh.triangles[0].v1;

    for (const auto &triangle : triangle_mesh.triangles) {
        UpdateMin(triangle, min);
        UpdateMax(triangle, max);
    }

    return AABB<T>{min, max};
}

template<class T>
T FindMin(const TriangleFaceWithPoints<T> &triangle, const AABBCoordinate &coordinate) {
    auto i = GetIndex(coordinate);
    std::vector<T>
        min_vector = {triangle.v1[i], triangle.v2[i], triangle.v3[i]};

    auto min_iterator = std::min_element(min_vector.begin(), min_vector.end());
    auto min_index = std::distance(min_vector.begin(), min_iterator);

    return min_vector[min_index];
}

template<class T>
T FindMax(const TriangleFaceWithPoints<T> &triangle, const AABBCoordinate &coordinate) {
    auto i = GetIndex(coordinate);
    std::vector<T> vector = {triangle.v1[i], triangle.v2[i], triangle.v3[i]};

    auto iterator = std::max_element(vector.begin(), vector.end());
    auto index = std::distance(vector.begin(), iterator);

    return vector[index];
}

template<class T>
void UpdateMin(const TriangleFaceWithPoints<T> &triangle, Eigen::Vector3<T> &min) {
    auto min_x = FindMin<T>(triangle, AABBCoordinate::X);
    if (min[0] > min_x)
        min[0] = min_x;

    auto min_y = FindMin<T>(triangle, AABBCoordinate::Y);
    if (min[1] > min_y)
        min[1] = min_y;

    auto min_z = FindMin<T>(triangle, AABBCoordinate::Z);
    if (min[2] > min_z)
        min[2] = min_z;
}

template<class T>
void UpdateMax(const TriangleFaceWithPoints<T> &triangle, Eigen::Vector3<T> &max) {
    auto max_x = FindMax<T>(triangle, AABBCoordinate::X);
    if (max[0] < max_x)
        max[0] = max_x;

    auto max_y = FindMax<T>(triangle, AABBCoordinate::Y);
    if (max[1] < max_y)
        max[1] = max_y;

    auto max_z = FindMax<T>(triangle, AABBCoordinate::Z);
    if (max[2] < max_z)
        max[2] = max_z;
}

unsigned int GetIndex(const AABBCoordinate &aabb_coordinate) {
    switch (aabb_coordinate) {
        case AABBCoordinate::X:return 0;
        case AABBCoordinate::Y:return 1;
        case AABBCoordinate::Z:return 2;
    }
    return -1;
}

template AABB<float> CreateAABB<float>(const TriangleMesh<float> &triangle_mesh);
template AABB<double> CreateAABB<double>(const TriangleMesh<double> &triangle_mesh);

template float FindMin<float>(const TriangleFaceWithPoints<float> &triangle,
                              const AABBCoordinate &coordinate);
template double FindMin<double>(const TriangleFaceWithPoints<double> &triangle,
                                const AABBCoordinate &coordinate);

template float FindMax<float>(const TriangleFaceWithPoints<float> &triangle,
                              const AABBCoordinate &coordinate);
template double FindMax<double>(const TriangleFaceWithPoints<double> &triangle,
                                const AABBCoordinate &coordinate);

template void UpdateMin<float>(const TriangleFaceWithPoints<float> &triangle,
                               Eigen::Vector3<float> &min);
template void UpdateMin<double>(const TriangleFaceWithPoints<double> &triangle,
                                Eigen::Vector3<double> &min);

template void UpdateMax<float>(const TriangleFaceWithPoints<float> &triangle,
                               Eigen::Vector3<float> &max);
template void UpdateMax<double>(const TriangleFaceWithPoints<double> &triangle,
                                Eigen::Vector3<double> &max);
}