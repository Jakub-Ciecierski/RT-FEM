#ifndef PROJECT_AABBALGORITHMS_H
#define PROJECT_AABBALGORITHMS_H

#include <RTFEM/DataTypes.h>

namespace rtfem {

template<class T>
struct AABB;

template<class T>
struct TriangleMesh;

template<class T>
struct TriangleFaceWithPoints;

enum class AABBCoordinate {
    X, Y, Z
};

template<class T>
AABB<T> CreateAABB(const TriangleMesh<T> &triangle_mesh);

template<class T>
T FindMin(const TriangleFaceWithPoints<T> &triangle, const AABBCoordinate &coordinate);

template<class T>
T FindMax(const TriangleFaceWithPoints<T> &triangle, const AABBCoordinate &coordinate);

template<class T>
void UpdateMin(const TriangleFaceWithPoints<T> &triangle, Eigen::Vector3<T> &min);

template<class T>
void UpdateMax(const TriangleFaceWithPoints<T> &triangle, Eigen::Vector3<T> &max);

unsigned int GetIndex(const AABBCoordinate &aabb_coordinate);

}

#endif //PROJECT_AABBALGORITHMS_H
