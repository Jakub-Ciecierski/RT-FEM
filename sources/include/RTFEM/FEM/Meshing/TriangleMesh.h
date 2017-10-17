#ifndef PROJECT_SURFACEMESH_H
#define PROJECT_SURFACEMESH_H

#include <RTFEM/DataTypes.h>

#include <vector>

namespace rtfem {

template<class T>
struct TriangleFaceWithPoints {
    Eigen::Vector3<T> v1;
    Eigen::Vector3<T> v2;
    Eigen::Vector3<T> v3;
};

template<class T>
struct TriangleMesh {
    std::vector<TriangleFaceWithPoints<T>> triangles;
};

struct TriangleFace {
    unsigned int v1;
    unsigned int v2;
    unsigned int v3;
};

/**
 * Represents a Triangle mesh.
 * The most common representation of surface meshes,
 * can be easily extracted from any renderable objects.
 */
template<class T>
struct TriangleMeshIndexed {
    std::vector<Eigen::Vector3<T>> points;
    std::vector<TriangleFace> triangles;
};

}

#endif //PROJECT_SURFACEMESH_H
