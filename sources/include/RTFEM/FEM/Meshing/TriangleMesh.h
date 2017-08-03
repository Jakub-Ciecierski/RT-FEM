#ifndef PROJECT_SURFACEMESH_H
#define PROJECT_SURFACEMESH_H

#include <RTFEM/DataTypes.h>

#include <vector>

namespace rtfem {

template<class T>
struct Triangle {
    Eigen::Vector3<T> v1;
    Eigen::Vector3<T> v2;
    Eigen::Vector3<T> v3;
};

/**
 * Represents a Triangle mesh.
 * The most common representation of surface meshes,
 * can be easily extracted from any renderable objects.
 */
template<class T>
struct TriangleMesh {
    std::vector<Triangle<T>> triangles;
};

}


#endif //PROJECT_SURFACEMESH_H
