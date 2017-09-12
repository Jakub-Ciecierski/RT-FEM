#ifndef PROJECT_INTERSECTIONS_H
#define PROJECT_INTERSECTIONS_H

#include <RTFEM/DataTypes.h>

namespace rtfem {

template<class T>
struct Triangle;

template<class T>
struct TriangleMesh;

template<class T>
struct Ray {
    Eigen::Vector3<T> origin;
    Eigen::Vector3<T> direction;
};

/**
 * Checks whether point is inside closed triangle_mesh
 * @tparam T
 * @param triangle_mesh
 * @return
 */
template<class T>
bool Contains(const Eigen::Vector3<T>,
              const TriangleMesh<T> &triangle_mesh);

/**
 * https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
 *
 * @tparam T
 * @param ray
 * @param triangle
 * @return
 */
template<class T>
bool Intersects(const Ray<T> &ray,
                const Triangle<T> &triangle);

}

#endif //PROJECT_INTERSECTIONS_H
