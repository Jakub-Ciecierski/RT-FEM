#include "RTFEM/FEM/Meshing/Intersections.h"

#include <RTFEM/FEM/Meshing/TriangleMesh.h>
#include <Eigen/Geometry>

namespace rtfem {

template<class T>
bool Contains(const Eigen::Vector3<T> point,
              const TriangleMesh<T>& triangle_mesh){
    Ray<T> ray{point, Eigen::Vector3<T>{1,0,0}};

    unsigned int intersection_count = 0;
    for(const auto& triangle : triangle_mesh.triangles){
        if(Intersects(ray, triangle))
            intersection_count++;
    }

    return ((intersection_count % 2) != 0);
}

template<class T>
bool Intersects(const Ray<T>& ray,
                const Triangle<T>& triangle){
    auto e1 = triangle.v2 - triangle.v1;
    auto e2 = triangle.v3 - triangle.v1;

    auto pvec = ray.direction.cross(e2);
    pvec.normalize();
    auto det = e1.dot(pvec);

    // Ray is parallel to plane
    if (det < 1e-8 && det > -1e-8) {
        return false;
    }

    auto inv_det = 1.0 / det;
    Eigen::Vector3<T> tvec = ray.origin - triangle.v1;
    tvec.normalize();
    auto u = tvec.dot(pvec) * inv_det;

    if (u < 0 || u > 1) {
        return false;
    }

    auto qvec = tvec.cross(e1);
    qvec.normalize();
    float v = ray.direction.dot(qvec) * inv_det;
    if (v < 0 || u + v > 1) {
        return false;
    }

    return true;
}

template struct Ray<float>;
template struct Ray<double>;

template bool Contains<float>(const Eigen::Vector3<float>,
                              const TriangleMesh<float>& triangle_mesh);
template bool Contains<double>(const Eigen::Vector3<double>,
                              const TriangleMesh<double>& triangle_mesh);

template bool Intersects<float>(const Ray<float> &ray,
                                const Triangle<float> &triangle);
template bool Intersects<double>(const Ray<double> &ray,
                                const Triangle<double> &triangle);

}