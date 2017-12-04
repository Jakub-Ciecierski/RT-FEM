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

template<class T>
struct TriangleFace {
    TriangleFace(
        unsigned int v1_,
        unsigned int v2_,
        unsigned int v3_) : v1(v1_), v2(v2_), v3(v3_),
                            is_boundary_face(false),
                            traction_force(0),
                            area(0),
                            normal(Eigen::Vector3<T>::Zero()),
                            B(0), C(0), D(0) {}

    TriangleFace(
        unsigned int v1_,
        unsigned int v2_,
        unsigned int v3_,
        bool is_boundary_face_) : v1(v1_), v2(v2_), v3(v3_),
                                  is_boundary_face(is_boundary_face_),
                                  traction_force(0),
                                  area(0),
                                  normal(Eigen::Vector3<T>::Zero()),
                                  B(0), C(0), D(0){}
    unsigned int v1;
    unsigned int v2;
    unsigned int v3;

    bool is_boundary_face;

    T traction_force;

    T area;
    Eigen::Vector3<T> normal;

    // Used to assemble traction force
    T B;
    T C;
    T D;

    bool operator==(const TriangleFace& other) const {
        std::vector<unsigned int> indices{v1, v2, v3};
        std::vector<unsigned int> other_indices{other.v1, other.v2, other.v3};

        for(unsigned int i = 0; i < 3; i++){
            bool valid_i_candidate = false;
            for(unsigned int j = 0; j < 3; j++){
                if(indices[i] == other_indices[j]){
                    valid_i_candidate = true;
                }
            }
            if(!valid_i_candidate)
                return false;
        }
        return true;
    }
};

/**
 * Represents a Triangle mesh.
 * The most common representation of surface meshes,
 * can be easily extracted from any renderable objects.
 */
template<class T>
struct TriangleMeshIndexed {
    std::vector<Eigen::Vector3<T>> points;
    std::vector<TriangleFace<T>> triangles;
};

}

#endif //PROJECT_SURFACEMESH_H
