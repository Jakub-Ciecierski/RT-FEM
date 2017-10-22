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

    bool operator==(const TriangleFace& other){
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
    std::vector<TriangleFace> triangles;
};

}

#endif //PROJECT_SURFACEMESH_H
