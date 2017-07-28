#ifndef PROJECT_SURFACEMESH_H
#define PROJECT_SURFACEMESH_H

#include <RTFEM/DataStructure/Vector3.h>
#include <vector>

namespace rtfem {

struct Triangle{
    rtfem::Vector3 v1;
    rtfem::Vector3 v2;
    rtfem::Vector3 v3;
};

/**
 * Represents a Triangle mesh.
 * The most common representation of surface meshes,
 * can be easily extracted from any renderable objects.
 */
struct TriangleMesh {
    std::vector<Triangle> triangles;
};

}


#endif //PROJECT_SURFACEMESH_H
