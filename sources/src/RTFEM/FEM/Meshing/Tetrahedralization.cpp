#include "RTFEM/FEM/Meshing/Tetrahedralization.h"

#include <RTFEM/DataStructure/Vector3.h>
#include <RTFEM/FEM/Meshing/TriangleMesh.h>
#include <RTFEM/FEM/FEMGeometry.h>

namespace rtfem {

FEMGeometry Tetrahedralization::Compute(const TriangleMesh &triangle_mesh,
                                        unsigned int vertex_count) {
    //auto

    return FEMGeometry();
}

std::vector<std::shared_ptr<Vector3>>
Tetrahedralization::GenerateRandomPointsInsideTriangleMesh(const TriangleMesh &triangle_mesh){

}

}