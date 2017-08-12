#include "RTFEM/FEM/Meshing/Tetrahedralization.h"

#include <RTFEM/FEM/Meshing/TriangleMesh.h>
#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/Meshing/AABBAlgorithms.h>

#include "tetgen.h"

namespace rtfem {

template<class T>
FEMGeometry<T> Tetrahedralization<T>::Compute(
        const TriangleMeshIndexed<T> &triangle_mesh){
    tetgenio tetgen_input, tetgen_output;

    constexpr unsigned int dimension = 3;
    tetgen_input.numberofpoints = triangle_mesh.points.size();
    tetgen_input.pointlist = new REAL[tetgen_input.numberofpoints * dimension];

    delete[] tetgen_input.pointlist;
}

template class Tetrahedralization<float>;
template class Tetrahedralization<double>;

}