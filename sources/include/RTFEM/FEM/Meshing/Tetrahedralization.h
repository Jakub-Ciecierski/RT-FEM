#ifndef PROJECT_TRIANGULATOR_H
#define PROJECT_TRIANGULATOR_H

#include <RTFEM/DataTypes.h>

#include <vector>
#include <memory>

#include <random>

namespace rtfem {

template<class T>
struct FEMGeometry;

template<class T>
struct TriangleMeshIndexed;


/**
 * Generates 3D Tetrahedron Mesh.
 *
 * Tetrahedralization is an offline tool used as a pre-processor to FEM Solver.
 */
template<class T>
class Tetrahedralization {
public:
    Tetrahedralization() = default;
    ~Tetrahedralization() = default;

    /**
     * Computes the Tetrahedralization.
     *
     * @param triangle_mesh
     * @return
     */
    FEMGeometry<T> Compute(const TriangleMeshIndexed<T> &triangle_mesh);

private:

};
}


#endif //PROJECT_TRIANGULATOR_H
