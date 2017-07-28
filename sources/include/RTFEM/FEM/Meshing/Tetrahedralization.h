#ifndef PROJECT_TRIANGULATOR_H
#define PROJECT_TRIANGULATOR_H

#include <RTFEM/DataTypes.h>

#include <vector>
#include <memory>

namespace rtfem {

template<class T>
struct FEMGeometry;

struct TriangleMesh;

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
     * @param triangle_mesh
     * @param vertex_count
     * Number of vertices in the FEMGeometry
     * @return
     */
    FEMGeometry<T> Compute(const TriangleMesh &triangle_mesh,
                           unsigned int vertex_count);
private:
    std::vector<std::shared_ptr<Eigen::Vector3<T>>>
    GenerateRandomPointsInsideTriangleMesh(const TriangleMesh &triangle_mesh);
};
}


#endif //PROJECT_TRIANGULATOR_H
