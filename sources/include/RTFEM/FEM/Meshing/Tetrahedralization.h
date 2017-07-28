#ifndef PROJECT_TRIANGULATOR_H
#define PROJECT_TRIANGULATOR_H

#include <vector>
#include <memory>

namespace rtfem {

struct FEMGeometry;
struct TriangleMesh;
struct Vector3;

/**
 * Generates 3D Tetrahedron Mesh.
 *
 * Tetrahedralization is an offline tool used as a pre-processor to FEM Solver.
 */
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
    FEMGeometry Compute(const TriangleMesh& triangle_mesh,
                        unsigned int vertex_count);
private:
    std::vector<std::shared_ptr<Vector3>>
    GenerateRandomPointsInsideTriangleMesh(const TriangleMesh &triangle_mesh);
};
}


#endif //PROJECT_TRIANGULATOR_H
