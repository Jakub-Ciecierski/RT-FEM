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
struct TriangleMesh;

template<class T>
struct AABB;

template<class T>
struct RandomDistributions;

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
     * @param vertex_count
     * Number of vertices in the FEMGeometry
     * @return
     */
    FEMGeometry<T> Compute(const TriangleMesh<T> &triangle_mesh,
                           unsigned int vertex_count);
private:
    std::vector<Eigen::Vector3<T>>
    GenerateRandomPointsInsideTriangleMesh(const TriangleMesh<T> &triangle_mesh,
                                           unsigned int vertex_count);

    Eigen::Vector3<T> GenereRandomValidPoint(const TriangleMesh<T> triangle_mesh,
                                             RandomDistributions<T>& random_distributions);

    RandomDistributions<T> CreateRandomDistributions(const AABB<T>& aabb);
};
}


#endif //PROJECT_TRIANGULATOR_H
