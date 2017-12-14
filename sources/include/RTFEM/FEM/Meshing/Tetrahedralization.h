#ifndef PROJECT_TRIANGULATOR_H
#define PROJECT_TRIANGULATOR_H

#include <RTFEM/DataTypes.h>

#include <vector>
#include <memory>

#include <random>

class tetgenio;
class tetgenbehavior;

namespace rtfem {

template<class T>
struct FEMGeometry;

template<class T>
struct TriangleMeshIndexed;

struct TetrahedralizationOptions {
    /**
     * No tetrahedra will be generated with volume greater than maximum_volume.
     * 0 indicates no constraint
     */
    float maximum_volume = 0;

    float max_radius_edge_ratio = 2.0f;
    int min_dihedral_angle_degree = 0;
};

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

    void SetOptions(const TetrahedralizationOptions &options);

    /**
     * Computes the Tetrahedralization.
     *
     * @param triangle_mesh
     * @return
     */
    FEMGeometry<T> Compute(
        const TriangleMeshIndexed<T> &triangle_mesh);

private:
    void SetupInput(const TriangleMeshIndexed<T> &triangle_mesh,
                    tetgenio &tetgen_input,
                    tetgenbehavior &tetgen_options);
    void SetupInputPoints(const TriangleMeshIndexed<T> &triangle_mesh,
                          tetgenio &tetgen_input);
    void SetupInputFacets(const TriangleMeshIndexed<T> &triangle_mesh,
                          tetgenio &tetgen_input);
    void SetupInputOptions(tetgenbehavior &tetgen_options);

    FEMGeometry<T> FetchOutput(tetgenio &tetgen_output);
    void FetchPoints(FEMGeometry<T> &fem_geometry,
                     tetgenio &tetgen_output);
    void FetchTetrahedra(FEMGeometry<T> &fem_geometry,
                         tetgenio &tetgen_output);
    std::vector<unsigned int> FetchTetrahedronFaces(
        FEMGeometry<T> &fem_geometry,
        unsigned int v1,
        unsigned int v2,
        unsigned int v3,
        unsigned int v4);
    void FetchBoundaryFaces(FEMGeometry<T> &fem_geometry,
                            tetgenio &tetgen_output);

    TetrahedralizationOptions options_;
};
}

#endif //PROJECT_TRIANGULATOR_H
