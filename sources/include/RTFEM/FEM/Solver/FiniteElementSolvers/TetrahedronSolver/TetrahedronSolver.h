#ifndef PROJECT_TETRAHEDRONSOLVER_H
#define PROJECT_TETRAHEDRONSOLVER_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/Matrix.h>
#include <memory>
#include <RTFEM/DataStructure/Vector3.h>
#include <RTFEM/FEM/Solver/FiniteElementSolver.h>
#include "RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver/TetrahedronSolverData.h"

namespace rtfem {

class Vertex;
class TetrahedronFiniteElement;

/**
 * Solver for Linear Tetrahedron (constant gradient of shape function).
 * The geometry matrix B is constant with respect to X
 * (and always will be constant, independent of material/strain equations)
 * .
 * The simplicity of Linear Tetrahedron allows for analytical integration.
 *
 * Solver Data:
 *
 *  Geometry matrix: [6 x 12] (Shape function gradient)
 *      Used in computing Stiffness Matrix
 *
 *  Force vector: [12x1]
 *      Body force (e.g. gravity) (f)
 *      Traction force (e.g. collision) (t)
 *
 * NOTE:
 *  1) Special care must be taken for Local vertex ordering (TetrahedronSolver does NOT do any ordering)
 *      Swap local ordering of vertices
 *      Pick a vertex as initial one 1.
 *      Pick a face containing the first three vertices 1,2,3 (including the initial one).
 *      The vertex 4 is the one opposite of chosen face.
 *      Number the first three vertices in counterclockwise wise order as seen from vertex 4
 *
 *  2) Research Materials:
 *      http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
 */
class TetrahedronSolver : public FiniteElementSolver{
public:
    TetrahedronSolver() = default;
    ~TetrahedronSolver() = default;

    virtual FiniteElementSolverData Solve(std::shared_ptr<FiniteElement> finite_element) override;

    virtual FiniteElementSolverData Solve(std::shared_ptr<FiniteElement> finite_element,
                                          const Eigen::Vector3<Float>& body_force,
                                          const TractionForce& traction_force) override;

    Eigen::Matrix<Float, TETRAHEDRON_JACOBIAN_MATRIX_N, TETRAHEDRON_JACOBIAN_MATRIX_N>
    SolveJacobianInverse(std::shared_ptr<FiniteElement> finite_element);

private:
    Edges ComputeEdgesCache(const Vertex &v1, const Vertex &v2,
                            const Vertex &v3, const Vertex &v4);

    FacesArea ComputeFacesArea(const Edges& edges);

    /**
     * Computes The coefficients on linear tetrahedron shape function.
     * fi(x) = (Ai + Bi*x1 + Ci*x2 + Di*x3) / 6*V
     *
     * @param v1
     * @param v2
     * @param v3
     * @param v4
     * @return
     */
    TetrahedronShapeFunctionCoefficients
    ComputeShapeFunctionCoefficients(const Vertex& v1, const Vertex& v2,
                                     const Vertex& v3, const Vertex& v4,
                                     const Edges& edges);
    Eigen::Vector4<Float> ComputeAi(const Vertex &v1, const Vertex &v2,
                                    const Vertex &v3, const Vertex &v4);

    Float ComputeVolume(const TetrahedronShapeFunctionCoefficients& coefficients);

    /**
     * Computes The [6x12] Geometric Matrix (B)
     *
     * @param coefficients
     * @param volume
     * @return
     */
    Eigen::Matrix<Float, TETRAHEDRON_GEOMETRIC_MATRIX_N, TETRAHEDRON_GEOMETRIC_MATRIX_M>
    ComputeGeometryMatrix(const TetrahedronShapeFunctionCoefficients& coefficients, Float volume);
    void AssemblyGeometryMatrix(
            Eigen::Matrix<Float, TETRAHEDRON_GEOMETRIC_MATRIX_N, TETRAHEDRON_GEOMETRIC_MATRIX_M> &B,
            UInt column_offset,
            Float b, Float c, Float d);

    /**
     * Computes Force vector [12x1] which is a combination of body force and traction force.
     *
     * @param shape_function_coefficients
     * @param body_force
     * @param traction_force
     * @return
     */
    Eigen::Vector<Float, TETRAHEDRON_FORCE_VECTOR_N>
    ComputeForceVector(const TetrahedronShapeFunctionCoefficients &shape_function_coefficients,
                       Float volume, const FacesArea &faces_area,
                       const Eigen::Vector3<Float> &body_force,
                       const TractionForce &traction_force);
    Eigen::Vector<Float, TETRAHEDRON_FORCE_VECTOR_N>
    ComputeBodyForceVector(Float volume,
                           const Eigen::Vector3<Float> &body_force);
    Eigen::Vector<Float, TETRAHEDRON_FORCE_VECTOR_N>
    ComputeTractionForceVector(const TetrahedronShapeFunctionCoefficients &shape_function_coefficients,
                               const FacesArea &faces_area,
                               const TractionForce &traction_force);
    Eigen::Vector<Float, TETRAHEDRON_FORCE_VECTOR_N>
    ComputeTractionForceVectorFace(UInt face_index, Float traction_force,
                                   Float area,
                                   Float B, Float C, Float D);

    Eigen::Matrix<Float, TETRAHEDRON_JACOBIAN_MATRIX_N, TETRAHEDRON_JACOBIAN_MATRIX_N>
    AssemblyJacobianInverse(const TetrahedronShapeFunctionCoefficients &coefficients,
                            Float volume);
};
}


#endif //PROJECT_TETRAHEDRONSOLVER_H
