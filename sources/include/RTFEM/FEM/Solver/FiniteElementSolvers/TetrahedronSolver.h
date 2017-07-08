#ifndef PROJECT_TETRAHEDRONSOLVER_H
#define PROJECT_TETRAHEDRONSOLVER_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/Matrix.h>
#include <memory>
#include <RTFEM/DataStructure/Vector3.h>
#include <RTFEM/FEM/Solver/FiniteElementSolver.h>

namespace rtfem {

class TetrahedronFiniteElement;

/**
 * Solver for Linear Tetrahedron (constant gradient of shape function).
 * The geometry matrix B is constant with respect to X
 * (and always will be constant, independent of material/strain equations
 *
 * Solver Data:
 *  Shape Matrix: [3 x 12]
 *  Used in computing Force vector
 *
 *  Geometry matrix: [6 x 12] (Shape function gradient)
 *  Used in computing Stiffness Matrix
 */
class TetrahedronSolver : public FiniteElementSolver{
public:
    TetrahedronSolver();

    ~TetrahedronSolver();

    virtual FiniteElementSolverData Solve(std::shared_ptr<FiniteElement> finite_element) override;

private:
    /**
     * Compute A/B/C/D coefficients of Shape Function.
     *
     * Determinant of 3x3 matrix made of column wise vertex coordinates.
     *
     * @param vertex_local_index
     * @param finite_element
     * @return
     */
    Float ComputeA(UInt vertex_local_index,
                   const std::shared_ptr<FiniteElement> finite_element);
    Float ComputeB(UInt vertex_local_index,
                   const std::shared_ptr<FiniteElement> finite_element);
    Float ComputeC(UInt vertex_local_index,
                   const std::shared_ptr<FiniteElement> finite_element);
    Float ComputeD(UInt vertex_local_index,
                   const std::shared_ptr<FiniteElement> finite_element);

    Float ComputeVolume(const std::shared_ptr<FiniteElement> finite_element);

    Float ComputeShapeFunctionValue(Float a, Float b,
                                    Float c, Float d,
                                    Float volume,
                                    const Vector3& x);

    /**
     * Assembles the constant [6x12] Geometry Matrix (B).
     * @param B
     * @param column_offset
     * @param b
     * @param c
     * @param d
     */
    void AssemblyGeometryMatrix(Matrix& B, UInt column_offset,
                                Float b, Float c, Float d);
};
}


#endif //PROJECT_TETRAHEDRONSOLVER_H
