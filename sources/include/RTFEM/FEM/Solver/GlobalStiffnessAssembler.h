#ifndef PROJECT_GLOBALSTIFFNESSASSEMBLER_H
#define PROJECT_GLOBALSTIFFNESSASSEMBLER_H

#include <memory>

#include <RTFEM/DataStructure/Matrix.h>

namespace rtfem {

class FEMModel;
class FiniteElement;

class FiniteElementSolver;
enum class FiniteElementType;

/**
 * Computes and Assembles Global Stiffness Matrix.
 *
 * Local Stiffness Matrix (k) is the stiffness of each element [3Ne x 3Ne]
 * e.g. for Tetrahedron (Ne = 4) thus: (dim = [12 x 12])
 *
 * Global Stiffness Matrix (K) is the stiffness of entire FEM Model [3N x 3N]
 * e.g. For 9 vertices (dim = [27 x 27])
 *
 * Partial Global Stiffness Matrix (Ke) is the matrix of dimension equal to Global Stiffness
 * but filled only with Local Stiffness data.
 */
class GlobalStiffnessAssembler {
public:
    GlobalStiffnessAssembler();

    ~GlobalStiffnessAssembler();

    /**
     * Computes Global Stiffness Matrix (K).
     *
     *      1) Computes Constitutive Matrix (C)
     *      2) Computes Geometry Matrix (B) for each Finite Element
     *      3) Computes Local Stiffness (k) for each Finite Element
     *          - Using Constitutive Matrix and Geometry Matrix.
     *      4) Assembles all Local Stiffness matrices into Global Stiffness Matrix (K)
     *
     * @param fem_model
     * @return
     */
    Matrix Compute(const std::shared_ptr<FEMModel> fem_model);

private:
    /**
     * Computes Isotropic Constitutive Matrix (C).
     * Used to compute Local Stiffness.
     *
     * In case of Linear Elasticity, Constitutive Matrix is constant.
     *
     * Isotropy allows Constitutive Matrix to remain unchanged
     * after orthogonal transformations
     *
     * @param fem_model
     * @return
     */
    Matrix ComputeConstitutiveMatrix(const std::shared_ptr<FEMModel> fem_model);

    /**
     * Computes Local Stiffness (k) of a given element.
     *
     * @param finite_element
     * @param constitutive_matrix
     * @return
     */
    Matrix
    ComputeLocalStiffness(const std::shared_ptr<FiniteElement> finite_element,
                          const Matrix &constitutive_matrix);

    /**
     * Computes Geometry Matrix (B).
     * Used to compute Local Stiffness
     *
     * @param type
     * @return
     */
    std::unique_ptr<FiniteElementSolver> GetFiniteElementSolver(const FiniteElementType& type);

    /**
     * Computes The Boolean Assembly Matrix (A) [3Ne x 3N].
     *
     * Maps Local Stiffness to Global Stiffness.
     *
     * @param finite_element
     * @param vertex_count
     * @return
     */
    Matrix ComputeBooleanAssemblyMatrix(const std::shared_ptr<FiniteElement> finite_element,
                                        UInt vertex_count);

    /**
     * Computes Partial Global Stiffness Matrix.
     * Uses Boolean Assembly Matrix to map a single Local Stiffness to Global Stiffness.
     *
     * @param boolean_assembly_matrix_A
     * @param local_stiffness_k
     * @return
     */
    Matrix ComputePartialGlobalStiffnessMatrix(const Matrix &boolean_assembly_matrix_A,
                                               const Matrix &local_stiffness_k);

};
}


#endif //PROJECT_GLOBALSTIFFNESSASSEMBLER_H
