#ifndef PROJECT_FEMDYNAMICASSEMBLER_H
#define PROJECT_FEMDYNAMICASSEMBLER_H

#include <RTFEM/FEM/Solver/FEMGlobalAssembler.h>

namespace rtfem {

template<class T>
class FEMGlobalDynamicAssembler : public FEMGlobalAssembler<T>{
public:
    FEMGlobalDynamicAssembler() = default;
    ~FEMGlobalDynamicAssembler() = default;

protected:
    virtual void ComputeAssemblerData(
        FEMGlobalAssemblerData<T> &fem_assembler_data,
        const FEMModel<T> &fem_model,
        Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N> &
        constitutive_matrix_C) override;

    virtual void ComputeAssemblerDataIteration(
            FEMGlobalAssemblerData<T> &fem_assembler_data,
            const FiniteElementSolverData<T>& finite_element_solver_data,
            const FEMModel<T> &fem_model,
            const Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N>
            &constitutive_matrix_C,
            const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
            boolean_assembly_matrix_A) override;

    virtual void ApplyBoundaryConditionsToFEM(
        FEMGlobalAssemblerData<T> &assembler_data,
        const BoundaryConditionContainer<T> &boundary_conditions) override;
private:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    ComputePartialGlobalMassMatrix(
            T density, T volume,
            const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
            boolean_assembly_matrix_A);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    ComputeLocalMassMatrix(T density, T volume);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    ComputeGlobalDampingMatrix(
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
        global_mass_matrix,
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
        global_stiffness_matrix,
        T damping_mass,
        T damping_stiffness);

};
}


#endif //PROJECT_FEMDYNAMICASSEMBLER_H

