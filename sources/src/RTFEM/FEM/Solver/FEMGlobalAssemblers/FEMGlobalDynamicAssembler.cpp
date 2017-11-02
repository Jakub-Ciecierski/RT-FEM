#include <RTFEM/FEM/Material.h>
#include <RTFEM/FEM/Solver/FiniteElementSolver.h>
#include "RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMGlobalDynamicAssembler.h"

namespace rtfem {

template<class T>
void FEMGlobalDynamicAssembler<T>::ComputeAssemblerDataIteration(
        FEMGlobalAssemblerData<T> &fem_assembler_data,
        const FiniteElementSolverData<T>& finite_element_solver_data,
        const FEMModel<T> &fem_model,
        const Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N> &
        constitutive_matrix_C,
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
        boolean_assembly_matrix_A) {
    FEMGlobalAssembler<T>::ComputeAssemblerDataIteration(fem_assembler_data,
                                                         finite_element_solver_data,
                                                         fem_model,
                                                         constitutive_matrix_C,
                                                         boolean_assembly_matrix_A);
    fem_assembler_data.global_mass +=
            ComputePartialGlobalMassMatrix(fem_model.material().density,
                                           finite_element_solver_data.volume,
                                           boolean_assembly_matrix_A);
}

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
FEMGlobalDynamicAssembler<T>::ComputePartialGlobalMassMatrix(
        T density, T volume,
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
        boolean_assembly_matrix_A){
    auto local_mass_matrix = ComputeLocalMassMatrix(density, volume);
};

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
FEMGlobalDynamicAssembler<T>::ComputeLocalMassMatrix(T density, T volume){

};

template<class T>
void FEMGlobalDynamicAssembler<T>::ApplyBoundaryConditions(
        FEMGlobalAssemblerData<T> &assembler_data,
        const BoundaryConditionContainer<T> &boundary_conditions){
    FEMGlobalAssembler<T>::ApplyBoundaryConditions(assembler_data,
                                                   boundary_conditions);
}



}