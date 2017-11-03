#include "RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMGlobalDynamicAssembler.h"

#include <RTFEM/FEM/Material.h>
#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/Solver/FiniteElementSolver.h>

namespace rtfem {

template<class T>
void FEMGlobalDynamicAssembler<T>::ComputeAssemblerData(
    FEMGlobalAssemblerData<T> &fem_assembler_data,
    const FEMModel<T> &fem_model,
    Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N> &
    constitutive_matrix_C){
    FEMGlobalAssembler<T>::ComputeAssemblerData(fem_assembler_data,
                                                fem_model,
                                                constitutive_matrix_C);

    fem_assembler_data.global_damping =
        ComputeGlobalDampingMatrix(fem_assembler_data.global_mass,
                                   fem_assembler_data.global_stiffness,
                                   fem_model.material().damping_mass,
                                   fem_model.material().damping_stiffness);
}

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

    return boolean_assembly_matrix_A.transpose() * local_mass_matrix
        * boolean_assembly_matrix_A;
};

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
FEMGlobalDynamicAssembler<T>::ComputeLocalMassMatrix(T density, T volume){
    Eigen::Matrix<T, 12, 12> mass_matrix;
    mass_matrix <<
                2,0,0,1,0,0,1,0,0,1,0,0,
                0,2,0,0,1,0,0,1,0,0,1,0,
                0,0,2,0,0,1,0,0,1,0,0,1,
                1,0,0,2,0,0,1,0,0,1,0,0,
                0,1,0,0,2,0,0,1,0,0,1,0,
                0,0,1,0,0,2,0,0,1,0,0,1,
                1,0,0,1,0,0,2,0,0,1,0,0,
                0,1,0,0,1,0,0,2,0,0,1,0,
                0,0,1,0,0,1,0,0,2,0,0,1,
                1,0,0,1,0,0,1,0,0,2,0,0,
                0,1,0,0,1,0,0,1,0,0,2,0,
                0,0,1,0,0,1,0,0,1,0,0,2;
    const T divider = 20;
    mass_matrix *= density * volume / divider;

    return mass_matrix;
};

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
FEMGlobalDynamicAssembler<T>::ComputeGlobalDampingMatrix(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
    global_mass_matrix,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
    global_stiffness_matrix,
    T damping_mass,
    T damping_stiffness){
    return
        global_mass_matrix * damping_mass +
        global_stiffness_matrix * damping_stiffness;
};

template
class FEMGlobalDynamicAssembler<double>;
template
class FEMGlobalDynamicAssembler<float>;

}