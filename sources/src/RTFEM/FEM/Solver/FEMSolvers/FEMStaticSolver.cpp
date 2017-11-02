#include "RTFEM/FEM/Solver/FEMSolvers/FEMStaticSolver.h"

#include <RTFEM/FEM/Solver/FEMGlobalAssembler.h>

#include <Eigen/LU>

namespace rtfem {

template<class T>
FEMSolverOutput<T> FEMStaticSolver<T>::Solve(const FEMModel<T>& fem_model){
    FEMGlobalAssembler<T> fem_assembler;
    auto fem_assembler_data = fem_assembler.Compute(fem_model);

    FEMSolverOutput<T> output;
    output.displacement = SolveSystemOfEquations(fem_assembler_data);

    return output;
}

template<class T>
Eigen::Vector<T, Eigen::Dynamic> FEMStaticSolver<T>::SolveSystemOfEquations(
        const FEMGlobalAssemblerData<T> &assembler_data) {
    return assembler_data.global_stiffness.fullPivLu().solve(
            assembler_data.global_force);
}

template
class FEMStaticSolver<double>;
template
class FEMStaticSolver<float>;

}