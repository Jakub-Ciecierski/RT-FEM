#include "RTFEM/FEM/Solver/FEMSolvers/FEMStaticSolver.h"

#include <RTFEM/FEM/Solver/FEMGlobalAssembler.h>

namespace rtfem {

template<class T>
FEMSolverOutput<T> FEMStaticSolver<T>::Solve(const FEMModel<T>& fem_model){
    FEMGlobalAssembler<T> fem_assembler;
    auto fem_assembler_data = fem_assembler.Compute(fem_model);

    FEMSolverOutput<T> output;

    output.displacement = this->SolveSystemOfEquations(
        fem_assembler_data.global_stiffness,
        fem_assembler_data.global_force);

    return output;
}

template
class FEMStaticSolver<double>;
template
class FEMStaticSolver<float>;

}