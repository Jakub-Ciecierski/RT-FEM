#include "RTFEM/FEM/Solver/FEMSolvers/FEMStaticSolver.h"

#include <RTFEM/FEM/Solver/FEMGlobalAssembler.h>

namespace rtfem {

template<class T>
FEMStaticSolver<T>::FEMStaticSolver(FEMModel<T> *fem_model) :
    FEMSolver<T>(fem_model) {}

template<class T>
FEMSolverOutput<T> FEMStaticSolver<T>::Solve(){
    FEMGlobalAssembler<T> fem_assembler;
    auto fem_assembler_data = fem_assembler.Compute(*this->fem_model_);

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