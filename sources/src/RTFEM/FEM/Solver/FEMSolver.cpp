#include "RTFEM/FEM/Solver/FEMSolver.h"

#include <RTFEM/FEM/Solver/FEMAssembler.h>

namespace rtfem {

template<class T>
FEMSolver<T>::FEMSolver(const ConstitutiveSolverType &&constitutive_solver_type,
                        const GeometrySolverType &&geometry_solver_type,
                        const AnalysisSolverType &&analysis_solver_type) :
    constitutive_solver_type_(constitutive_solver_type),
    geometry_solver_type_(geometry_solver_type),
    analysis_solver_type_(analysis_solver_type) {}

template<class T>
void FEMSolver<T>::Solve(const std::shared_ptr<FEMModel<T>> fem_model) {
    // TODO pick solver based on this types
    FEMAssembler<T> fem_assembler;
    auto fem_assembler_data = fem_assembler.Compute(fem_model);

    displacements_ = SolveSystemOfEquations(fem_assembler_data);
}

template<class T>
Eigen::Vector<T, Eigen::Dynamic> FEMSolver<T>::SolveSystemOfEquations(
    const FEMAssemblerData<T> &assembler_data) {
    /*
    return assembler_data.global_stiffness.lu().solve(
            assembler_data.global_force);*/
}

template
class FEMSolver<double>;
template
class FEMSolver<float>;

}