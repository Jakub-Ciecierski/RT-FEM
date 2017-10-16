#include "RTFEM/FEM/Solver/FEMSolver.h"

#include <RTFEM/FEM/Solver/FEMAssembler.h>
#include <RTFEM/FEM/FEMModel.h>

#include <Eigen/LU>

#include <iostream>

namespace rtfem {

template<class T>
FEMSolver<T>::FEMSolver() :
    constitutive_solver_type_(ConstitutiveSolverType::LinearElastic),
    geometry_solver_type_(GeometrySolverType::Linear),
    analysis_solver_type_(AnalysisSolverType::Static) {}

template<class T>
FEMSolverOutput<T> FEMSolver<T>::Solve(const FEMModel<T>& fem_model) {
    // TODO pick solver based on this types
    FEMAssembler<T> fem_assembler;
    auto fem_assembler_data = fem_assembler.Compute(fem_model);

    FEMSolverOutput<T> output;
    output.displacement = SolveSystemOfEquations(fem_assembler_data);

    return output;
}

template<class T>
Eigen::Vector<T, Eigen::Dynamic> FEMSolver<T>::SolveSystemOfEquations(
    const FEMAssemblerData<T> &assembler_data) {

    return assembler_data.global_stiffness.fullPivLu().solve(
            assembler_data.global_force);
}

template
class FEMSolver<double>;
template
class FEMSolver<float>;

}