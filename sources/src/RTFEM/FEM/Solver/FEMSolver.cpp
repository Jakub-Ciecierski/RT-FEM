#include "RTFEM/FEM/Solver/FEMSolver.h"

#include <RTFEM/DataStructure/Matrix.h>
#include <RTFEM/Memory/UniquePointer.h>

namespace rtfem {

FEMSolver::FEMSolver(const ConstitutiveSolverType &&constitutive_solver_type,
                     const GeometrySolverType &&geometry_solver_type,
                     const AnalysisSolverType &&analysis_solver_type) :
        constitutive_solver_type_(constitutive_solver_type),
        geometry_solver_type_(geometry_solver_type),
        analysis_solver_type_(analysis_solver_type) {}

FEMSolver::~FEMSolver() {}

void FEMSolver::Solve(const std::shared_ptr<FEMModel> fem_model) {
    // TODO pick solver based on this types

    //auto GlobalStiffness = ComputeGlobalStiffness();
    //auto GlobalForces = ComputeGlobalForces();
}

}