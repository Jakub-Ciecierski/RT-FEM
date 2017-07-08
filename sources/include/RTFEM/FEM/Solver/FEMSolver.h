#ifndef PROJECT_FEMSOLVER_H
#define PROJECT_FEMSOLVER_H

#include <RTFEM/FEM/Solver/FEMSolverTypes.h>

#include <memory>

namespace rtfem {

class FEMModel;
class Matrix;

class FEMSolver {
public:
    FEMSolver(
            const ConstitutiveSolverType&& constitutive_solver_type,
            const GeometrySolverType&& geometry_solver_type,
            const AnalysisSolverType&& analysis_solver_type);
    ~FEMSolver();

    const ConstitutiveSolverType& constitutive_solver_type(){
        return constitutive_solver_type_;
    }
    const GeometrySolverType& geometry_solver_type(){
        return geometry_solver_type_;
    }
    const AnalysisSolverType& analysis_solver_type(){
        return analysis_solver_type_;
    }

    void Solve(const std::shared_ptr<FEMModel> fem_model);

private:
    ConstitutiveSolverType constitutive_solver_type_;
    GeometrySolverType geometry_solver_type_;
    AnalysisSolverType analysis_solver_type_;

    std::unique_ptr<Matrix> global_stiffness_;
};
}


#endif //PROJECT_FEMSOLVER_H
