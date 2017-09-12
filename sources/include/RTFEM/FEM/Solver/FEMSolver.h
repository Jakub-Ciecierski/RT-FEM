#ifndef PROJECT_FEMSOLVER_H
#define PROJECT_FEMSOLVER_H

#include <RTFEM/FEM/Solver/FEMSolverTypes.h>
#include <RTFEM/DataTypes.h>

#include <memory>

namespace rtfem {

template<class T>
class FEMModel;

template<class T>
class FEMAssemblerData;

template<class T>
class FEMSolver {
public:
    FEMSolver(
            const ConstitutiveSolverType&& constitutive_solver_type,
            const GeometrySolverType&& geometry_solver_type,
            const AnalysisSolverType&& analysis_solver_type);
    ~FEMSolver() = default;

    const ConstitutiveSolverType& constitutive_solver_type(){
        return constitutive_solver_type_;
    }
    const GeometrySolverType& geometry_solver_type(){
        return geometry_solver_type_;
    }
    const AnalysisSolverType& analysis_solver_type(){
        return analysis_solver_type_;
    }

    void Solve(const std::shared_ptr<FEMModel<T>> fem_model);

private:
    Eigen::Vector<T, Eigen::Dynamic> SolveSystemOfEquations(
            const FEMAssemblerData<T>& assembler_data);

    Eigen::Vector<T, Eigen::Dynamic> displacements_;

    ConstitutiveSolverType constitutive_solver_type_;
    GeometrySolverType geometry_solver_type_;
    AnalysisSolverType analysis_solver_type_;
};
}


#endif //PROJECT_FEMSOLVER_H
