#ifndef PROJECT_FEMSOLVER_H
#define PROJECT_FEMSOLVER_H

#include <RTFEM/FEM/Solver/FEMSolverTypes.h>
#include <RTFEM/DataTypes.h>

#include <memory>
#include <vector>

namespace rtfem {

template<class T>
class FEMModel;

template<class T>
class FEMAssemblerData;

template<class T>
struct FEMSolverOutput {
    Eigen::Vector<T, Eigen::Dynamic> displacement;
};

template<class T>
class FEMSolver {
public:
    FEMSolver();
    ~FEMSolver() = default;

    const ConstitutiveSolverType &constitutive_solver_type() const {
        return constitutive_solver_type_;
    }
    const GeometrySolverType &geometry_solver_type() const {
        return geometry_solver_type_;
    }
    const AnalysisSolverType &analysis_solver_type() const {
        return analysis_solver_type_;
    }

    void constitutive_solver_type(const ConstitutiveSolverType &type) {
        constitutive_solver_type_ = type;
    }
    void geometry_solver_type(const GeometrySolverType &type) {
        geometry_solver_type_ = type;
    }
    void analysis_solver_type(const AnalysisSolverType &type) {
        analysis_solver_type_ = type;
    }

    FEMSolverOutput<T> Solve(const FEMModel<T>& fem_model);

private:
    Eigen::Vector<T, Eigen::Dynamic> SolveSystemOfEquations(
        const FEMAssemblerData<T> &assembler_data);

    ConstitutiveSolverType constitutive_solver_type_;
    GeometrySolverType geometry_solver_type_;
    AnalysisSolverType analysis_solver_type_;
};
}

#endif //PROJECT_FEMSOLVER_H
