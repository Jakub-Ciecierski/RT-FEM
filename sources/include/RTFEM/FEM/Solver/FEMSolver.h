#ifndef PROJECT_FEMSOLVER_H
#define PROJECT_FEMSOLVER_H

#include <RTFEM/DataTypes.h>

#include <memory>
#include <vector>
#include <RTFEM/Timer.h>
#include <RTFEM/FEM/Solver/FEMGlobalAssembler.h>

namespace rtfem {

template<class T>
class FEMModel;

class FEMSolverTimer : public Timer {
public:
    double TotalTime() const {
        return reassemble_forces_time +
               rhs_solve_time +
               boundary_conditions_solve_time +
               cuda_solve_time +
               integration_solve_time +
               reset_force_time;
    }

    double reassemble_forces_time = 0;
    double rhs_solve_time = 0;
    double boundary_conditions_solve_time = 0;
    double cuda_solve_time = 0;
    double integration_solve_time = 0;
    double reset_force_time = 0;
};

template<class T>
struct FEMSolverOutput {
    Eigen::Vector<T, Eigen::Dynamic> displacement;

    FEMSolverTimer timer;
};

enum class FEMSolverType {
    CPU, GPU
};

template<class T>
class FEMSolver {
public:
    FEMSolver(FEMModel<T>* fem_model);
    virtual ~FEMSolver() = default;

    virtual FEMSolverOutput<T> Solve() = 0;

    void type(const FEMSolverType& type) { type_ = type; }
protected:
    Eigen::Vector<T, Eigen::Dynamic> SolveSystemOfEquations(
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
        const Eigen::Vector<T, Eigen::Dynamic>& b);

    FEMModel<T>* fem_model_;

    FEMSolverType type_;
};
}

#endif //PROJECT_FEMSOLVER_H
