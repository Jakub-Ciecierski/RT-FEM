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
struct FEMSolverOutput {
    Eigen::Vector<T, Eigen::Dynamic> displacement;
};

template<class T>
class FEMSolver {
public:
    FEMSolver() = default;
    virtual ~FEMSolver() = default;

    virtual FEMSolverOutput<T> Solve(const FEMModel<T>& fem_model) = 0;

private:
};
}

#endif //PROJECT_FEMSOLVER_H
