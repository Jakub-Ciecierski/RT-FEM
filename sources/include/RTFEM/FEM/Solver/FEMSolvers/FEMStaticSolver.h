#ifndef PROJECT_FEMLINEARSOLVER_H
#define PROJECT_FEMLINEARSOLVER_H

#include <RTFEM/FEM/Solver/FEMSolver.h>

namespace rtfem {

template<class T>
class FEMGlobalAssemblerData;

template<class T>
class FEMStaticSolver : public FEMSolver<T>{
public:
    FEMStaticSolver() = default;
    ~FEMStaticSolver() = default;

    virtual FEMSolverOutput<T> Solve(const FEMModel<T>& fem_model) override;
private:

};
}


#endif //PROJECT_FEMLINEARSOLVER_H
