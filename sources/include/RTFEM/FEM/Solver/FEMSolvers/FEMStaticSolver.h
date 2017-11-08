#ifndef PROJECT_FEMLINEARSOLVER_H
#define PROJECT_FEMLINEARSOLVER_H

#include <RTFEM/FEM/Solver/FEMSolver.h>

namespace rtfem {

template<class T>
class FEMGlobalAssemblerData;

template<class T>
class FEMStaticSolver : public FEMSolver<T>{
public:
    FEMStaticSolver(FEMModel<T>* fem_model);
    ~FEMStaticSolver() = default;

    virtual FEMSolverOutput<T> Solve() override;
private:

};
}


#endif //PROJECT_FEMLINEARSOLVER_H
