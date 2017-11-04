#ifndef PROJECT_FEMDYNAMICSOLVER_H
#define PROJECT_FEMDYNAMICSOLVER_H

#include <RTFEM/FEM/Solver/FEMSolver.h>
#include <RTFEM/FEM/Solver/FEMGlobalAssembler.h>

namespace rtfem {

template<class T>
class FEMDynamicSolver : public FEMSolver<T>{
public:

    FEMDynamicSolver(T delta_time);
    ~FEMDynamicSolver() = default;

    const FEMSolverOutput<T>& solver_output(){return solver_output_;}

    T delta_time(){return delta_time_;}
    T total_time(){return total_time_;}

    virtual FEMSolverOutput<T> Solve(const FEMModel<T> &fem_model) override;

    void RunIteration();

private:
    void ExplicitNewton();
    void ImplicitNewton();

    FEMSolverOutput<T> solver_output_;

    Eigen::Vector<T, Eigen::Dynamic> displacement_velocity_current_;
    Eigen::Vector<T, Eigen::Dynamic> displacement_acceleration_current_;

    FEMGlobalAssemblerData<T> fem_assembler_data_;

    T delta_time_;
    T total_time_;
};
}


#endif //PROJECT_FEMDYNAMICSOLVER_H
