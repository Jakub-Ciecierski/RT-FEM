#ifndef PROJECT_FEMDYNAMICSOLVER_H
#define PROJECT_FEMDYNAMICSOLVER_H

#include <RTFEM/FEM/Solver/FEMSolver.h>
#include <RTFEM/FEM/Solver/FEMGlobalAssembler.h>

namespace rtfem {

template<class T>
class FEMDynamicSolver : public FEMSolver<T>{
public:

    FEMDynamicSolver();
    ~FEMDynamicSolver() = default;

    const FEMSolverOutput<T>& solver_output(){return solver_output_;}

    T total_time(){return total_time_;}

    virtual FEMSolverOutput<T> Solve(const FEMModel<T> &fem_model) override;

    void RunIteration(T delta_time);

private:
    void ExplicitNewton(T delta_time);
    void ImplicitNewton(T delta_time);

    FEMSolverOutput<T> solver_output_;

    Eigen::Vector<T, Eigen::Dynamic> displacement_velocity_current_;
    Eigen::Vector<T, Eigen::Dynamic> displacement_acceleration_current_;

    FEMGlobalAssemblerData<T> fem_assembler_data_;

    T total_time_;
};
}


#endif //PROJECT_FEMDYNAMICSOLVER_H
