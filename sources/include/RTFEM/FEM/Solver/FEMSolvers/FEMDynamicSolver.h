#ifndef PROJECT_FEMDYNAMICSOLVER_H
#define PROJECT_FEMDYNAMICSOLVER_H

#include <RTFEM/FEM/Solver/FEMSolver.h>
#include <RTFEM/FEM/Solver/FEMGlobalAssembler.h>
#include "RTFEM/GPU/LinearSolver/GPULULinearSolver.cuh"
#include "RTFEM/GPU/LinearSolver/GPUSparseCGLinearSolver.cuh"
#include "RTFEM/GPU/LinearSolver/GPUSparseCGPrecondLinearSolver.cuh"
#include "RTFEM/GPU/GPUMVMultiplication.cuh"
#include "RTFEM/GPU/GPUMVSparseMultiplication.cuh"

namespace rtfem {

enum class LinearSystemSolverType{
    LU, CG, CG_PreCond
};

template<class T>
struct LinearSystemSolvers{
    LinearSystemSolverType type;

    std::unique_ptr<GPUSparseLinearSolver<T>> gpu_sparse_linear_solver_;
    std::unique_ptr<GPUSparseLinearSolver<T>> gpu_sparse_precond_linear_solver_;
    std::unique_ptr<GPULULinearSolver<T>> gpu_linear_solver_;
};

/**
 * http://www.sciencedirect.com/science/article/pii/S0045794915001479
 * @tparam T
 */
template<class T>
class FEMDynamicSolver : public FEMSolver<T>{
public:

    FEMDynamicSolver(FEMModel<T>* fem_model,
                     LinearSystemSolverType type);
    ~FEMDynamicSolver() = default;

    const FEMSolverOutput<T>* solver_output(){return &solver_output_;}

    T total_time(){return total_time_;}

    const FEMSolverTimer& timer(){return timer_;}

    virtual FEMSolverOutput<T> Solve() override;

    void RunIteration(T delta_time);

private:
    void InitAssembly();
    void InitDisplacementData();
    void InitPreSolveLHS(
            const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& global_mass);
    void InitPreSolveRHS(
            const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& global_mass
    );

    void ReassembleForces();
    void SolveForDisplacements(T delta_time);
    void ResetForces();

    void ImplicitNewtonGPU(T delta_time);
    void SolveRHSGPU(T delta_time,
                     Eigen::Vector<T, Eigen::Dynamic>& global_force);
    void SolveBoundaryConditions(Eigen::Vector<T, Eigen::Dynamic>& rhs);
    void SolveLinearSystem(const Eigen::Vector<T, Eigen::Dynamic>& rhs,
                           Eigen::Vector<T, Eigen::Dynamic>& velocity);
    void SolveIntegration(
        T delta_time,
        const Eigen::Vector<T, Eigen::Dynamic>& new_velocity);

    void ImplicitNewtonCPU(T delta_time);

    LinearSystemSolvers<T> solvers_;

    GPUMVSparseMultiplication<T> gpu_mv_sparse_rhs_mass_;
    GPUMVSparseMultiplication<T> gpu_mv_sparse_rhs_stiffness_;

    FEMSolverOutput<T> solver_output_;

    // Saved for CPU solver only
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> left_hand_side_bc_;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> left_hand_side_no_bc_;

    Eigen::Vector<T, Eigen::Dynamic> displacement_velocity_current_;
    Eigen::Vector<T, Eigen::Dynamic> displacement_acceleration_current_;
    FEMGlobalAssemblerData<T> fem_assembler_data_;

    T total_time_;

    FEMSolverTimer timer_;
};
}


#endif //PROJECT_FEMDYNAMICSOLVER_H

