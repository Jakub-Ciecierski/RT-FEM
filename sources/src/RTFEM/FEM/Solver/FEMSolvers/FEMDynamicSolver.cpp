#include "RTFEM/FEM/Solver/FEMSolvers/FEMDynamicSolver.h"

#include <RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMGlobalDynamicAssembler.h>
#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMFastForceAssembler.h>
#include "RTFEM/GPU/LinearSolver/GPULULinearSolver.cuh"
#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/DataStructure/Dense2SparseMatrix.h>

namespace rtfem {

template<class T>
FEMDynamicSolver<T>::FEMDynamicSolver(FEMModel<T>* fem_model,
                                      LinearSystemSolverType type) :
    FEMSolver<T>(fem_model),
    fem_assembler_data_(FEMGlobalAssemblerData<T>{0}),
    total_time_(0){
    solvers_ = LinearSystemSolvers<T>{
            type,
            rtfem::make_unique<GPUSparseCGLinearSolver<T>>(),
            rtfem::make_unique<GPUSparseCGPrecondLinearSolver<T>>(),
            rtfem::make_unique<GPULULinearSolver<T>>()};
}

template<class T>
FEMSolverOutput<T> FEMDynamicSolver<T>::Solve(){
    InitAssembly();
    auto global_mass = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
            fem_assembler_data_.global_mass_diagonal);
    InitDisplacementData();
    InitPreSolveLHS(global_mass);
    InitPreSolveRHS(global_mass);

    total_time_ = 0;

    return solver_output_;
}

template<class T>
void FEMDynamicSolver<T>::InitAssembly(){
    FEMGlobalDynamicAssembler<T> fem_assembler;
    fem_assembler_data_ = fem_assembler.Compute(*this->fem_model_);
}

template<class T>
void FEMDynamicSolver<T>::InitDisplacementData(){
    auto n = fem_assembler_data_.global_stiffness.rows();
    solver_output_.displacement = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);
    displacement_velocity_current_ = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);
    displacement_acceleration_current_ = Eigen::Vector<T,
            Eigen::Dynamic>::Zero(n);
}

template<class T>
void FEMDynamicSolver<T>::InitPreSolveLHS(
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& global_mass){
    T delta_time = 1.0 / 60.0;
    T delta_time_sqr = delta_time * delta_time;

    left_hand_side_bc_ =
            global_mass
            + delta_time * fem_assembler_data_.global_damping
            + delta_time_sqr * fem_assembler_data_.global_stiffness;
    left_hand_side_no_bc_ = left_hand_side_bc_;

    FEMGlobalDynamicAssembler<T> fem_assembler;
    fem_assembler.ApplyBoundaryConditionsMatrix(
            left_hand_side_bc_, this->fem_model_->boundary_conditions()
    );

    switch(solvers_.type) {
        case LinearSystemSolverType::CG: {
            Dense2SparseMatrix<T> dense2sparse;
            auto sparse_lhs = dense2sparse.Transform(left_hand_side_bc_,
                                                     MatrixType::General);
            solvers_.gpu_sparse_linear_solver_->PreSolve(sparse_lhs);
            break;
        }
        case LinearSystemSolverType::CG_PreCond:{
            Dense2SparseMatrix<T> dense2sparse;
            auto sparse_lhs = dense2sparse.Transform(left_hand_side_bc_,
                                                     MatrixType::General);
            solvers_.gpu_sparse_precond_linear_solver_->PreSolve(sparse_lhs);
            break;
        }
        case LinearSystemSolverType::LU: {
            solvers_.gpu_linear_solver_->PreSolve(left_hand_side_bc_.data(),
                                                  left_hand_side_bc_.rows());
            break;
        }
    }
}

template<class T>
void FEMDynamicSolver<T>::InitPreSolveRHS(
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& global_mass){
    Dense2SparseMatrix<T> dense2sparse;

    auto mass_sparse = dense2sparse.Transform(global_mass, MatrixType::General);
    gpu_mv_sparse_rhs_mass_.PreSolve(mass_sparse);

    auto stiffness_sparse =
            dense2sparse.Transform(fem_assembler_data_.global_stiffness,
                                   MatrixType::General);
    gpu_mv_sparse_rhs_stiffness_.PreSolve(stiffness_sparse);
}

template<class T>
void FEMDynamicSolver<T>::RunIteration(T delta_time){
    timer_ = FEMSolverTimer{};

    ReassembleForces();
    SolveForDisplacements(delta_time);
    ResetForces();

    solver_output_.timer = timer_;
}

template<class T>
void FEMDynamicSolver<T>::ReassembleForces(){
    timer_.Start();
    FEMFastForceAssembler<T> force_assembler;
    force_assembler.Assemble(this->fem_model_->fem_geometry(),
                             this->fem_model_->total_body_force(),
                             this->fem_model_->material(),
                             fem_assembler_data_.global_force);
    timer_.reassemble_forces_time = timer_.Stop();
}

template<class T>
void FEMDynamicSolver<T>::SolveForDisplacements(T delta_time){
    switch(this->type_){
        case FEMSolverType::GPU:{
            ImplicitNewtonGPU(delta_time);
            break;
        }
        case FEMSolverType::CPU:{
            ImplicitNewtonCPU(delta_time);
            break;
        }
    }

    total_time_ += delta_time;
}

template<class T>
void FEMDynamicSolver<T>::ResetForces(){
    timer_.Start();

    this->fem_model_->ResetBodyForce();
    this->fem_model_->ResetTractionForces();

    timer_.reset_force_time = timer_.Stop();
}

template<class T>
void FEMDynamicSolver<T>::ImplicitNewtonGPU(T delta_time){
    SolveRHSGPU(delta_time, fem_assembler_data_.global_force);
    auto& rhs = fem_assembler_data_.global_force;
    SolveBoundaryConditions(rhs);
    SolveLinearSystem(rhs, displacement_velocity_current_);
    SolveIntegration(delta_time, displacement_velocity_current_);
}

template<class T>
void FEMDynamicSolver<T>::SolveRHSGPU(
    T delta_time,
    Eigen::Vector<T, Eigen::Dynamic>& global_force){
    timer_.Start();
    gpu_mv_sparse_rhs_mass_.Solve(
            displacement_velocity_current_.data(),
            1.0,
            global_force.data(),
            delta_time);
    gpu_mv_sparse_rhs_stiffness_.Solve(
            solver_output_.displacement.data(),
            delta_time,
            global_force.data(),
            1.0);

    timer_.rhs_solve_time = timer_.Stop();
}

template<class T>
void FEMDynamicSolver<T>::SolveBoundaryConditions(
    Eigen::Vector<T, Eigen::Dynamic>& rhs){
    timer_.Start();
    FEMGlobalDynamicAssembler<T> fem_assembler;
    fem_assembler.ApplyBoundaryConditionsVector(
        left_hand_side_no_bc_, rhs, this->fem_model_->boundary_conditions()
    );
    timer_.boundary_conditions_solve_time = timer_.Stop();
}

template<class T>
void FEMDynamicSolver<T>::SolveLinearSystem(
    const Eigen::Vector<T, Eigen::Dynamic>& rhs,
    Eigen::Vector<T, Eigen::Dynamic>& velocity){
    timer_.Start();

    switch(solvers_.type){
        case LinearSystemSolverType::CG:
            solvers_.gpu_sparse_linear_solver_->Solve(
                    rhs.data(), velocity.data());
            break;
        case LinearSystemSolverType::CG_PreCond:
            solvers_.gpu_sparse_precond_linear_solver_->Solve(
                    rhs.data(), velocity.data());
            break;
        case LinearSystemSolverType::LU:
            solvers_.gpu_linear_solver_->Solve(
                    rhs.data(), rhs.size(), velocity.data());
            break;
    }

    timer_.cuda_solve_time = timer_.Stop();
}

template<class T>
void FEMDynamicSolver<T>::SolveIntegration(
    T delta_time,
    const Eigen::Vector<T, Eigen::Dynamic>& new_velocity){
    timer_.Start();
    const auto& initial_positions = fem_assembler_data_.global_position;
    auto current_positions = solver_output_.displacement + initial_positions;
    auto new_positions = current_positions
        + (delta_time * displacement_velocity_current_);
    solver_output_.displacement = new_positions - initial_positions;
    timer_.integration_solve_time = timer_.Stop();
}

template<class T>
void FEMDynamicSolver<T>::ImplicitNewtonCPU(T delta_time){
    auto global_mass = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
        fem_assembler_data_.global_mass_diagonal);
    Eigen::Vector<T, Eigen::Dynamic> rhs =
        delta_time * fem_assembler_data_.global_force
        + global_mass * displacement_velocity_current_
        + delta_time * (fem_assembler_data_.global_stiffness
            * solver_output_.displacement);

    SolveBoundaryConditions(rhs);

    displacement_velocity_current_ = this->SolveSystemOfEquations(
            left_hand_side_bc_, rhs);

    const auto& initial_positions = fem_assembler_data_.global_position;
    auto current_positions = solver_output_.displacement + initial_positions;
    auto new_positions = current_positions
        + (delta_time * displacement_velocity_current_);
    solver_output_.displacement = new_positions - initial_positions;
}

template
class FEMDynamicSolver<double>;
template
class FEMDynamicSolver<float>;

template
class LinearSystemSolvers<double>;
template
class LinearSystemSolvers<float>;

}