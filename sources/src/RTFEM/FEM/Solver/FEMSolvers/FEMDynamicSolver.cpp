#include "RTFEM/FEM/Solver/FEMSolvers/FEMDynamicSolver.h"

#include <RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMGlobalDynamicAssembler.h>
#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMFastForceAssembler.h>
#include "RTFEM/GPU/LinearSolver/GPULinearSolver.cuh"
#include <RTFEM/Memory/UniquePointer.h>

#include <iostream>

namespace rtfem {

template<class T>
FEMDynamicSolver<T>::FEMDynamicSolver(FEMModel<T>* fem_model) :
    FEMSolver<T>(fem_model),
    fem_assembler_data_(FEMGlobalAssemblerData<T>{0}),
    total_time_(0){}

template<class T>
FEMSolverOutput<T> FEMDynamicSolver<T>::Solve(){
    FEMGlobalDynamicAssembler<T> fem_assembler;
    fem_assembler_data_ = fem_assembler.Compute(*this->fem_model_);

    auto n = fem_assembler_data_.global_stiffness.rows();
    solver_output_.displacement = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);
    displacement_velocity_current_ = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);
    displacement_acceleration_current_ = Eigen::Vector<T,
                                                       Eigen::Dynamic>::Zero(n);

    T delta_time = 1.0 / 60.0;
    T delta_time_sqr = delta_time * delta_time;
    auto global_mass = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
        fem_assembler_data_.global_mass_diagonal);
    left_hand_side_ =
        global_mass
            + delta_time * fem_assembler_data_.global_damping
            + delta_time_sqr * fem_assembler_data_.global_stiffness;

    gpu_linear_solver_.PreSolve(left_hand_side_.data(), n);

    gpu_multiplication_rhs_mass_.PreSolve(
        global_mass.data(),
        global_mass.rows());

    gpu_multiplication_rhs_stiffness_.PreSolve(
        fem_assembler_data_.global_stiffness.data(),
        fem_assembler_data_.global_stiffness.rows());

    total_time_ = 0;

    return solver_output_;
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
    timer_.Start();
    gpu_multiplication_rhs_mass_.Solve(
            displacement_velocity_current_.data(),
            1.0,
            fem_assembler_data_.global_force.data(),
            delta_time);
    gpu_multiplication_rhs_stiffness_.Solve(
            solver_output_.displacement.data(),
            delta_time,
            fem_assembler_data_.global_force.data(),
            1.0);
    auto& rhs = fem_assembler_data_.global_force;
    timer_.rhs_solve_time = timer_.Stop();

    timer_.Start();
    FEMGlobalDynamicAssembler<T> fem_assembler;
    fem_assembler.ApplyBoundaryConditions(
            left_hand_side_, rhs, this->fem_model_->boundary_conditions()
    );
    timer_.boundary_conditions_solve_time = timer_.Stop();

    timer_.Start();
    gpu_linear_solver_.Solve(rhs.data(), rhs.size(),
                             displacement_velocity_current_.data());
    timer_.cuda_solve_time = timer_.Stop();

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

    FEMGlobalDynamicAssembler<T> fem_assembler;
    fem_assembler.ApplyBoundaryConditions(
        left_hand_side_, rhs, this->fem_model_->boundary_conditions()
    );

    displacement_velocity_current_ = this->SolveSystemOfEquations(
            left_hand_side_, rhs);

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

}