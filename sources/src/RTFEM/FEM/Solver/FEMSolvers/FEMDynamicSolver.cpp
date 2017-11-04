#include "RTFEM/FEM/Solver/FEMSolvers/FEMDynamicSolver.h"

#include <RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMGlobalDynamicAssembler.h>

namespace rtfem {

template<class T>
FEMDynamicSolver<T>::FEMDynamicSolver(T delta_time) :
    fem_assembler_data_(FEMGlobalAssemblerData<T>{0}),
    delta_time_(delta_time),
    total_time_(0) {}

template<class T>
FEMSolverOutput<T> FEMDynamicSolver<T>::Solve(const FEMModel<T> &fem_model){
    FEMGlobalDynamicAssembler<T> fem_assembler;
    fem_assembler_data_ = fem_assembler.Compute(fem_model);

    auto n = fem_assembler_data_.global_stiffness.rows();
    solver_output_.displacement = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);
    displacement_velocity_current_ = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);
    displacement_acceleration_current_ = Eigen::Vector<T,
                                                       Eigen::Dynamic>::Zero(n);

    total_time_ = 0;

    return solver_output_;
}

template<class T>
void FEMDynamicSolver<T>::RunIteration(){
    ImplicitNewton();

    total_time_ += delta_time_;
}

template<class T>
void FEMDynamicSolver<T>::ExplicitNewton(){
    displacement_velocity_current_ = displacement_velocity_current_ +
        (delta_time_ * displacement_acceleration_current_);
    solver_output_.displacement = solver_output_.displacement +
        (delta_time_ * displacement_velocity_current_);

    auto right_hand_side =
        fem_assembler_data_.global_force
            - (fem_assembler_data_.global_damping * displacement_velocity_current_)
            - (fem_assembler_data_.global_stiffness * solver_output_.displacement);

    displacement_acceleration_current_ = this->SolveSystemOfEquations(
        fem_assembler_data_.global_mass,
        right_hand_side);
}

template<class T>
void FEMDynamicSolver<T>::ImplicitNewton(){
    auto delta_time_sqr = delta_time_ * delta_time_;

    auto lhs = fem_assembler_data_.global_mass
        + delta_time_ * fem_assembler_data_.global_damping
        + delta_time_sqr * fem_assembler_data_.global_stiffness;

    auto rhs = delta_time_ * fem_assembler_data_.global_force
        + fem_assembler_data_.global_mass * displacement_velocity_current_
        + delta_time_ * (fem_assembler_data_.global_stiffness
            * solver_output_.displacement);

    displacement_velocity_current_ = this->SolveSystemOfEquations(lhs, rhs);

    solver_output_.displacement = solver_output_.displacement
        + delta_time_ * displacement_velocity_current_;
}

template
class FEMDynamicSolver<double>;
template
class FEMDynamicSolver<float>;

}