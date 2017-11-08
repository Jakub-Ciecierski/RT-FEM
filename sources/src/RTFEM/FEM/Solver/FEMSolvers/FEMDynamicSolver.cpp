#include "RTFEM/FEM/Solver/FEMSolvers/FEMDynamicSolver.h"

#include <RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMGlobalDynamicAssembler.h>
#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/FEMModel.h>
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

    total_time_ = 0;

    return solver_output_;
}

template<class T>
void FEMDynamicSolver<T>::RunIteration(T delta_time){
    ReassembleForces();
    SolveForDisplacements(delta_time);
    ResetForces();
}

template<class T>
void FEMDynamicSolver<T>::ReassembleForces(){
    FEMGlobalDynamicAssembler<T> fem_assembler;
    fem_assembler_data_ = fem_assembler.Compute(*this->fem_model_);
}

template<class T>
void FEMDynamicSolver<T>::SolveForDisplacements(T delta_time){
    ImplicitNewton(delta_time);
    total_time_ += delta_time;
}

template<class T>
void FEMDynamicSolver<T>::ResetForces(){
    this->fem_model_->ResetBodyForce();
    this->fem_model_->ResetTractionForces();
}

template<class T>
void FEMDynamicSolver<T>::ExplicitNewton(T delta_time){
    displacement_velocity_current_ = displacement_velocity_current_ +
        (delta_time * displacement_acceleration_current_);
    solver_output_.displacement = solver_output_.displacement +
        (delta_time * displacement_velocity_current_);

    auto right_hand_side =
        fem_assembler_data_.global_force
            - (fem_assembler_data_.global_damping * displacement_velocity_current_)
            - (fem_assembler_data_.global_stiffness * solver_output_.displacement);

    displacement_acceleration_current_ = this->SolveSystemOfEquations(
        fem_assembler_data_.global_mass,
        right_hand_side);
}

template<class T>
void FEMDynamicSolver<T>::ImplicitNewton(T delta_time){
    auto delta_time_sqr = delta_time * delta_time;

    auto lhs = fem_assembler_data_.global_mass
        + delta_time * fem_assembler_data_.global_damping
        + delta_time_sqr * fem_assembler_data_.global_stiffness;

    auto rhs = delta_time * fem_assembler_data_.global_force
        + fem_assembler_data_.global_mass * displacement_velocity_current_
        + delta_time * (fem_assembler_data_.global_stiffness
            * solver_output_.displacement);

    displacement_velocity_current_ = this->SolveSystemOfEquations(lhs, rhs);

    solver_output_.displacement = solver_output_.displacement
        + delta_time * displacement_velocity_current_;
}

template
class FEMDynamicSolver<double>;
template
class FEMDynamicSolver<float>;

}