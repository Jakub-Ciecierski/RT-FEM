#include "RTFEM/FEM/Solver/FEMGlobalAssemblers/FEMGlobalDynamicAssembler.h"

#include <RTFEM/FEM/Material.h>
#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/Solver/FiniteElementSolver.h>
#include "RTFEM/GPU/GPUMMMultiplication.cuh"

namespace rtfem {

template<class T>
void FEMGlobalDynamicAssembler<T>::ComputeAssemblerData(
    FEMGlobalAssemblerData<T> &fem_assembler_data,
    FEMModel<T> &fem_model,
    Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N> &
    constitutive_matrix_C,
    bool force_only){
    FEMGlobalAssembler<T>::ComputeAssemblerData(fem_assembler_data,
                                                fem_model,
                                                constitutive_matrix_C,
                                                force_only);
    this->timer_.Start();
    if(!force_only){
        fem_assembler_data.global_damping =
                ComputeGlobalDampingMatrix(fem_assembler_data.global_mass,
                                           fem_assembler_data.global_stiffness,
                                           fem_model.material().damping_mass,
                                           fem_model.material().damping_stiffness);
    }

    this->timer_.partial_global_damping_time +=
            this->timer_.Stop();
    if(!force_only){
        AssembleGlobalPositionVector(fem_model.fem_geometry().vertices,
                                     fem_assembler_data.global_position);
    }
}

template<class T>
void FEMGlobalDynamicAssembler<T>::ComputeAssemblerDataIteration(
        FEMGlobalAssemblerData<T> &fem_assembler_data,
        const FiniteElementSolverData<T>& finite_element_solver_data,
        const FEMModel<T> &fem_model,
        const Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N> &
        constitutive_matrix_C,
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
        boolean_assembly_matrix_A,
        bool force_only) {
    FEMGlobalAssembler<T>::ComputeAssemblerDataIteration(
            fem_assembler_data,
            finite_element_solver_data,
            fem_model,
            constitutive_matrix_C,
            boolean_assembly_matrix_A,
            force_only);
    this->timer_.Start();
    if(!force_only){
        ComputePartialGlobalMassMatrix(fem_model.material().density,
                                       finite_element_solver_data.volume,
                                       boolean_assembly_matrix_A,
                                       fem_assembler_data.global_mass);
    }
    this->timer_.partial_global_mass_time +=
            this->timer_.Stop();
}

template<class T>
void
FEMGlobalDynamicAssembler<T>::ComputePartialGlobalMassMatrix(
        T density, T volume,
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
        boolean_assembly_matrix_A,
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& global_mass){
    auto local_mass_matrix = ComputeLocalMassMatrix(density, volume);
    Timer timer;
    timer.Start();

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic,
            Eigen::ColMajor> A_T = boolean_assembly_matrix_A.transpose();

    this->timer_.partial_global_mass_time_transpose += timer.Stop();

    int m = boolean_assembly_matrix_A.rows();
    int k = boolean_assembly_matrix_A.cols();
    int n = local_mass_matrix.cols();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> C
            = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(k, n);

    /**
     * A - k x m
     * B - m x n
     * C - k x n
     */
    timer.Start();
    GPUMMMultiplication<T> gpu_mm;
    gpu_mm.Solve(
            A_T.data(),
            local_mass_matrix.data(),
            C.data(),
            1, 0,
            k, m, n,
            MatrixOperation::None,
            MatrixOperation::None);
    this->timer_.partial_global_mass_time_cuda1 += timer.Stop();

/**
     * A - k x n
     * B - m x k
     * C - k x k
     */
    timer.Start();
    gpu_mm.Solve(
            C.data(),
            boolean_assembly_matrix_A.data(),
            global_mass.data(),
            1, 1,
            k, n, k,
            MatrixOperation::None,
            MatrixOperation::None);
    this->timer_.partial_global_mass_time_cuda2 += timer.Stop();
    /*
    return boolean_assembly_matrix_A.transpose() * local_mass_matrix
        * boolean_assembly_matrix_A;*/
};

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
FEMGlobalDynamicAssembler<T>::ComputeLocalMassMatrix(T density, T volume){
    Eigen::Matrix<T, 12, 12> mass_matrix;
    mass_matrix <<
                2,0,0,1,0,0,1,0,0,1,0,0,
                0,2,0,0,1,0,0,1,0,0,1,0,
                0,0,2,0,0,1,0,0,1,0,0,1,
                1,0,0,2,0,0,1,0,0,1,0,0,
                0,1,0,0,2,0,0,1,0,0,1,0,
                0,0,1,0,0,2,0,0,1,0,0,1,
                1,0,0,1,0,0,2,0,0,1,0,0,
                0,1,0,0,1,0,0,2,0,0,1,0,
                0,0,1,0,0,1,0,0,2,0,0,1,
                1,0,0,1,0,0,1,0,0,2,0,0,
                0,1,0,0,1,0,0,1,0,0,2,0,
                0,0,1,0,0,1,0,0,1,0,0,2;
    const T divider = 20;
    mass_matrix *= density * volume / divider;

    return mass_matrix;
};

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
FEMGlobalDynamicAssembler<T>::ComputeGlobalDampingMatrix(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
    global_mass_matrix,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
    global_stiffness_matrix,
    T damping_mass,
    T damping_stiffness){
    return
        global_mass_matrix * damping_mass +
        global_stiffness_matrix * damping_stiffness;
};

template<class T>
void FEMGlobalDynamicAssembler<T>::ApplyBoundaryConditionsToFEM(
    FEMGlobalAssemblerData<T> &assembler_data,
    const BoundaryConditionContainer<T> &boundary_conditions){

}

template<class T>
void FEMGlobalDynamicAssembler<T>::AssembleGlobalPositionVector(
    const std::vector<std::shared_ptr<Vertex<T>>>& vertices,
    Eigen::Vector<T, Eigen::Dynamic>& global_position){
    for(unsigned int i = 0; i < vertices.size(); i++){
        auto& vertex = vertices[i];
        global_position[i*DIMENSION_COUNT + 0] = vertex->x();
        global_position[i*DIMENSION_COUNT + 1] = vertex->y();
        global_position[i*DIMENSION_COUNT + 2] = vertex->z();
    }
};

template
class FEMGlobalDynamicAssembler<double>;
template
class FEMGlobalDynamicAssembler<float>;

}