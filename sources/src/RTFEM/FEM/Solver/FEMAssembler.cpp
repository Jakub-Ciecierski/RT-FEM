#include "RTFEM/FEM/Solver/FEMAssembler.h"

#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/Solver/FiniteElementSolver.h>
#include <RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver/TetrahedronSolver.h>
#include <RTFEM/FEM/FiniteElements/FiniteElementType.h>
#include <RTFEM/FEM/FiniteElement.h>
#include <RTFEM/FEM/Vertex.h>

#include <iostream>

namespace rtfem {

template<class T>
FEMAssemblerData<T> FEMAssembler<T>::Compute(const std::shared_ptr<FEMModel<T>> fem_model) {
    if(fem_model->finite_elements().size() == 0)
        throw std::invalid_argument("FEMModel contains no Finite Elements!");
    auto vertex_count = fem_model->VertexCount();
    auto global_dof_count = DIMENSION_COUNT*vertex_count;

    auto constitutive_matrix_C = ComputeConstitutiveMatrix(fem_model);

    FEMAssemblerData<T> fem_assembler_data(global_dof_count);
    for(const auto& finite_element : fem_model->finite_elements()){
        auto finite_element_solver = GetFiniteElementSolver(finite_element->type());
        auto finite_element_solver_data = finite_element_solver->Solve(finite_element);

        auto boolean_assembly_matrix_A = ComputeBooleanAssemblyMatrix(finite_element, vertex_count);

        auto partial_global_stiffness_matrix_Ke =
                ComputePartialGlobalStiffnessMatrix(finite_element_solver_data.geometry_matrix,
                                                    constitutive_matrix_C,
                                                    boolean_assembly_matrix_A,
                                                    finite_element_solver_data.volume);
        auto partial_global_force_vector_Q = ComputePartialGlobalForceVector(
                finite_element_solver_data.force_vector,
                boolean_assembly_matrix_A);

        fem_assembler_data.global_stiffness += partial_global_stiffness_matrix_Ke;
        fem_assembler_data.global_force += partial_global_force_vector_Q;
    }

    return fem_assembler_data;
}

template<class T>
Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N>
FEMAssembler<T>::ComputeConstitutiveMatrix(const std::shared_ptr<FEMModel<T>> fem_model){
    Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N>
            constitutive_matrix = Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N>::Zero();

    auto& material = fem_model->material();
    T p = material.poisson_coefficient;
    T e = material.young_modulus;
    T a = 1.0 - p;
    T b = (1.0 - 2.0*p) / 2.0;
    T c = e / ((1.0 + p) * (1.0 - 2.0*p));

    constitutive_matrix(0 ,0) = a;
    constitutive_matrix(0 ,1) = p;
    constitutive_matrix(0 ,2) = p;
    constitutive_matrix(1 ,0) = p;
    constitutive_matrix(1 ,1) = a;
    constitutive_matrix(1 ,2) = p;
    constitutive_matrix(2 ,0) = p;
    constitutive_matrix(2 ,1) = p;
    constitutive_matrix(2 ,2) = a;
    constitutive_matrix(3 ,3) = b;
    constitutive_matrix(4 ,4) = b;
    constitutive_matrix(5 ,5) = b;

    constitutive_matrix *= c;

    return constitutive_matrix;
}

template<class T>
std::unique_ptr<FiniteElementSolver<T>>
FEMAssembler<T>::GetFiniteElementSolver(const FiniteElementType& type){
    switch(type){
        case FiniteElementType::Tetrahedron:
            return rtfem::make_unique<TetrahedronSolver<T>>();
    }

    return nullptr;
}

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
FEMAssembler<T>::ComputeBooleanAssemblyMatrix(const std::shared_ptr<FiniteElement<T>> finite_element,
                                           unsigned int vertex_count) {
    auto local_vertex_count = finite_element->GetVertexCount();

    unsigned int n = DIMENSION_COUNT * local_vertex_count;
    unsigned int m = DIMENSION_COUNT * vertex_count;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A
            = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, m);

    for (unsigned int i = 0; i < local_vertex_count; i++) {
        const auto &vertex = finite_element->vertices()[i];
        auto global_id = vertex->id();

        auto global_column_start = global_id * DIMENSION_COUNT;
        auto local_row_start = i * DIMENSION_COUNT;

        A(local_row_start + 0, global_column_start + 0) = 1;
        A(local_row_start + 1, global_column_start + 1) = 1;
        A(local_row_start + 2, global_column_start + 2) = 1;
    }
    return A;
}

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
FEMAssembler<T>::ComputePartialGlobalStiffnessMatrix(
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &geometry_matrix,
        const Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N>& constitutive_matrix_C,
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& boolean_assembly_matrix_A,
        T volume) {
    auto local_stiffness_k = ComputeLocalStiffness(geometry_matrix,
                                                   constitutive_matrix_C,
                                                   volume);

    return boolean_assembly_matrix_A.transpose() * local_stiffness_k * boolean_assembly_matrix_A;
}

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
FEMAssembler<T>::ComputeLocalStiffness(
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &geometry_matrix,
        const Eigen::Matrix<T, CONSTITUTIVE_MATRIX_N, CONSTITUTIVE_MATRIX_N> &constitutive_matrix,
        T volume) {
    auto &C = constitutive_matrix;
    auto &B = geometry_matrix;
    auto BT = geometry_matrix.transpose();

    auto local_stiffness = volume * (BT * C * B);

    return local_stiffness;
}

template<class T>
Eigen::Vector<T, Eigen::Dynamic>
FEMAssembler<T>::ComputePartialGlobalForceVector(
        const Eigen::Vector<T, Eigen::Dynamic> &force_vector,
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &boolean_assembly_matrix_A) {
    return boolean_assembly_matrix_A.transpose() * force_vector;
}

template class FEMAssembler<double>;
template class FEMAssembler<float>;

template struct FEMAssemblerData<double>;
template struct FEMAssemblerData<float>;

}