#include "RTFEM/FEM/Solver/FEMAssembler.h"

#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/DataStructure/Matrix.h>
#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/FEM/Solver/FiniteElementSolver.h>
#include <RTFEM/FEM/FiniteElements/FiniteElementType.h>
#include <RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver.h>
#include <RTFEM/FEM/FiniteElement.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/Math/MatrixMath.h>

namespace rtfem {

constexpr UInt dimensions = 3;

FEMAssembler::FEMAssembler() {}

FEMAssembler::~FEMAssembler() {}

FEMAssemblerData FEMAssembler::Compute(const std::shared_ptr<FEMModel> fem_model) {
    if(fem_model->finite_elements().size() == 0)
        throw std::invalid_argument("FEMModel contains no Finite Elements!");
    auto vertex_count = fem_model->VertexCount();
    auto global_dof_count = dimensions*vertex_count;

    auto constitutive_matrix_C = ComputeConstitutiveMatrix(fem_model);

    FEMAssemblerData fem_assembler_data(global_dof_count);
    for(const auto& finite_element : fem_model->finite_elements()){
        auto finite_element_solver = GetFiniteElementSolver(finite_element->type());
        auto finite_element_solver_data = finite_element_solver->Solve(finite_element);

        auto boolean_assembly_matrix_A = ComputeBooleanAssemblyMatrix(finite_element, vertex_count);

        auto partial_global_stiffness_matrix_Ke =
                ComputePartialGlobalStiffnessMatrix(finite_element_solver_data.geometry_matrix,
                                                    constitutive_matrix_C,
                                                    boolean_assembly_matrix_A,
                                                    finite_element_solver_data.volume);
        auto partial_global_force_vector_Q = ComputePartialGlobalForceVector(finite_element_solver_data.force_vector,
                                                                             boolean_assembly_matrix_A);

        fem_assembler_data.global_stiffness += partial_global_stiffness_matrix_Ke;
        fem_assembler_data.global_force += partial_global_force_vector_Q;
    }

    return fem_assembler_data;
}

Matrix FEMAssembler::ComputeConstitutiveMatrix(const std::shared_ptr<FEMModel> fem_model){
    constexpr UInt constitutive_matrix_size = 6;
    Matrix constitutive_matrix(constitutive_matrix_size, constitutive_matrix_size);
    auto& material = fem_model->material();

    Float p = material.poisson_coefficient;
    Float e = material.young_modulus;

    Float a = (Float)1.0 - p;
    Float b = ((Float)1.0 - (Float)2.0*p) / (Float)2.0;
    Float c = e / (((Float)1.0 + p) * ((Float)1 - (Float)2*p));

    constitutive_matrix[0][0] = a;
    constitutive_matrix[0][1] = p;
    constitutive_matrix[0][2] = p;
    constitutive_matrix[1][0] = p;
    constitutive_matrix[1][1] = a;
    constitutive_matrix[1][2] = p;
    constitutive_matrix[2][0] = p;
    constitutive_matrix[2][1] = p;
    constitutive_matrix[2][2] = a;
    constitutive_matrix[3][3] = b;
    constitutive_matrix[4][4] = b;
    constitutive_matrix[5][5] = b;

    constitutive_matrix = constitutive_matrix * c;

    return constitutive_matrix;
}

std::unique_ptr<FiniteElementSolver>
FEMAssembler::GetFiniteElementSolver(const FiniteElementType& type){
    switch(type){
        case FiniteElementType::Tetrahedron:
            return rtfem::make_unique<TetrahedronSolver>();
    }

    return nullptr;
}

Matrix FEMAssembler::ComputeBooleanAssemblyMatrix(const std::shared_ptr<FiniteElement> finite_element,
                                                  UInt vertex_count) {
    auto local_vertex_count = finite_element->GetVertexCount();

    Matrix A(dimensions * local_vertex_count, dimensions * vertex_count);

    for (UInt i = 0; i < local_vertex_count; i++) {
        const auto &vertex = finite_element->vertices()[i];
        auto global_id = vertex->id();

        auto global_column_start = global_id * dimensions;
        auto local_row_start = i * dimensions;

        A[local_row_start + 0][global_column_start + 0] = 1;
        A[local_row_start + 1][global_column_start + 1] = 1;
        A[local_row_start + 2][global_column_start + 2] = 1;
    }
    return A;
}

Matrix FEMAssembler::ComputePartialGlobalStiffnessMatrix(const Matrix &geometry_matrix,
                                                         const Matrix &constitutive_matrix_C,
                                                         const Matrix &boolean_assembly_matrix_A,
                                                         Float volume) {
    auto local_stiffness_k = ComputeLocalStiffness(geometry_matrix,
                                                   constitutive_matrix_C,
                                                   volume);
    return MatrixMath().Transpose(boolean_assembly_matrix_A) * local_stiffness_k * boolean_assembly_matrix_A;
}

Matrix
FEMAssembler::ComputeLocalStiffness(const Matrix &geometry_matrix,
                                    const Matrix &constitutive_matrix,
                                    Float volume) {
    auto &C = constitutive_matrix;
    auto &B = geometry_matrix;
    auto BT = MatrixMath().Transpose(B);

    auto local_stiffness = volume * (BT * C * B);

    return local_stiffness;
}

Matrix FEMAssembler::ComputePartialGlobalForceVector(const Matrix &force_vector,
                                                     const Matrix &boolean_assembly_matrix_A) {
    MatrixMath matrix_math;
    return matrix_math.Transpose(boolean_assembly_matrix_A) * force_vector;
}

}