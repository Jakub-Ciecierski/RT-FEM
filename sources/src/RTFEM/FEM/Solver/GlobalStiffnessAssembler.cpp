#include "RTFEM/FEM/Solver/GlobalStiffnessAssembler.h"

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

GlobalStiffnessAssembler::GlobalStiffnessAssembler() {}

GlobalStiffnessAssembler::~GlobalStiffnessAssembler() {}

Matrix GlobalStiffnessAssembler::Compute(const std::shared_ptr<FEMModel> fem_model) {
    if(fem_model->finite_elements().size() == 0)
        throw std::invalid_argument("FEMModel contains no Finite Elements!");
    const UInt dim = 3;
    auto vertex_count = fem_model->VertexCount();
    auto global_stiffness_dimension = dim*vertex_count;

    Matrix global_stiffness_K(global_stiffness_dimension,
                              global_stiffness_dimension);

    auto constitutive_matrix_C = ComputeConstitutiveMatrix(fem_model);
    for(const auto& finite_element : fem_model->finite_elements()){
        auto local_stiffness_k = ComputeLocalStiffness(finite_element, constitutive_matrix_C);
        auto boolean_assembly_matrix_A = ComputeBooleanAssemblyMatrix(finite_element, vertex_count);
        auto partial_global_stiffness_matrix_Ke = ComputePartialGlobalStiffnessMatrix(boolean_assembly_matrix_A,
                                                                                      local_stiffness_k);
        global_stiffness_K = global_stiffness_K + partial_global_stiffness_matrix_Ke;
    }

    return global_stiffness_K;
}

Matrix GlobalStiffnessAssembler::ComputeConstitutiveMatrix(const std::shared_ptr<FEMModel> fem_model){
    Matrix constitutive_matrix(6,6);
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

Matrix
GlobalStiffnessAssembler::ComputeLocalStiffness(const std::shared_ptr<FiniteElement> finite_element,
                                                const Matrix &constitutive_matrix) {
    auto finite_element_solver = GetFiniteElementSolver(finite_element->type());
    auto finite_element_solver_data = finite_element_solver->Solve(finite_element);

    auto& C = constitutive_matrix;
    auto& B = finite_element_solver_data.geometry_matrix;
    auto BT = MatrixMath().Transpose(B);

    auto local_stiffness = BT * C * B;

    return local_stiffness;
}

std::unique_ptr<FiniteElementSolver>
GlobalStiffnessAssembler::GetFiniteElementSolver(const FiniteElementType& type){
    switch(type){
        case FiniteElementType::Tetrahedron:
            return rtfem::make_unique<TetrahedronSolver>();
    }

    return nullptr;
}

Matrix GlobalStiffnessAssembler::ComputeBooleanAssemblyMatrix(const std::shared_ptr<FiniteElement> finite_element,
                                                              UInt vertex_count){
    auto local_vertex_count = finite_element->GetVertexCount();

    const UInt dim = 3;
    Matrix A(dim*local_vertex_count, dim*vertex_count);

    for(UInt i = 0; i < local_vertex_count; i++){
        const auto& vertex = finite_element->vertices()[i];
        auto global_id = vertex->id();

        auto global_column_start = global_id * (Float)dim;
        auto local_row_start = i * (Float)dim;

        for(UInt k = 0; k < dim; k++){
            for(UInt l = 0; l < dim; l++){
                A[local_row_start + k][global_column_start + l] = 1;
            }
        }
    }
    return A;
}

Matrix GlobalStiffnessAssembler::ComputePartialGlobalStiffnessMatrix(const Matrix& boolean_assembly_matrix_A,
                                                                     const Matrix& local_stiffness_k){
    return MatrixMath().Transpose(boolean_assembly_matrix_A) * local_stiffness_k * boolean_assembly_matrix_A;
}

}