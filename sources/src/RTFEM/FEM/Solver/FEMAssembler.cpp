#include "RTFEM/FEM/Solver/FEMAssembler.h"

namespace rtfem {

FEMAssembler::FEMAssembler() {}

FEMAssembler::~FEMAssembler() {}

Matrix FEMAssembler::ComputeGlobalStiffness(const std::shared_ptr<FEMModel> fem_model) {
    return Matrix(0, 0);
}

}