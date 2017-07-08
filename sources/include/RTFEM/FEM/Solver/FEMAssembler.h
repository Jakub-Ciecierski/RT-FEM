#ifndef PROJECT_FEMASSEMBLER_H
#define PROJECT_FEMASSEMBLER_H

#include <RTFEM/DataStructure/Matrix.h>
#include <memory>

namespace rtfem {

class FEMModel;

class FEMAssembler {
public:
    FEMAssembler();
    ~FEMAssembler();

    Matrix ComputeGlobalStiffness(const std::shared_ptr<FEMModel> fem_model);
};
}


#endif //PROJECT_FEMASSEMBLER_H
