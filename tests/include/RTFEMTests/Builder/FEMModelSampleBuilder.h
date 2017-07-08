#ifndef PROJECT_FEMMODELTESTBUILDER_H
#define PROJECT_FEMMODELTESTBUILDER_H

#include <memory>
#include <RTFEM/DataTypes.h>

namespace rtfem {
class FEMModel;
}
class FEMModelSampleBuilder {
public:
    FEMModelSampleBuilder();

    ~FEMModelSampleBuilder();

    const rtfem::UInt finite_element_count() const {return finite_element_count_;}
    const rtfem::UInt vertex_count() const {return vertex_count_;}

    std::shared_ptr<rtfem::FEMModel> CreateRandomFEMModel();
private:
    const rtfem::UInt finite_element_count_ = 4;
    const rtfem::UInt vertex_count_ = 9;
};


#endif //PROJECT_FEMMODELTESTBUILDER_H
