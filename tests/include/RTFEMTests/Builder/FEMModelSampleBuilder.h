#ifndef PROJECT_FEMMODELTESTBUILDER_H
#define PROJECT_FEMMODELTESTBUILDER_H

#include <RTFEM/DataTypes.h>

#include <memory>

namespace rtfem {
template<class T>
class FEMModel;
}
class FEMModelSampleBuilder {
public:
    FEMModelSampleBuilder() = default;
    ~FEMModelSampleBuilder() = default;

    const unsigned int finite_element_count() const { return finite_element_count_; }
    const unsigned int vertex_count() const { return vertex_count_; }

    std::shared_ptr<rtfem::FEMModel<double>> CreateRandomFEMModel();
private:
    const unsigned int finite_element_count_ = 4;
    const unsigned int vertex_count_ = 9;
};

#endif //PROJECT_FEMMODELTESTBUILDER_H
