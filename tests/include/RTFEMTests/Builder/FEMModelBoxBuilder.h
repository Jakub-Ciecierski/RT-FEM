#ifndef PROJECT_FEMMODELBOXBUILDER_H
#define PROJECT_FEMMODELBOXBUILDER_H

#include <memory>

namespace rtfem {
template<class T>
class FEMModel;
}

class FEMModelBoxBuilder {
public:
    FEMModelBoxBuilder() = default;
    ~FEMModelBoxBuilder() = default;

    std::shared_ptr<rtfem::FEMModel<double>> Create();
private:
};


#endif //PROJECT_FEMMODELBOXBUILDER_H
