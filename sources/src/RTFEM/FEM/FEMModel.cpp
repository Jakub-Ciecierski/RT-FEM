#include "RTFEM/FEM/FEMModel.h"

namespace rtfem {

FEMModel::FEMModel(std::vector<std::shared_ptr<rtfem::FiniteElement>>& finite_elements,
                   std::vector<std::shared_ptr<rtfem::Vertex>>& vertices,
                   const Material&& material) :
        finite_elements_(finite_elements),
        vertices_(vertices),
        material_(material){
}

FEMModel::~FEMModel() {}

UInt FEMModel::VertexCount(){
    return vertices_.size();
}

UInt FEMModel::FiniteElementCount(){
    return finite_elements_.size();
}

}