#include "RTFEM/FEM/FEMModel.h"

namespace rtfem {

template<class T>
FEMModel<T>::FEMModel(std::vector<std::shared_ptr<rtfem::FiniteElement<T>>>& finite_elements,
                   std::vector<std::shared_ptr<rtfem::Vertex<T>>>& vertices,
                   const Material<T>&& material) :
        finite_elements_(finite_elements),
        vertices_(vertices),
        material_(material){
}

template<class T>
unsigned int FEMModel<T>::VertexCount(){
    return vertices_.size();
}

template<class T>
unsigned int FEMModel<T>::FiniteElementCount(){
    return finite_elements_.size();
}

template class FEMModel<double>;
template class FEMModel<float>;

}