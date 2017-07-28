#include "RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h"

namespace rtfem {

template<class T>
TetrahedronFiniteElement<T>::TetrahedronFiniteElement(std::shared_ptr<Vertex<T>> vertex1,
                                                      std::shared_ptr<Vertex<T>> vertex2,
                                                      std::shared_ptr<Vertex<T>> vertex3,
                                                      std::shared_ptr<Vertex<T>> vertex4) :
        FiniteElement<T>(std::move(FiniteElementType::Tetrahedron)) {
    this->vertices_.resize(vertex_count);
    this->vertices_[0] = vertex1;
    this->vertices_[1] = vertex2;
    this->vertices_[2] = vertex3;
    this->vertices_[3] = vertex4;
}

template<class T>
TetrahedronFiniteElement<T>::~TetrahedronFiniteElement() {}

template<class T>
unsigned int TetrahedronFiniteElement<T>::GetVertexCount() const {
    return 4;
}

template class TetrahedronFiniteElement<double>;
template class TetrahedronFiniteElement<float>;

}