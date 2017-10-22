#include "RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h"

namespace rtfem {

template<class T>
TetrahedronFiniteElement<T>::TetrahedronFiniteElement(
    unsigned int vertex1,
    unsigned int vertex2,
    unsigned int vertex3,
    unsigned int vertex4) :
    FiniteElement<T>(std::move(FiniteElementType::Tetrahedron)) {

    this->vertices_indices_.resize(vertex_count);
    this->vertices_indices_[0] = vertex1;
    this->vertices_indices_[1] = vertex2;
    this->vertices_indices_[2] = vertex3;
    this->vertices_indices_[3] = vertex4;
}

template<class T>
TetrahedronFiniteElement<T>::TetrahedronFiniteElement(unsigned int vertex1,
                                                      unsigned int vertex2,
                                                      unsigned int vertex3,
                                                      unsigned int vertex4,
                                                      unsigned int face1,
                                                      unsigned int face2,
                                                      unsigned int face3,
                                                      unsigned int face4) :
    FiniteElement<T>(std::move(FiniteElementType::Tetrahedron)){
    this->vertices_indices_.resize(vertex_count);
    this->vertices_indices_[0] = vertex1;
    this->vertices_indices_[1] = vertex2;
    this->vertices_indices_[2] = vertex3;
    this->vertices_indices_[3] = vertex4;

    this->faces_indices_.resize(face_count);
    this->faces_indices_[0] = face1;
    this->faces_indices_[1] = face2;
    this->faces_indices_[2] = face3;
    this->faces_indices_[3] = face4;
}

template<class T>
unsigned int TetrahedronFiniteElement<T>::GetVertexCount() const {
    return vertex_count;
}

template<class T>
unsigned int TetrahedronFiniteElement<T>::GetFaceCount() const {
    return face_count;
}

template
class TetrahedronFiniteElement<double>;
template
class TetrahedronFiniteElement<float>;

}