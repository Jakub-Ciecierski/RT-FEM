#include "RTFEM/FEM/Vertex.h"

namespace rtfem {

template<class T>
Vertex<T>::Vertex(unsigned int id, const Eigen::Vector3<T> &coordinates)
        : id_(id), coordinates_(coordinates) {}

template class Vertex<double>;
template class Vertex<float>;

}