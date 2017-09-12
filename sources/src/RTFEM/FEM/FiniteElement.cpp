#include "RTFEM/FEM/FiniteElement.h"

namespace rtfem {

template<class T>
FiniteElement<T>::FiniteElement(const FiniteElementType &&type) : type_(type) {}

template
class FiniteElement<double>;
template
class FiniteElement<float>;

}
