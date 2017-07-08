#include "RTFEM/FEM/FiniteElement.h"

namespace rtfem{

FiniteElement::FiniteElement(const FiniteElementType&& type) : type_(type){}
FiniteElement::~FiniteElement(){}

}
