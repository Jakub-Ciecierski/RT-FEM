#include "RTFEM/FEM/BoundaryConditionContainer.h"

#include <RTFEM/FEM/BoundaryCondition.h>

namespace rtfem {

template<class T>
bool BoundaryConditionContainer<T>::AddBoundaryCondition(
    const BoundaryCondition<T> &boundary_condition) {
    if(ExistsBoundaryCondition(boundary_condition))
        return false;
    boundary_conditions_.push_back(boundary_condition);
    return true;
}

template<class T>
void BoundaryConditionContainer<T>::RemoveBoundaryCondition(
    const BoundaryCondition<T> &boundary_condition) {
    for (unsigned int i = 0; i < boundary_conditions_.size(); i++) {
        if (boundary_conditions_[i].vertex_id == boundary_condition.vertex_id) {
            boundary_conditions_.erase(boundary_conditions_.begin() + i);
        }
    }
}

template<class T>
bool BoundaryConditionContainer<T>::ExistsBoundaryCondition(
    const BoundaryCondition<T> &boundary_condition) {
    for (unsigned int i = 0; i < boundary_conditions_.size(); i++) {
        if (boundary_conditions_[i].vertex_id == boundary_condition.vertex_id) {
            return true;
        }
    }
    return false;
}

template<class T>
BoundaryCondition<T>* BoundaryConditionContainer<T>::begin(){
    return &(boundary_conditions_[0]);
}

template<class T>
BoundaryCondition<T>* BoundaryConditionContainer<T>::end(){
    return &(boundary_conditions_[boundary_conditions_.size()]);
}

template<class T>
const BoundaryCondition<T>* BoundaryConditionContainer<T>::begin() const{
    return &(boundary_conditions_[0]);
}

template<class T>
const BoundaryCondition<T>* BoundaryConditionContainer<T>::end() const {
    return &(boundary_conditions_[boundary_conditions_.size()]);
}

template
class BoundaryConditionContainer<double>;
template
class BoundaryConditionContainer<float>;

}