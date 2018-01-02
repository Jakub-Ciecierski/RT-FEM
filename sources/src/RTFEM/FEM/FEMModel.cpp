#include "RTFEM/FEM/FEMModel.h"

#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/BoundaryCondition.h>
#include <RTFEM/Memory/UniquePointer.h>

namespace rtfem {

template<class T>
FEMModel<T>::FEMModel(const FEMGeometry<T>& fem_geometry) :
    fem_geometry_(fem_geometry),
    material_(Material<T>{80000, 0.3, 1}),
    static_body_force_(Eigen::Vector3<T>::Zero()),
    dynamic_body_force_(Eigen::Vector3<T>::Zero()){}

template<class T>
void FEMModel<T>::SetStaticBodyForce(const Eigen::Vector3<T> &body_force){
    static_body_force_ = body_force;
    dynamic_body_force_ = static_body_force_;
}

template<class T>
void FEMModel<T>::AddDynamicBodyForce(const Eigen::Vector3<T> &body_force){
    dynamic_body_force_ += body_force;
}

template<class T>
void FEMModel<T>::ResetBodyForce(){
    dynamic_body_force_ = static_body_force_;
}

template<class T>
void FEMModel<T>::ResetTractionForces(){
    for(auto& triangle_face : fem_geometry_.triangle_faces){
        triangle_face.traction_force = triangle_face.constant_traction_force;
    }
}

template
class FEMModel<double>;
template
class FEMModel<float>;

}