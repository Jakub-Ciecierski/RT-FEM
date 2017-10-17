#include "RTFEM/FEM/FEMModel.h"

#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/BoundaryCondition.h>
#include <RTFEM/Memory/UniquePointer.h>

namespace rtfem {

template<class T>
FEMModel<T>::FEMModel() :
    material_(Material<T>{80000, 0.3, 1}),
    body_force_(Eigen::Vector3<T>(0, 0, 0)){}

template<class T>
FEMModel<T>::FEMModel(const FEMGeometry<T>& fem_geometry,
                      const Material<T> &&material) :
    fem_geometry_(fem_geometry),
    material_(material),
    body_force_(Eigen::Vector3<T>(0, 0, 0)){}

template<class T>
void FEMModel<T>::SetBodyForce(const Eigen::Vector3<T> &body_force){
    body_force_ = body_force;
}

template
class FEMModel<double>;
template
class FEMModel<float>;

}