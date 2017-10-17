#include "RTFEM/FEM/FEMModel.h"

#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/BoundaryConditionContainer.h>
#include <RTFEM/FEM/BoundaryCondition.h>
#include <RTFEM/Memory/UniquePointer.h>

namespace rtfem {

template<class T>
FEMModel<T>::FEMModel() :
    fem_geometry_(nullptr),
    material_(Material<T>{80000, 0.3, 1}),
    boundary_conditions_(
        std::move(rtfem::make_unique<BoundaryConditionContainer<T>>())),
    body_force_(Eigen::Vector3<T>(0, 0, 0)){}

template<class T>
FEMModel<T>::FEMModel(std::unique_ptr<FEMGeometry<T>> fem_geometry,
                      const Material<T> &&material) :
    fem_geometry_(std::move(fem_geometry)),
    material_(material),
    boundary_conditions_(
        std::move(rtfem::make_unique<BoundaryConditionContainer<T>>())),
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