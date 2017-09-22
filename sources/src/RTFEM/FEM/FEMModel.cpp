#include "RTFEM/FEM/FEMModel.h"

#include <RTFEM/FEM/FEMGeometry.h>

namespace rtfem {

template<class T>
FEMModel<T>::FEMModel(std::unique_ptr<FEMGeometry<T>> fem_geometry,
                      const Material<T> &&material) :
    fem_geometry_(std::move(fem_geometry)),
    material_(material),
    body_force_(Eigen::Vector3<T>(0, 0, 0)) {}

template<class T>
void FEMModel<T>::AddBoundaryCondition(
    const BoundaryCondition<T> &boundary_condition) {
    boundary_conditions_.push_back(boundary_condition);
}

template
class FEMModel<double>;
template
class FEMModel<float>;

}