#include "RTFEM/FEM/FEMModel.h"

#include <RTFEM/FEM/FEMGeometry.h>

namespace rtfem {

template<class T>
FEMModel<T>::FEMModel(std::unique_ptr<FEMGeometry<T>> fem_geometry,
                      const Material<T>&& material) :
        fem_geometry_(std::move(fem_geometry)),
        material_(material){}

template class FEMModel<double>;
template class FEMModel<float>;

}