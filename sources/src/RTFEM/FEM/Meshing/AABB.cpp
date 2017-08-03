#include "RTFEM/FEM/Meshing/AABB.h"

namespace rtfem {

template<class T>
AABB<T>::AABB(const Eigen::Vector3<T>& min,
              const Eigen::Vector3<T>& max) :
        min_(min), max_(max){}

template struct AABB<float>;
template struct AABB<double>;

}