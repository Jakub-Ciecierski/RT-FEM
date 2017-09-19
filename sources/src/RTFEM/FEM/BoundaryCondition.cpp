#include "RTFEM/FEM/BoundaryCondition.h"

namespace rtfem {

template
struct BoundaryCondition<double>;
template
struct BoundaryCondition<float>;

}