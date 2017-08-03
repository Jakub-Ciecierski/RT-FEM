#ifndef PROJECT_DATAPRECISION_H
#define PROJECT_DATAPRECISION_H

#include <Eigen/Core>

namespace Eigen {

template<class T>
using Vector3 = Matrix<T, 3, 1>;

template<class T>
using Vector4 = Matrix<T, 4, 1>;

template<class T, int N>
using Vector = Matrix<T, N, 1>;

}

#endif //PROJECT_DATAPRECISION_H
