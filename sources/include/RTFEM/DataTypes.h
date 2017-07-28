#ifndef PROJECT_DATAPRECISION_H
#define PROJECT_DATAPRECISION_H

#include <Eigen/Core>

namespace rtfem {

//using Float = double;

//using UInt = unsigned int;
//using Int = int;

}

namespace Eigen {

template<class T>
using Vector3 = Matrix<T, 3, 1>;

template<class T>
using Vector4 = Matrix<T, 4, 1>;

template<class T, int N>
using Vector = Matrix<T, N, 1>;

}

#endif //PROJECT_DATAPRECISION_H
