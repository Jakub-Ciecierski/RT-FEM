#ifndef PROJECT_VECTORMATH_H
#define PROJECT_VECTORMATH_H

#include <RTFEM/DataTypes.h>

namespace rtfem {

class Vector3;

class VectorMath {
public:
    VectorMath();

    ~VectorMath();

    Vector3 Cross(const Vector3& v1, const Vector3& v2);

    Float Magnitude(const Vector3& v);
};
}


#endif //PROJECT_VECTORMATH_H
