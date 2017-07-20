#ifndef PROJECT_VECTOR3_H
#define PROJECT_VECTOR3_H

#include <RTFEM/DataTypes.h>

namespace rtfem {

/**
 * Three dimensional vector
 */
struct Vector3 {
public:
    Vector3();
    Vector3(Float x, Float y, Float z);

    ~Vector3();

    bool operator==(const Vector3 &rhs) const;

    bool operator!=(const Vector3 &rhs) const;

    Float x;
    Float y;
    Float z;
};
}


#endif //PROJECT_VECTOR3_H
