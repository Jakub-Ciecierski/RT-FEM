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
    Vector3(double x, double y, double z);

    ~Vector3();

    bool operator==(const Vector3 &rhs) const;

    bool operator!=(const Vector3 &rhs) const;

    double x;
    double y;
    double z;
};
}

#endif //PROJECT_VECTOR3_H
