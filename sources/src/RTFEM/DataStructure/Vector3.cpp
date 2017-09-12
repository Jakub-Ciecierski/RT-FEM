#include "RTFEM/DataStructure/Vector3.h"

namespace rtfem {
Vector3::Vector3() : x(0), y(0), z(0) {}

Vector3::Vector3(double x, double y, double z) : x(x), y(y), z(z) {}

Vector3::~Vector3() {}

bool Vector3::operator==(const Vector3 &rhs) const {
    return x == rhs.x &&
        y == rhs.y &&
        z == rhs.z;
}

bool Vector3::operator!=(const Vector3 &rhs) const {
    return !(rhs == *this);
}

}
