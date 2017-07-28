#include "RTFEM/Math/VectorMath.h"

#include "RTFEM/DataStructure/Vector3.h"

#include <cmath>

namespace rtfem {

VectorMath::VectorMath() {}

VectorMath::~VectorMath() {}

Vector3 VectorMath::Cross(const Vector3& v1, const Vector3& v2){
    Vector3 vec(
            (v1.y * v2.z) - (v1.z * v2.y),
            (v1.z * v2.x) - (v1.x * v2.z),
            (v1.x * v2.y) - (v1.y * v2.x)
    );
    return vec;
}

double VectorMath::Magnitude(const Vector3 &v) {
    return std::sqrt(v.x * v.x + v.y*v.y + v.z * v.z);
}

}