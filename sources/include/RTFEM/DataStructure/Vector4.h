#ifndef PROJECT_VECTOR4_H
#define PROJECT_VECTOR4_H

#include <RTFEM/DataTypes.h>

namespace rtfem {

struct Vector4 {
public:
    Vector4();
    Vector4(double x, double y, double z, double w);

    ~Vector4();

    double x;
    double y;
    double z;
    double w;
};
}

#endif //PROJECT_VECTOR4_H
