#ifndef PROJECT_VECTOR4_H
#define PROJECT_VECTOR4_H

#include <RTFEM/DataTypes.h>

namespace rtfem {

struct Vector4 {
public:
    Vector4();
    Vector4(Float x, Float y, Float z, Float w);

    ~Vector4();

    Float x;
    Float y;
    Float z;
    Float w;
};
}


#endif //PROJECT_VECTOR4_H
