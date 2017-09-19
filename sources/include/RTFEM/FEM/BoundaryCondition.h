#ifndef PROJECT_BOUNDARYCONDITION_H
#define PROJECT_BOUNDARYCONDITION_H

#include <RTFEM/DataTypes.h>

namespace rtfem {

template<class T>
struct BoundaryCondition {
public:
    unsigned int vertex_id;
    Eigen::Vector3<T> value;
};

}

#endif //PROJECT_BOUNDARYCONDITION_H
