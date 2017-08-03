#ifndef PROJECT_RANDOMDISTRIBUTIONS_H
#define PROJECT_RANDOMDISTRIBUTIONS_H

#include "RTFEM/FEM/Meshing/AABB.h"

#include <random>

namespace rtfem {

template<class T>
class RandomDistributions {
public:
    RandomDistributions(const AABB<T> &aabb);

    T GenerateX();
    T GenerateY();
    T GenerateZ();

private:
    std::mt19937 random_generator_;

    std::uniform_real_distribution<> x_distribution_;
    std::uniform_real_distribution<> y_distribution_;
    std::uniform_real_distribution<> z_distribution_;
};

}


#endif //PROJECT_RANDOMDISTRIBUTIONS_H
