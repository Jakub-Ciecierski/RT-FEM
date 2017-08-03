#include "RTFEM/FEM/Meshing/RandomDistributions.h"

namespace rtfem {

template<class T>
RandomDistributions<T>::RandomDistributions(const AABB<T> &aabb) :
        x_distribution_(std::uniform_real_distribution<>(aabb.Min()[0],
                                                        aabb.Max()[0])),
        y_distribution_(std::uniform_real_distribution<>(aabb.Min()[1],
                                                        aabb.Max()[1])),
        z_distribution_(std::uniform_real_distribution<>(aabb.Min()[2],
                                                        aabb.Max()[2])) {
    std::random_device random_device;
    random_generator_.seed(random_device());
}

template<class T>
T RandomDistributions<T>::GenerateX() {
    return x_distribution_(random_generator_);
}

template<class T>
T RandomDistributions<T>::GenerateY() {
    return y_distribution_(random_generator_);
}

template<class T>
T RandomDistributions<T>::GenerateZ() {
    return z_distribution_(random_generator_);
}

template class RandomDistributions<float>;
template class RandomDistributions<double>;

}