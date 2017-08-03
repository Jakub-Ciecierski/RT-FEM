#ifndef PROJECT_AABB_H
#define PROJECT_AABB_H

#include <RTFEM/DataTypes.h>

namespace rtfem {

template<class T>
struct AABB {
public:
    AABB(const Eigen::Vector3<T>& min,
         const Eigen::Vector3<T>& max);
    ~AABB() = default;

    const Eigen::Vector3<T>& Min() const {return min_;}
    const Eigen::Vector3<T>& Max() const {return max_;}

private:
    Eigen::Vector3<T> min_;
    Eigen::Vector3<T> max_;
};

}

#endif //PROJECT_AABB_H
