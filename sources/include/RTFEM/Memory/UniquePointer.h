#ifndef PROJECT_UNIQUEPTR_H
#define PROJECT_UNIQUEPTR_H

#include <memory>

// namespace makes it compatible with c++14
namespace rtfem {

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args &&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

}

#endif //PROJECT_UNIQUEPTR_H
