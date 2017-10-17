#ifndef PROJECT_BOUNDARYCONDITIONCONTAINER_H
#define PROJECT_BOUNDARYCONDITIONCONTAINER_H

#include <vector>

namespace rtfem {

template <class T>
struct BoundaryCondition;

template<class T>
class BoundaryConditionContainer {
public:
    BoundaryConditionContainer() = default;
    ~BoundaryConditionContainer() = default;

    bool AddBoundaryCondition(
        const BoundaryCondition<T> &boundary_condition);
    void RemoveBoundaryCondition(
        const BoundaryCondition<T> &boundary_condition);
    bool ExistsBoundaryCondition(
        const BoundaryCondition<T> &boundary_condition);

    void Clear();

    unsigned int Size() const;

    BoundaryCondition<T>* begin();
    BoundaryCondition<T>* end();

    const BoundaryCondition<T>* begin() const;
    const BoundaryCondition<T>* end() const;

    const BoundaryCondition<T>& operator[](int i ) const;
    BoundaryCondition<T>& operator[](int i );
private:
    std::vector<BoundaryCondition<T>> boundary_conditions_;
};

}

#endif //PROJECT_BOUNDARYCONDITIONCONTAINER_H
