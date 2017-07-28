#ifndef PROJECT_TETRAHEDRONFINITEELEMENT_H
#define PROJECT_TETRAHEDRONFINITEELEMENT_H

#include <RTFEM/FEM/FiniteElement.h>

#include <memory>

namespace rtfem {

template<class T>
class Vertex;

template<class T>
class TetrahedronFiniteElement : public FiniteElement<T> {
public:
    TetrahedronFiniteElement(std::shared_ptr<Vertex<T>> vertex1,
                             std::shared_ptr<Vertex<T>> vertex2,
                             std::shared_ptr<Vertex<T>> vertex3,
                             std::shared_ptr<Vertex<T>> vertex4);

    ~TetrahedronFiniteElement();

    unsigned int GetVertexCount() const override;
private:
    const unsigned int vertex_count = 4;
};
}


#endif //PROJECT_TETRAHEDRONFINITEELEMENT_H
