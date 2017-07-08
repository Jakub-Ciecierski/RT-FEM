#ifndef PROJECT_TETRAHEDRONFINITEELEMENT_H
#define PROJECT_TETRAHEDRONFINITEELEMENT_H

#include <RTFEM/FEM/FiniteElement.h>

#include <memory>

namespace rtfem {

class Vertex;

class TetrahedronFiniteElement : public FiniteElement {
public:
    TetrahedronFiniteElement(std::shared_ptr<Vertex> vertex1,
                             std::shared_ptr<Vertex> vertex2,
                             std::shared_ptr<Vertex> vertex3,
                             std::shared_ptr<Vertex> vertex4);

    ~TetrahedronFiniteElement();

    UInt GetVertexCount() const override;
private:
    const UInt vertex_count = 4;
};
}


#endif //PROJECT_TETRAHEDRONFINITEELEMENT_H
