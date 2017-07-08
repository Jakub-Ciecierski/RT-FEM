#ifndef PROJECT_FINITEELEMENT_H
#define PROJECT_FINITEELEMENT_H

#include <RTFEM/FEM/FiniteElements/FiniteElementType.h>
#include <RTFEM/DataTypes.h>

#include <vector>
#include <memory>

namespace rtfem {

class Vertex;

/**
 * Abstract class for Finite Element.
 */
class FiniteElement {
public:
    FiniteElement(const FiniteElementType&& type);
    virtual ~FiniteElement();

    const FiniteElementType& type() const {return type_;}

    const std::vector<std::shared_ptr<Vertex>>& vertices(){return vertices_;}

    virtual UInt GetVertexCount() const = 0;

protected:
    std::vector<std::shared_ptr<Vertex>> vertices_;
private:
    FiniteElementType type_;
};
}


#endif //PROJECT_FINITEELEMENT_H
