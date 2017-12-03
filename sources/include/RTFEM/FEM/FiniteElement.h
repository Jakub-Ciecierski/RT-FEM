#ifndef PROJECT_FINITEELEMENT_H
#define PROJECT_FINITEELEMENT_H

#include <RTFEM/FEM/FiniteElements/FiniteElementType.h>
#include <RTFEM/DataTypes.h>

#include <vector>
#include <memory>

namespace rtfem {

template<class T>
class Vertex;

/**
 * Abstract class for Finite Element.
 */
template<class T>
class FiniteElement {
public:
    FiniteElement(const FiniteElementType &&type);
    virtual ~FiniteElement() = default;

    const FiniteElementType &type() const { return type_; }

    const std::vector<unsigned int> &vertices_indices() {
        return
            vertices_indices_;
    }

    const std::vector<unsigned int>& faces_indices(){
        return faces_indices_;
    }

    T volume() { return volume_; }

    virtual unsigned int GetVertexCount() const = 0;
    virtual unsigned int GetFaceCount() const = 0;

    void SetVolume(T volume);
protected:
    std::vector<unsigned int> vertices_indices_;
    std::vector<unsigned int> faces_indices_;

    T volume_;
private:
    FiniteElementType type_;
};
}

#endif //PROJECT_FINITEELEMENT_H
