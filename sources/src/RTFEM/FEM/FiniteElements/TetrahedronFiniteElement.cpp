#include "RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h"

namespace rtfem {

TetrahedronFiniteElement::TetrahedronFiniteElement(std::shared_ptr<Vertex> vertex1,
                                                   std::shared_ptr<Vertex> vertex2,
                                                   std::shared_ptr<Vertex> vertex3,
                                                   std::shared_ptr<Vertex> vertex4) :
        FiniteElement(std::move(FiniteElementType::Tetrahedron)){
    vertices_.resize(vertex_count);
    vertices_[0] = vertex1;
    vertices_[1] = vertex2;
    vertices_[2] = vertex3;
    vertices_[3] = vertex4;
}

TetrahedronFiniteElement::~TetrahedronFiniteElement() {}

UInt rtfem::TetrahedronFiniteElement::GetVertexCount() const {
    return 4;
}
}