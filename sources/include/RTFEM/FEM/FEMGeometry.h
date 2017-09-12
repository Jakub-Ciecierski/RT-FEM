#ifndef PROJECT_FEMGEOMETRY_H
#define PROJECT_FEMGEOMETRY_H

#include "vector"
#include "memory"

namespace rtfem {

template<class T>
class FiniteElement;

template<class T>
class Vertex;

/**
 * A list of vertices and finite elements that make up a 3D FEM Geometry
 */
template<class T>
struct FEMGeometry {
    std::vector<std::shared_ptr<Vertex<T>>> vertices;
    std::vector<std::shared_ptr<FiniteElement<T>>> finite_elements;
};
}


#endif //PROJECT_FEMGEOMETRY_H
