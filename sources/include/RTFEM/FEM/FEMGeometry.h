#ifndef PROJECT_FEMGEOMETRY_H
#define PROJECT_FEMGEOMETRY_H

#include "vector"
#include "memory"

namespace rtfem {

class FiniteElement;
class Vertex;

/**
 * A list of vertices and finite elements that make up a 3D FEM Geometry
 */
struct FEMGeometry {
    std::vector<std::shared_ptr<FiniteElement>> finite_elements;
    std::vector<std::shared_ptr<Vertex>> vertices;
};
}


#endif //PROJECT_FEMGEOMETRY_H
