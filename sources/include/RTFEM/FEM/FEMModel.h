#ifndef PROJECT_FEM_H
#define PROJECT_FEM_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/FEM/Material.h>

#include <vector>
#include <memory>

namespace rtfem {

class Vertex;
class FiniteElement;

/**
 * Main class representing the FEM Model to be computed by FEM Solvers.
 *
 * Contains needed geometry (finite elements and vertices)
 * and constitutive coefficients (materials).
 *
 * FEMModel contains data of a single connected object.
 */
class FEMModel {
public:
    FEMModel(std::vector<std::shared_ptr<FiniteElement>>& finite_elements,
             std::vector<std::shared_ptr<Vertex>>& vertices,
             const Material&& material);
    ~FEMModel();

    const std::vector<std::shared_ptr<FiniteElement>>& finite_elements() const {
        return finite_elements_;
    }
    const std::vector<std::shared_ptr<Vertex>>& vertices() const {
        return vertices_;
    }

    Material& material(){return material_;}

    UInt VertexCount();
    UInt FiniteElementCount();

private:
    std::vector<std::shared_ptr<FiniteElement>> finite_elements_;
    std::vector<std::shared_ptr<Vertex>> vertices_;

    Material material_;
};

}


#endif //PROJECT_FEM_H
