#ifndef PROJECT_FEM_H
#define PROJECT_FEM_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/FEM/Material.h>

#include <vector>
#include <memory>

namespace rtfem {

template<class T>
class Vertex;

template<class T>
class FiniteElement;

/**
 * Main class representing the FEM Model to be computed by FEM Solvers.
 *
 * Contains needed geometry (finite elements and vertices)
 * and constitutive coefficients (materials).
 *
 * FEMModel contains data of a single connected object.
 */
template<class T>
class FEMModel {
public:
    FEMModel(std::vector<std::shared_ptr<FiniteElement<T>>>& finite_elements,
             std::vector<std::shared_ptr<Vertex<T>>>& vertices,
             const Material<T>&& material);
    ~FEMModel() = default;

    const std::vector<std::shared_ptr<FiniteElement<T>>>& finite_elements() const {
        return finite_elements_;
    }
    const std::vector<std::shared_ptr<Vertex<T>>>& vertices() const {
        return vertices_;
    }

    Material<T>& material(){return material_;}

    unsigned int VertexCount();
    unsigned int FiniteElementCount();

private:
    std::vector<std::shared_ptr<FiniteElement<T>>> finite_elements_;
    std::vector<std::shared_ptr<Vertex<T>>> vertices_;

    Material<T> material_;
};

}


#endif //PROJECT_FEM_H
