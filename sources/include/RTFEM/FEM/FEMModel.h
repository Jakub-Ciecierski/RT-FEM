#ifndef PROJECT_FEM_MODEL_H
#define PROJECT_FEM_MODEL_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/FEM/Material.h>

#include <vector>
#include <memory>

namespace rtfem {

template<class T>
class Vertex;

template<class T>
class FiniteElement;

template<class T>
struct FEMGeometry;

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
    FEMModel(std::unique_ptr<FEMGeometry<T>> fem_geometry,
             const Material <T> &&material);
    ~FEMModel() = default;

    const FEMGeometry<T> &fem_geometry() const { return *fem_geometry_; }

    Material <T> &material() { return material_; }

private:
    std::unique_ptr<FEMGeometry<T>> fem_geometry_;

    Material <T> material_;
};

}

#endif //PROJECT_FEM_MODEL_H
