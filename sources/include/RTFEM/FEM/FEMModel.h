#ifndef PROJECT_FEM_MODEL_H
#define PROJECT_FEM_MODEL_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/FEM/Material.h>
#include <RTFEM/FEM/BoundaryCondition.h>

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
    FEMModel();

    FEMModel(std::unique_ptr<FEMGeometry<T>> fem_geometry,
             const Material<T> &&material);
    ~FEMModel() = default;

    const FEMGeometry<T> &fem_geometry() const { return *fem_geometry_; }

    void fem_geometry(std::unique_ptr<FEMGeometry<T>> fem_geometry) {
        fem_geometry_ = std::move(fem_geometry);
    }

    Material<T> &material() { return material_; }
    void material(Material<T> material) { material_ = material; }

    const std::vector<BoundaryCondition<T>> &boundary_conditions() const {
        return boundary_conditions_;
    }

    const Eigen::Vector3<T> &body_force() const {
        return body_force_;
    }

    bool AddBoundaryCondition(const BoundaryCondition<T> &boundary_condition);
    void RemoveBoundaryCondition(
        const BoundaryCondition<T> &boundary_condition);
    bool ExistsBoundaryCondition(
        const BoundaryCondition<T> &boundary_condition);
    /**
     * Sets Body Force (e.g. gravity) to the entire model.
     * Body Force is added to each finite element.
     *
     * @param body_force
     */
    void SetBodyForce(const Eigen::Vector3<T> &body_force);

private:
    std::unique_ptr<FEMGeometry<T>> fem_geometry_;

    Material<T> material_;

    std::vector<BoundaryCondition<T>> boundary_conditions_;

    Eigen::Vector3<T> body_force_;
};

}

#endif //PROJECT_FEM_MODEL_H
