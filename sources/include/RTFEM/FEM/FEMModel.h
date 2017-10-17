#ifndef PROJECT_FEM_MODEL_H
#define PROJECT_FEM_MODEL_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/FEM/Material.h>
#include <RTFEM/FEM/BoundaryConditionContainer.h>
#include <RTFEM/FEM/FEMGeometry.h>

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
    FEMModel(const FEMGeometry<T>& fem_geometry);
    ~FEMModel() = default;

    FEMGeometry<T> &fem_geometry() { return fem_geometry_; }
    const FEMGeometry<T> &fem_geometry() const { return fem_geometry_; }

    void fem_geometry(const FEMGeometry<T>& fem_geometry) {
        fem_geometry_ = fem_geometry;
    }

    const Material<T> &material() const { return material_; }
    void material(const Material<T>& material) { material_ = material; }

    BoundaryConditionContainer<T> &boundary_conditions() {
        return boundary_conditions_;
    }
    const BoundaryConditionContainer<T> &boundary_conditions() const {
        return boundary_conditions_;
    }
    void boundary_conditions(
        const BoundaryConditionContainer<T>& boundary_conditions) {
        boundary_conditions_ = boundary_conditions;
    }

    const Eigen::Vector3<T> &body_force() const {
        return body_force_;
    }

    /**
     * Sets Body Force (e.g. gravity) to the entire model.
     * Body Force is added to each finite element.
     *
     * @param body_force
     */
    void SetBodyForce(const Eigen::Vector3<T> &body_force);

private:
    FEMGeometry<T> fem_geometry_;

    Material<T> material_;

    BoundaryConditionContainer<T> boundary_conditions_;

    Eigen::Vector3<T> body_force_;
};

}

#endif //PROJECT_FEM_MODEL_H
