#ifndef PROJECT_FEMFASTFORCEASSEMBLER_H
#define PROJECT_FEMFASTFORCEASSEMBLER_H

#include <RTFEM/FEM/FEMGeometry.h>

namespace rtfem {

template<class T>
class Material;

template<class T>
class Vertex;

template<class T>
class FEMFastForceAssembler {
public:

    FEMFastForceAssembler() = default;
    ~FEMFastForceAssembler() = default;

    void Assemble(FEMGeometry<T>& fem_geometry,
                  const Eigen::Vector3<T>& body_force,
                  const Material<T>& material,
                  Eigen::Vector<T,Eigen::Dynamic>& global_force);
private:
    void ResetGlobalForce(
        Eigen::Vector<T, Eigen::Dynamic> &global_force);

    void AddTractionForces(
        FiniteElement<T>& finite_element,
        FEMGeometry<T>& fem_geometry,
        Eigen::Vector<T, Eigen::Dynamic> &global_force);
    void AddTractionForce(
        const TriangleFace<T>& triangle_face,
        Eigen::Vector<T, Eigen::Dynamic> &global_force);
    void AddTractionForceToVertex(
        unsigned int start_index,
        const T& x_value,
        const T& y_value,
        const T& z_value,
        Eigen::Vector<T, Eigen::Dynamic> &global_force);

    void AddBodyForce(
        FiniteElement<T>& finite_element,
        const Material<T>& material,
        const FEMGeometry<T>& fem_geometry,
        const Eigen::Vector3<T> &body_force,
        Eigen::Vector<T, Eigen::Dynamic> &global_force);
    void AddBodyForceToVertex(
        const Vertex<T> &vertex,
        const Eigen::Vector3<T> &body_force,
        T multiplier,
        Eigen::Vector<T, Eigen::Dynamic> &global_force);
};
}


#endif //PROJECT_FEMFASTFORCEASSEMBLER_H
