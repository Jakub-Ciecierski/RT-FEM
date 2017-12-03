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
    void AddBodyForce(
            const Vertex<T>& vertex,
            const Eigen::Vector3<T>& body_force,
            T multiplier,
            Eigen::Vector<T, Eigen::Dynamic> &global_force);
};
}


#endif //PROJECT_FEMFASTFORCEASSEMBLER_H
