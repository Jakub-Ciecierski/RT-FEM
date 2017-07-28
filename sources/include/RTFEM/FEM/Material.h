#ifndef PROJECT_MATERIAL_H
#define PROJECT_MATERIAL_H

#include <RTFEM/DataTypes.h>

namespace rtfem {

/**
 * One Material per FEMModel
 *
 * TODO:
 *  Bounds of material properties (e.g. poisson [0, 0.5]
 *
 * TODO:
 *  Material could be bound to each FiniteElement seperatly,
 *  allowing for 'illusion' of composite materials.
 *  That would require computing
 *  Constitutive Matrix for each FiniteElement.
 *
 */
template<class T>
struct Material {
    T young_modulus;
    T poisson_coefficient;

    T density;
};

}

#endif //PROJECT_MATERIAL_H
