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
struct Material {
    Float young_modulus;
    Float poisson_coefficient;
};

}

#endif //PROJECT_MATERIAL_H
