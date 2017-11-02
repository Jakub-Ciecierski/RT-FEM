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
    Material() : young_modulus(0),
                 poisson_coefficient(0),
                 density(1),
                 damping1(0),
                 damping2(0) {}

    Material(T young_modulus_,
             T poisson_coefficient_) : young_modulus(young_modulus_),
                                       poisson_coefficient(
                                               poisson_coefficient_),
                                       density(1),
                                       damping1(0),
                                       damping2(0) {}

    Material(T young_modulus_,
             T poisson_coefficient_,
             T density_) : young_modulus(young_modulus_),
                           poisson_coefficient(poisson_coefficient_),
                           density(density_),
                           damping1(0),
                           damping2(0) {}

    T young_modulus;
    T poisson_coefficient;

    T density;

    T damping1;
    T damping2;
};

}

#endif //PROJECT_MATERIAL_H
