#ifndef PROJECT_FINITEELEMENTSOLVER_H
#define PROJECT_FINITEELEMENTSOLVER_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/Matrix.h>

#include <memory>

namespace rtfem {

class FiniteElement;
class Vector3;

 /**
 *  Contains:
 *      Volume (V)
 *      Geometry Matrix (B)
 *      Force vector (p)
 *
 *  Coordinates:
 *      x2 is assumed to point 'up'
 */
struct FiniteElementSolverData{
    FiniteElementSolverData() :
            volume(0),
            geometry_matrix(0,0),
            force_vector(0,0){}
    Float volume;

    // Used to compute Global Stiffness
    Matrix geometry_matrix;

    // Used to compute Global Force
    Matrix force_vector;
};


/**
 * Computes FiniteElementSolverData for a given FiniteElement.
 */
class FiniteElementSolver {
public:
    FiniteElementSolver();

    virtual ~FiniteElementSolver();

    virtual FiniteElementSolverData Solve(std::shared_ptr<FiniteElement> finite_element) = 0;

    virtual FiniteElementSolverData Solve(std::shared_ptr<FiniteElement> finite_element,
                                          const Vector3& body_force,
                                          const Vector3& traction_force) = 0;
};
}


#endif //PROJECT_FINITEELEMENTSOLVER_H
