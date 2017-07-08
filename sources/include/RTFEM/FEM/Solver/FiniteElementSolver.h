#ifndef PROJECT_FINITEELEMENTSOLVER_H
#define PROJECT_FINITEELEMENTSOLVER_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/Matrix.h>

#include <memory>

namespace rtfem {

class FiniteElement;

 /**
  * Contains:
 *      Volume.
 *      Shape Matrix
 *      Geometry Matrix
 *
 *  Coordinates:
 *      x2 is assumed to point 'up'
 */
struct FiniteElementSolverData{
    FiniteElementSolverData() :
            volume(0),
            shape_matrix(0,0),
            geometry_matrix(0,0){}
    Float volume;

    Matrix shape_matrix;

    Matrix geometry_matrix;
};


/**
 * Computes FiniteElementSolverData for a given FiniteElement.
 */
class FiniteElementSolver {
public:
    FiniteElementSolver();

    virtual ~FiniteElementSolver();

    virtual FiniteElementSolverData Solve(std::shared_ptr<FiniteElement> finite_element) = 0;
};
}


#endif //PROJECT_FINITEELEMENTSOLVER_H
