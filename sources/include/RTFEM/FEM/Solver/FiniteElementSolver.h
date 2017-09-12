#ifndef PROJECT_FINITEELEMENTSOLVER_H
#define PROJECT_FINITEELEMENTSOLVER_H

#include <RTFEM/DataTypes.h>

#include <memory>
#include <vector>

namespace rtfem {

template<class T>
class FiniteElement;

template<class T>
class Vertex;

/**
*  Contains:
*      Volume (V)
*      Geometry Matrix (B)
*      Force vector (p)
*
*  Coordinates:
*      x2 is assumed to point 'up'
*/
template<class T>
struct FiniteElementSolverData {
    FiniteElementSolverData() :
        volume(0) {}
    T volume;

    // Used to compute Global Stiffness
    //Matrix geometry_matrix;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> geometry_matrix;

    // Used to compute Global Force
    //Matrix force_vector;
    Eigen::Matrix<T, Eigen::Dynamic, 1> force_vector;
};

/**
 * Pressure forces directed along unit normal of i-th face
 * (traction_force_face1, i == 1)
 */
template<class T>
struct TractionForce {
    T force_face1 = 0;
    T force_face2 = 0;
    T force_face3 = 0;
    T force_face4 = 0;
};

/**
 * Computes FiniteElementSolverData for a given FiniteElement.
 */
template<class T>
class FiniteElementSolver {
public:
    FiniteElementSolver() = default;
    virtual ~FiniteElementSolver() = default;

    virtual FiniteElementSolverData<T> Solve(
        std::shared_ptr<FiniteElement<T>> finite_element,
        const std::vector<std::shared_ptr<Vertex<T>>> &vertices) = 0;

    virtual FiniteElementSolverData<T> Solve(
        std::shared_ptr<FiniteElement<T>> finite_element,
        const std::vector<std::shared_ptr<Vertex<T>>> &vertices,
        const Eigen::Vector3<T> &body_force,
        const TractionForce<T> &traction_force) = 0;
};
}

#endif //PROJECT_FINITEELEMENTSOLVER_H
