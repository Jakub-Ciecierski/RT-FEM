#include <RTFEM/FEM/Solver/FiniteElementSolver.h>

namespace rtfem {

template
class FiniteElementSolver<double>;
template
class FiniteElementSolver<float>;

template
struct TractionForce<double>;
template
struct TractionForce<float>;

template
struct FiniteElementSolverData<double>;
template
struct FiniteElementSolverData<float>;

}