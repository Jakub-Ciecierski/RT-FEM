#ifndef PROJECT_FEMSOLVERTYPES_H
#define PROJECT_FEMSOLVERTYPES_H

namespace rtfem{

/**
 * The Constitutive material model
 */
enum class ConstitutiveSolverType{
    LinearElastic
};

/**
 * The strain tensor representation.
 */
enum class GeometrySolverType{
    Linear, NonLinear
};

/**
 * Static or Dynamic Analysis
 */
enum class AnalysisSolverType{
    Static, Dynamic
};

}

#endif //PROJECT_FEMSOLVERTYPES_H
