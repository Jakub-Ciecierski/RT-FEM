#ifndef PROJECT_TETRAHEDRONSOLVER_H
#define PROJECT_TETRAHEDRONSOLVER_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/Matrix.h>
#include <memory>
#include <RTFEM/DataStructure/Vector3.h>
#include <RTFEM/FEM/Solver/FiniteElementSolver.h>

namespace rtfem {

class Vertex;
class Vector4;

class TetrahedronFiniteElement;

/**
 * fi(x) = (Ai + Bi*x1 + Ci*x2 + Di*x3) / 6*V
 */
struct TetrahedronShapeFunctionCoefficients{
    // A1 = 6*V0i
    Float A1;
    Float A2;
    Float A3;
    Float A4;

    Float B1;
    Float B2;
    Float B3;
    Float B4;

    Float C1;
    Float C2;
    Float C3;
    Float C4;

    Float D1;
    Float D2;
    Float D3;
    Float D4;
};

/**
 * xij = xi - xj
 */
struct Edges {
    Float x32 = 0;
    Float x34 = 0;
    Float x43 = 0;
    Float x14 = 0;
    Float x21 = 0;
    Float x31 = 0;
    Float x24 = 0;
    Float x42 = 0;
    Float x13 = 0;
    Float x12 = 0;

    Float z43 = 0;
    Float z31 = 0;
    Float z32 = 0;
    Float z24 = 0;
    Float z34 = 0;
    Float z13 = 0;
    Float z14 = 0;
    Float z21 = 0;
    Float z42 = 0;
    Float z12 = 0;

    Float y42 = 0;
    Float y31 = 0;
    Float y24 = 0;
    Float y13 = 0;
    Float y32 = 0;
    Float y34 = 0;
    Float y14 = 0;
    Float y12 = 0;
    Float y43 = 0;
    Float y21 = 0;
};

struct FacesArea {
    Float area1;
    Float area2;
    Float area3;
    Float area4;
};

/**
 * Solver for Linear Tetrahedron (constant gradient of shape function).
 * The geometry matrix B is constant with respect to X
 * (and always will be constant, independent of material/strain equations)
 * .
 * The simplicity of Linear Tetrahedron allows for analytical integration.
 *
 * Solver Data:
 *
 *  Geometry matrix: [6 x 12] (Shape function gradient)
 *      Used in computing Stiffness Matrix
 *
 *  Force vector: [12x1]
 *      Body force (e.g. gravity) (f)
 *      Traction force (e.g. collision) (t)
 *
 * NOTE:
 *  1) Special care must be taken for Local vertex ordering (TetrahedronSolver does NOT do any ordering)
 *      Swap local ordering of vertices
 *      Pick a vertex as initial one 1.
 *      Pick a face containing the first three vertices 1,2,3 (including the initial one).
 *      The vertex 4 is the one opposite of chosen face.
 *      Number the first three vertices in counterclockwise wise order as seen from vertex 4
 *
 *  2) Research Materials:
 *      http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
 */
class TetrahedronSolver : public FiniteElementSolver{
public:
    TetrahedronSolver();

    ~TetrahedronSolver();

    virtual FiniteElementSolverData Solve(std::shared_ptr<FiniteElement> finite_element) override;

    virtual FiniteElementSolverData Solve(std::shared_ptr<FiniteElement> finite_element,
                                          const Vector3& body_force,
                                          const TractionForce& traction_force) override;

    Matrix SolveJacobianInverse(std::shared_ptr<FiniteElement> finite_element);

private:
    Edges ComputeEdgesCache(const Vertex &v1, const Vertex &v2,
                                 const Vertex &v3, const Vertex &v4);

    /**
     * Computes The coefficients on linear tetrahedron shape function.
     * fi(x) = (Ai + Bi*x1 + Ci*x2 + Di*x3) / 6*V
     *
     * @param v1
     * @param v2
     * @param v3
     * @param v4
     * @return
     */
    TetrahedronShapeFunctionCoefficients
    ComputeShapeFunctionCoefficients(const Vertex& v1, const Vertex& v2,
                                     const Vertex& v3, const Vertex& v4,
                                     const Edges& edges);
    Vector4 ComputeAi(const Vertex &v1, const Vertex &v2,
                      const Vertex &v3, const Vertex &v4);

    /**
     * Computes The [6x12] Geometric Matrix (B)
     *
     * @param coefficients
     * @param volume
     * @return
     */
    Matrix ComputeGeometryMatrix(const TetrahedronShapeFunctionCoefficients& coefficients, Float volume);
    void AssemblyGeometryMatrix(Matrix& B, UInt column_offset,
                                Float b, Float c, Float d);

    Float ComputeVolume(const TetrahedronShapeFunctionCoefficients& coefficients);

    Matrix AssemblyJacobianInverse(const TetrahedronShapeFunctionCoefficients& coefficients,
                                   Float volume);

    /**
     * Computes Force vector [12x1] which is a combination of body force and traction force.
     *
     * @param shape_function_coefficients
     * @param body_force
     * @param traction_force
     * @return
     */
    Matrix ComputeForceVector(const TetrahedronShapeFunctionCoefficients& shape_function_coefficients,
                              Float volume, const FacesArea& faces_area,
                              const Vector3& body_force,
                              const TractionForce& traction_force);

    Matrix ComputeBodyForceVector(Float volume,
                                  const Vector3& body_force);

    Matrix ComputeTractionForceVector(const TetrahedronShapeFunctionCoefficients& shape_function_coefficients,
                                      const FacesArea& faces_area,
                                      const TractionForce& traction_force);

    Matrix ComputeTractionForceVectorFace(UInt face_index, Float traction_force,
                                          Float area,
                                          Float B, Float C, Float D);
    FacesArea ComputeFacesArea(const Edges& edges);
};
}


#endif //PROJECT_TETRAHEDRONSOLVER_H
