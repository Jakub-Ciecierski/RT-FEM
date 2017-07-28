#ifndef PROJECT_TETRAHEDRONSOLVERDATA_H
#define PROJECT_TETRAHEDRONSOLVERDATA_H

namespace rtfem {

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

constexpr int TETRAHEDRON_DOF_COUNT = 12;
constexpr int TETRAHEDRON_JACOBIAN_MATRIX_N = 4;
constexpr int TETRAHEDRON_GEOMETRIC_MATRIX_N = 6;
constexpr int TETRAHEDRON_GEOMETRIC_MATRIX_M = TETRAHEDRON_DOF_COUNT;
constexpr int TETRAHEDRON_FORCE_VECTOR_N = TETRAHEDRON_DOF_COUNT;

}

#endif //PROJECT_TETRAHEDRONSOLVERDATA_H
