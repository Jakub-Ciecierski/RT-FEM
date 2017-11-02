#ifndef PROJECT_TETRAHEDRONSOLVERDATA_H
#define PROJECT_TETRAHEDRONSOLVERDATA_H

namespace rtfem {

/**
 * fi(x) = (Ai + Bi*x1 + Ci*x2 + Di*x3) / 6*V
 */
template<class T>
struct TetrahedronShapeFunctionCoefficients {
    // A1 = 6*V0i
    T A1;
    T A2;
    T A3;
    T A4;

    T B1;
    T B2;
    T B3;
    T B4;

    T C1;
    T C2;
    T C3;
    T C4;

    T D1;
    T D2;
    T D3;
    T D4;
};

/**
 * xij = xi - xj
 */
template<class T>
struct Edges {
    T x32 = 0;
    T x34 = 0;
    T x43 = 0;
    T x14 = 0;
    T x21 = 0;
    T x31 = 0;
    T x24 = 0;
    T x42 = 0;
    T x13 = 0;
    T x12 = 0;

    T z43 = 0;
    T z31 = 0;
    T z32 = 0;
    T z24 = 0;
    T z34 = 0;
    T z13 = 0;
    T z14 = 0;
    T z21 = 0;
    T z42 = 0;
    T z12 = 0;

    T y42 = 0;
    T y31 = 0;
    T y24 = 0;
    T y13 = 0;
    T y32 = 0;
    T y34 = 0;
    T y14 = 0;
    T y12 = 0;
    T y43 = 0;
    T y21 = 0;
};

template<class T>
struct FacesArea {
    T area1;
    T area2;
    T area3;
    T area4;
};

constexpr int TETRAHEDRON_DOF_COUNT = 12;
constexpr int TETRAHEDRON_JACOBIAN_MATRIX_N = 4;

constexpr int TETRAHEDRON_GEOMETRIC_MATRIX_N = 6;
constexpr int TETRAHEDRON_GEOMETRIC_MATRIX_M = TETRAHEDRON_DOF_COUNT;

constexpr int TETRAHEDRON_FORCE_VECTOR_N = TETRAHEDRON_DOF_COUNT;

constexpr int TETRAHEDRON_SHAPE_MATRIX_N = 3;
constexpr int TETRAHEDRON_SHAPE_MATRIX_M = TETRAHEDRON_DOF_COUNT;
}

#endif //PROJECT_TETRAHEDRONSOLVERDATA_H
