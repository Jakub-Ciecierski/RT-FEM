#include <RTFEM/FEM/Meshing/TriangleMesh.h>

namespace rtfem {

template
struct Triangle<double>;
template
struct Triangle<float>;

template
struct TriangleMesh<double>;
template
struct TriangleMesh<float>;

template
struct TriangleMeshIndexed<double>;
template
struct TriangleMeshIndexed<float>;

}
