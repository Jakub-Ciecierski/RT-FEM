#include <RTFEM/FEM/Meshing/TriangleMesh.h>

namespace rtfem {

template
struct TriangleFaceWithPoints<double>;
template
struct TriangleFaceWithPoints<float>;

template
struct TriangleMesh<double>;
template
struct TriangleMesh<float>;

template
struct TriangleMeshIndexed<double>;
template
struct TriangleMeshIndexed<float>;

template
struct TriangleFace<double>;
template
struct TriangleFace<float>;

}
