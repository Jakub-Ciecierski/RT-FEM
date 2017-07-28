#include "RTFEM/FEM/Vertex.h"

namespace rtfem {

Vertex::Vertex(UInt id, const Eigen::Vector3<Float>& coordinates)
        : id_(id), coordinates_(coordinates){}

}