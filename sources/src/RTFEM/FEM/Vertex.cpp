#include "RTFEM/FEM/Vertex.h"

rtfem::Vertex::Vertex(rtfem::UInt id,
                      const rtfem::Vector3 &&cooridnates) :
    id_(id), coordinates_(cooridnates){}

rtfem::Vertex::~Vertex() {

}
