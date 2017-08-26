#ifndef PROJECT_TRIANGLEMESHBUILDER_H
#define PROJECT_TRIANGLEMESHBUILDER_H

#include <memory>

namespace rtfem {
template<class T>
struct TriangleMeshIndexed;
}

class TriangleMeshBuilder {
public:
    TriangleMeshBuilder() = default;

    ~TriangleMeshBuilder() = default;

    std::unique_ptr<rtfem::TriangleMeshIndexed<float>> BuildCube();
};



#endif //PROJECT_TRIANGLEMESHBUILDER_H
