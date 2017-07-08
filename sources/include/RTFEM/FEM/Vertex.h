#ifndef PROJECT_VERTEX_H
#define PROJECT_VERTEX_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/Vector3.h>

namespace rtfem {

class Vertex {
public:
    Vertex(UInt id, const Vector3&& cooridnates);
    ~Vertex();

    UInt id() const {return id_;}
    const Vector3& coordinates() const {return coordinates_;}

private:
    /**
     * ID is assumed to be continues integer numbers [0,N]
     */
    UInt id_;

    // TODO local/global ?
    Vector3 coordinates_;
};
}


#endif //PROJECT_VERTEX_H
