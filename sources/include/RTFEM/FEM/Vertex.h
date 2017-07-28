#ifndef PROJECT_VERTEX_H
#define PROJECT_VERTEX_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/Vector3.h>

#include <Eigen/Core>

namespace rtfem {

class Vertex {
public:
    Vertex(UInt id, const Eigen::Vector3<Float>& coordinates);

    ~Vertex() = default;

    UInt id() const {return id_;}
    const Eigen::Vector3<Float>& coordinates() const {return coordinates_;}
    Float x() const {return coordinates_(0);}
    Float y() const {return coordinates_(1);}
    Float z() const {return coordinates_(2);}

private:
    /**
     * ID is assumed to be continues integer numbers [0,N]
     */
    UInt id_;

    // TODO local/global ?
    //Vector3 coordinates_;
    Eigen::Matrix<Float, 3, 1> coordinates_;
};
}


#endif //PROJECT_VERTEX_H
