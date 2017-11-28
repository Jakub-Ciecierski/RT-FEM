#ifndef PROJECT_DEVICEPROPERTIES_H
#define PROJECT_DEVICEPROPERTIES_H

#include <ostream>

#include "cuda_runtime.h"

namespace rtfem {

class DeviceProperties{
public:
    DeviceProperties(const cudaDeviceProp& properties);
    ~DeviceProperties() = default;

    const cudaDeviceProp& properties() {return properties_;}

    friend std::ostream& operator<<(std::ostream& os,
                                    const DeviceProperties& device_properties);

private:
    cudaDeviceProp properties_;
};

}

#endif //PROJECT_DEVICEPROPERTIES_H
