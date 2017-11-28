#include "RTFEM/GPU/DeviceProperties/DeviceProperties.cuh"

namespace rtfem {

DeviceProperties::DeviceProperties(const cudaDeviceProp& properties) :
    properties_(properties){}

std::ostream& operator<<(std::ostream& os,
                         const DeviceProperties& device_properties){
    os << "GPU Device Properties" << std::endl;
    os << "Name: " << device_properties.properties_.name << std::endl;

    double global_memory_giga_bytes =
        device_properties.properties_.totalGlobalMem / 1000000000.0;
    os << "Total Global Memory [GB]: "
       << global_memory_giga_bytes
       << std::endl;

    return os;
}

}