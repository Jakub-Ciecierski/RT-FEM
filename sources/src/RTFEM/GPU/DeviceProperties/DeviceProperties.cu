#include "RTFEM/GPU/DeviceProperties/DeviceProperties.cuh"

#include <assert.h>

namespace rtfem {

DeviceProperties::DeviceProperties(const cudaDeviceProp& properties) :
    properties_(properties){}

double DeviceProperties::GetUsedMegaBytesCount(){
    size_t free_byte ;

    size_t total_byte ;

    auto cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
    assert(cudaSuccess == cuda_status);

    double used = (double)total_byte - (double)free_byte;

    return used / 1024.0 / 1024.0;
}

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