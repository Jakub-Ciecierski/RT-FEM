cmake_minimum_required(VERSION 3.3)

set(LIB_NAME "rtfem")
project(${LIB_NAME})

#---------------------------------
# MACROS
#---------------------------------

MACRO(FIND_INCLUDE_DIR return_list regex)
    FILE(GLOB_RECURSE new_list ${regex})
    SET(dir_list "")
    FOREACH (file_path ${new_list})
        string(REGEX MATCH ^.*/include dir ${file_path})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH ()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

#---------------------------------
# COMPILE FLAGS
#---------------------------------

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
set(CMAKE_NVCC_FLAGS "${CMAKE_NVCC_FLAGS} -std=c++11 -arch=compute_61 -code=sm_61 -Wall")

#---------------------------------
# DEPENDENCIES
#---------------------------------

include_directories(../external/tetgen1.5.0)
include_directories(../external/eigen/)

# CMAKE should locate cuda using CUDA_BIN_PATH by default, but it doesn't :(
set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_BIN_PATH})
find_package(CUDA REQUIRED)

#---------------------------------
# PROJECT
#---------------------------------

set(INC_DIR include)
set(SRC_DIR src)

# INCLUDE AUTOMATIC SEARCH
FIND_INCLUDE_DIR(INC_DIRS ${INC_DIR}/*.h)
include_directories(${INC_DIR})

# SOURCES AUTOMATIC SEARCH
file(GLOB_RECURSE SRC_FILES ${SRC_DIR}/*.cpp)
file(GLOB_RECURSE INC_FILES ${INC_DIR}/*.h)

file(GLOB_RECURSE CUDA_SRC_FILES ${SRC_DIR}/*.cu)
file(GLOB_RECURSE CUDA_INC_FILES ${INC_DIR}/*.cuh)

#---------------------------------

if (NOT TARGET ${LIB_NAME})
    CUDA_ADD_LIBRARY(${LIB_NAME} STATIC
            ${SRC_FILES}
            ${INC_FILES}
            ${CUDA_SRC_FILES}
            ${CUDA_INC_FILES}
            OPTIONS -arch=compute_61 -code=sm_61)
endif ()

set_target_properties(${LIB_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY lib)
