cmake_minimum_required(VERSION 3.3)

set(APP_NAME "rtfem_unittests")
project(${APP_NAME})

set(INC_DIR include)
set(SRC_DIR src)
set(DEPS_DIR dependencies)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

#---------------------------------
# MACROS
#---------------------------------

MACRO(FIND_INCLUDE_DIR return_list regex)
    FILE(GLOB_RECURSE new_list ${regex})
    SET(dir_list "")
    FOREACH(file_path ${new_list})
        string(REGEX MATCH ^.*/include dir ${file_path})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

# INCLUDE AUTOMATIC SEARCH
FIND_INCLUDE_DIR(INC_DIRS ./include/*.h)
include_directories(${INC_DIR} ${INC_DIR})

# SOURCES AUTOMATIC SEARCH
file(GLOB_RECURSE SRC_FILES src/*.cpp)
set(SOURCE_FILES )

add_executable(${APP_NAME} ${SOURCE_FILES} ${SRC_FILES} ${INC_DIR} ${INC_DIR})

#----------------------------------
# EXTERNAL
#----------------------------------

if (NOT TARGET gtest)
    add_subdirectory(external/googletest)
endif()
include_directories(external/googletest/include)

include_directories(../external/tetgen1.5.0)
include_directories(../external/eigen/)

include_directories(../sources/include)

#---------------------------------
# LINK
#---------------------------------

target_link_libraries(${APP_NAME} rtfem)
target_link_libraries(${APP_NAME} tet)
target_link_libraries(${APP_NAME} gtest gtest_main)