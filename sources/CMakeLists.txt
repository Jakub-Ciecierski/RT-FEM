cmake_minimum_required(VERSION 3.3)
set(LIB_NAME "rtfem")
project(${LIB_NAME})

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


#---------------------------------
# COMPILE FLAGS
#---------------------------------

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

#---------------------------------
# DEPENDENCIES
#---------------------------------

#include_directories()

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
set(SPAM)

#---------------------------------

if (NOT TARGET ${LIB_NAME})
    add_library(${LIB_NAME} STATIC ${SRC_FILES} ${INC_FILES} ${SPAM})
endif()
set_target_properties(${LIB_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY lib)
