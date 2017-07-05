cmake_minimum_required(VERSION 3.3)

SET(BUILD_UNITTESTS ON)

add_subdirectory(sources)

if (${BUILD_UNITTESTS})
    add_subdirectory(tests)
endif()

