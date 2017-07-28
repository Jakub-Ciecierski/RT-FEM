cmake_minimum_required(VERSION 3.3)

SET(BUILD_UNITTESTS ON)

add_subdirectory(external)

add_subdirectory(sources)

# UnitTests
if (${BUILD_UNITTESTS})
    add_subdirectory(tests)
endif()

