cmake_minimum_required(VERSION 2.8.9)

project(fluidsim)
set(CMAKE_BUILD_TYPE Release)

if (MSVC)
    set(CMAKE_CXX_FLAGS "/O2 /W4 /EHsc")
else()
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -Wall")
endif()

include_directories(src)
file(GLOB SOURCES "src/*.cpp")

add_executable(fluidsim ${SOURCES})

file(COPY "${CMAKE_SOURCE_DIR}/sample_meshes/" DESTINATION "${CMAKE_BINARY_DIR}/sample_meshes")
