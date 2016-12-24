cmake_minimum_required(VERSION 2.8.9)

project(fluidsim)
set(CMAKE_BUILD_TYPE Release)

if (MSVC)
    add_definitions(-DCL_USE_DEPRECATED_OPENCL_1_1_APIS -DNOMINMAX -D_CRT_SECURE_NO_WARNINGS)
endif()

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /W4 /EHsc /wd4244")
elseif(MINGW)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall")
elseif(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -fPIC")
endif()

include_directories(src)
file(GLOB SOURCES "src/*.cpp")

add_library(objects OBJECT ${SOURCES})
add_executable(fluidsim $<TARGET_OBJECTS:objects>)
