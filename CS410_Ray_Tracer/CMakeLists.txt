cmake_minimum_required(VERSION 3.17)
project(P5)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")

include_directories(${CMAKE_SOURCE_DIR}/eigen-3.3.8)

add_executable(raytracer test.cpp Camera.h Camera.cpp Light.cpp Light.h Sphere.cpp Sphere.h Ray.cpp Ray.h Triangle.cpp Triangle.h Obj.cpp Obj.h)