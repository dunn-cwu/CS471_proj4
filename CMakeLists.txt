cmake_minimum_required(VERSION 3.1.0)

# Use C++ 11 standard
set(CMAKE_CXX_STANDARD 11)

project(cs471_proj1)

#Specify include dir
include_directories(include)

# Add all *.cpp files in ../src
file(GLOB SOURCES "src/*.cpp")

add_executable(cs471_proj1.out ${SOURCES})

