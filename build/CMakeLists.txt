cmake_minimum_required(VERSION 2.8.9)
project(cs471_proj1)

#Specify include dir
include_directories(../include)

# Add all *.cpp files in ../src
file(GLOB SOURCES "../src/*.cpp")

add_executable(cs471_proj1.out ${SOURCES})

