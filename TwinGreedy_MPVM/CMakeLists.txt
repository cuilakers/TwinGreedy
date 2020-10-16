cmake_minimum_required(VERSION 3.15)
project(TwinGreedy)

set(CMAKE_CXX_STANDARD 14)
include_directories ("${PROJECT_SOURCE_DIR}/sfmt")


add_executable(TwinGreedy main.cc allocator.h utils.cc TimGraph.cc anyoption.cc sfmt/SFMT.c allocator.cc)
