cmake_minimum_required(VERSION 3.10)

# Projektname
project(fib)

# C++-Standard setzen
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add includes
include_directories(include)
# Quellverzeichnis angeben
add_executable(MyProject src/main.cpp)
