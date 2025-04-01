cmake_minimum_required(VERSION 3.31)
cmake_policy(SET CMP0135 NEW)

# Projektname
project(fib)

# C++-Standard setzen
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  DOWNLOAD_EXTRACT_TIMESTAMP true
)
FetchContent_MakeAvailable(googletest)

# add includes
include_directories(include)
# Quellverzeichnis angeben
add_executable(MyProject src/main.cpp)

enable_testing()

add_executable(
  test_stuff
  src/test_stuff.cc
)

target_link_libraries(
  test_stuff
  GTest::gtest_main
)

include(GoogleTest)

set(CMAKE_CXX_COMPILER /usr/bin/gcc)
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_EXE_LINKER_FLAGS "-lstdc++")

target_link_libraries(MyProject m)
target_link_libraries(test_stuff m)


gtest_discover_tests(test_stuff)

