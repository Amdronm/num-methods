cmake_minimum_required(VERSION 3.24.0)
project(gauss_ldlt VERSION 0.1.0 LANGUAGES C CXX)

add_compile_options(-Wall -Wextra -Wpedantic -Werror)
add_compile_options(-O2)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(gauss_ldlt main.cpp)

