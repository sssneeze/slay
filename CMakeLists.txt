cmake_minimum_required(VERSION 3.10)
project(SLAY LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -O0)
    else()
        add_compile_options(-O2)
    endif()
endif()


include(${CMAKE_SOURCE_DIR}/cmake/src_include.cmake)


include(${CMAKE_SOURCE_DIR}/cmake/tests_utils.cmake)


add_executable(ut_tests tests/ut_test.cpp)
target_include_directories(ut_tests PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/tests/utils
)
target_link_libraries(ut_tests PUBLIC
    vec_ops
    plot_matrix
    csr_matrix
    simple_iteration
    jacobi
    gauss_seidel
    qr_decomposition
    thomas_algorithm
    methods_runner
)


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")