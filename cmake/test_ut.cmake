
add_library(methods_runner INTERFACE)
target_include_directories(methods_runner INTERFACE 
    ${CMAKE_SOURCE_DIR}/tests/utils
    ${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(methods_runner INTERFACE
    simple_iteration
    jacobi
    gauss_seidel
    plot_matrix
    vec_ops
)
