add_library(methods_runner INTERFACE)
target_include_directories(methods_runner INTERFACE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/tests
)
target_link_libraries(methods_runner INTERFACE
    vec_ops
    plot_matrix
    simple_iteration
    jacobi
    gauss_seidel
    symmetric_gauss_seidel
    chebyshev_acceleration
    symmetric_chebyshev_acceleration
)

add_executable(symmetric_acceleration_test tests/methods/symmetric_acceleration_main.cpp)
target_include_directories(symmetric_acceleration_test PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/tests
)
target_link_libraries(symmetric_acceleration_test PRIVATE
    jacobi gauss_seidel symmetric_gauss_seidel
    symmetric_chebyshev_acceleration plot_matrix vec_ops
)


add_executable(ut_tests tests/ut_test.cpp)
target_include_directories(ut_tests PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/tests
)
target_link_libraries(ut_tests PRIVATE
    vec_ops plot_matrix csr_matrix
    simple_iteration jacobi gauss_seidel symmetric_gauss_seidel
    qr_decomposition thomas_algorithm
    chebyshev_acceleration symmetric_chebyshev_acceleration
    methods_runner
)