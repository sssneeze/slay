# Vector operations
add_library(vec_ops src/vector/operations_vector.cpp)
target_include_directories(vec_ops PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Matrices
add_library(plot_matrix src/matrix/plot_matrix.cpp)
target_include_directories(plot_matrix PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(plot_matrix PUBLIC vec_ops)

add_library(csr_matrix src/matrix/csr_matrix.cpp)
target_include_directories(csr_matrix PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(csr_matrix PUBLIC vec_ops)

# Basic iterative methods
add_library(simple_iteration src/methods/simple_iteration.cpp)
target_include_directories(simple_iteration PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(simple_iteration PUBLIC plot_matrix vec_ops)

add_library(jacobi src/methods/jacobi.cpp)
target_include_directories(jacobi PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(jacobi PUBLIC plot_matrix vec_ops)

add_library(gauss_seidel src/methods/gauss_seidel.cpp)
target_include_directories(gauss_seidel PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(gauss_seidel PUBLIC plot_matrix vec_ops)

# Symmetric methods & Chebyshev
add_library(symmetric_gauss_seidel src/methods/symmetric_gauss_seidel.cpp)
target_include_directories(symmetric_gauss_seidel PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(symmetric_gauss_seidel PUBLIC plot_matrix vec_ops)

add_library(chebyshev_acceleration src/methods/chebyshev_acceleration.cpp)
target_include_directories(chebyshev_acceleration PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(chebyshev_acceleration PUBLIC plot_matrix vec_ops)

add_library(symmetric_chebyshev_acceleration src/methods/symmetric_chebyshev_acceleration.cpp)
target_include_directories(symmetric_chebyshev_acceleration PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(symmetric_chebyshev_acceleration PUBLIC plot_matrix vec_ops)

# Direct methods
add_library(qr_decomposition src/methods/householder_alg.cpp)
target_include_directories(qr_decomposition PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(qr_decomposition PUBLIC plot_matrix vec_ops)

add_library(thomas_algorithm src/methods/progonka.cpp)
target_include_directories(thomas_algorithm PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(thomas_algorithm PUBLIC vec_ops)