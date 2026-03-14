
add_library(vec_ops src/vector/operations_vector.cpp)
target_include_directories(vec_ops PUBLIC ${CMAKE_SOURCE_DIR}/include)


add_library(plot_matrix src/matrix/plot_matrix.cpp)
target_include_directories(plot_matrix PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(plot_matrix PUBLIC vec_ops)

add_library(csr_matrix src/matrix/csr_matrix.cpp)
target_include_directories(csr_matrix PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(csr_matrix PUBLIC vec_ops)


add_library(simple_iteration src/methods/simple_iteration.cpp)
target_include_directories(simple_iteration PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(simple_iteration PUBLIC plot_matrix vec_ops)

add_library(jacobi src/methods/jacobi.cpp)
target_include_directories(jacobi PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(jacobi PUBLIC plot_matrix vec_ops)

add_library(gauss_seidel src/methods/gauss_seidel.cpp)
target_include_directories(gauss_seidel PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(gauss_seidel PUBLIC plot_matrix vec_ops)

add_library(qr_decomposition src/methods/qr_decomposition.cpp)
target_include_directories(qr_decomposition PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(qr_decomposition PUBLIC plot_matrix vec_ops)

add_library(thomas_algorithm src/methods/thomas_algorithm.cpp)
target_include_directories(thomas_algorithm PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(thomas_algorithm PUBLIC vec_ops)


add_executable(task_2 src/task_2/main.cpp)
target_include_directories(task_2 PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(task_2 PUBLIC 
    simple_iteration 
    jacobi 
    gauss_seidel
    plot_matrix
    vec_ops
)

add_executable(task_3 src/task_3/main.cpp)
target_include_directories(task_3 PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(task_3 PUBLIC
    qr_decomposition
    thomas_algorithm
    plot_matrix
    vec_ops
)