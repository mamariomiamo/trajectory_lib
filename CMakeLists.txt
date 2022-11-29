cmake_minimum_required(VERSION 3.10)
project(trajectory_lib)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-std=c++17")

add_subdirectory(
    include/third-party/osqp
)

# add_subdirectory(
# include/third-party/eigen-3.4.0
# )
find_package(Eigen3 REQUIRED)

add_library(min_snap_traj_gen src/min_snap.cpp)

target_include_directories(min_snap_traj_gen PUBLIC
    include
    include/third-party
    include/third-party/eigen-3.4.0
    ${EIGEN_INCLUDE_DIRS})
target_link_libraries(min_snap_traj_gen
    Eigen3::Eigen
    osqp
)

add_executable(main src/main.cpp)
target_link_libraries(main
    min_snap_traj_gen
)

add_executable(main_back src/main_back.cpp)
target_include_directories(main_back PUBLIC
    include/third-party
    ${EIGEN_INCLUDE_DIRS})
target_link_libraries(main_back
    Eigen3::Eigen
    osqp
)

add_executable(main_test src/main_test.cpp)
target_link_libraries(main_test
    min_snap_traj_gen
)