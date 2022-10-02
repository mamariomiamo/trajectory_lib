cmake_minimum_required(VERSION 3.10)
project(trajectory_lib)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-std=c++11")

add_subdirectory(
    include/third-party/matplotlib-cpp
)

add_subdirectory(
    include/third-party/osqp
)

find_package(Eigen3 REQUIRED)

add_executable(matplotlib-minimal src/matplotlib-minimal.cpp)
target_include_directories(matplotlib-minimal PUBLIC
    #${PYTHON_INCLUDE_DIRS}
    include/third-party/matplotlib-cpp)
target_link_libraries(matplotlib-minimal
    # ${PYTHON_LIBRARIES}
    matplotlib_cpp)

add_executable(main src/main.cpp)
target_include_directories(main PUBLIC
    include/third-party
    include/third-party/matplotlib-cpp
    ${EIGEN_INCLUDE_DIRS})
target_link_libraries(main
    matplotlib_cpp
    Eigen3::Eigen
    osqp
    )

    add_executable(main_back src/main_back.cpp)
    target_include_directories(main_back PUBLIC
        include/third-party
        include/third-party/matplotlib-cpp
        ${EIGEN_INCLUDE_DIRS})
    target_link_libraries(main_back
        matplotlib_cpp
        Eigen3::Eigen
        osqp
        )