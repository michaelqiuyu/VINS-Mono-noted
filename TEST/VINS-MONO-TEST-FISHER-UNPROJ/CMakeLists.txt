cmake_minimum_required(VERSION 3.17)
project(VINS-MONO-TEST-FISHER-UNPROJ)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenCV 3.4.5)
include_directories(${OpenCV_INCLUDE_DIR})

find_package(Eigen3 REQUIRED)

add_executable(unproj main.cpp)
target_link_libraries(unproj ${OpenCV_LIBS})