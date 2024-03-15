cmake_minimum_required(VERSION 3.10)
project(casadi_mpc)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-std=c++14)

find_package(casadi REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Python3 COMPONENTS Development NumPy)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
)

add_library(casadi_mpc INTERFACE)
target_include_directories(casadi_mpc INTERFACE include)

catkin_package(
 INCLUDE_DIRS include
  LIBRARIES casadi_mpc
 CATKIN_DEPENDS roscpp std_msgs
 DEPENDS system_lib
)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
)

add_executable(inverted_pendulum_mpc src/inverted_pendulum_mpc.cpp)
target_link_libraries(inverted_pendulum_mpc ${catkin_LIBRARIES} casadi_mpc casadi Python3::Python Python3::NumPy)
add_dependencies(inverted_pendulum_mpc ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

add_executable(cartpole_nmpc src/cartpole_nmpc.cpp)
target_link_libraries(cartpole_nmpc ${catkin_LIBRARIES} casadi_mpc casadi Python3::Python Python3::NumPy)
add_dependencies(cartpole_nmpc ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})