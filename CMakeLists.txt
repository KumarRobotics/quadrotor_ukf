cmake_minimum_required(VERSION 2.8.3)
project(quadrotor_ukf)
add_definitions(-Wall -std=c++11 -Werror=return-type ${PCL_DEFINITIONS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
#set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-msse2")


find_package(cmake_modules REQUIRED)
find_package(catkin REQUIRED COMPONENTS roscpp sensor_msgs)
find_package(Eigen REQUIRED)


catkin_package(CATKIN_DEPENDS roscpp sensor_msgs)

set(UKF_SOURCE_FILES     
		src/quadrotor_ukf.cpp
		src/vio_utils.cpp

)


set(UKF_HEADER_FILES  
		include/quadrotor_ukf/quadrotor_ukf.h
		include/quadrotor_ukf/vio_utils.h 

)

include_directories(
  include
  ${EIGEN_INCLUDE_DIR}
#  LIBRARIES imu_3dm_gx3
#  CATKIN_DEPENDS other_catkin_pkg
#  DEPENDS system_lib
)

###########
## Build ##
###########
add_library(ukf_filter_lib ${UKF_SOURCE_FILES})
add_executable(quadrotor_ukf  src/quadrotor_ukf_node.cpp)
 set_property(TARGET quadrotor_ukf APPEND_STRING PROPERTY COMPILE_FLAGS "-O3 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsee -funsafe-loop-optimizations -fno-signed-zeros -fno-math-errno -funroll-loops")
target_link_libraries(quadrotor_ukf ukf_filter_lib ${catkin_LIBRARIES})


 install(PROGRAMS
  scripts/run_bag
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
