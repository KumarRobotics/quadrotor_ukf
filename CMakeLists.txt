cmake_minimum_required(VERSION 2.8.3)
project(quadrotor_ukf)
add_definitions(-Wall -std=c++11 -Werror=return-type ${PCL_DEFINITIONS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
#set(CMAKE_BUILD_TYPE Release)

# Set build flags, set ARM_ARCHITECTURE environment variable on Odroid
SET(CMAKE_CXX_FLAGS "-Wall -D_LINUX -D_REENTRANT -march=native -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unknown-pragmas")

IF( ${TARGET_ARCH} MATCHES "armv7")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpu=neon -march=armv7-a")
ELSE()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmmx -msse -msse -msse2 -msse3 -mssse3")
ENDIF()
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
ELSE()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF()
SET(CMAKE_CXX_FLAGS_R)

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
target_link_libraries(quadrotor_ukf ukf_filter_lib ${catkin_LIBRARIES})


 install(PROGRAMS
  scripts/run_bag
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
