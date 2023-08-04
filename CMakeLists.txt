cmake_minimum_required(VERSION 2.8.3)
project(quadrotor_ukf)

# set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

add_definitions(-Wall -std=c++11 -Werror=return-type)

IF(${TARGET_ARCH} MATCHES "armv7")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpu=neon -march=armv7-a")
ENDIF()
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
ELSE()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF()

find_package(catkin REQUIRED COMPONENTS roscpp sensor_msgs cmake_modules geometry_msgs tf2_ros)
find_package(Eigen REQUIRED)

include_directories(include ${catkin_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIR})

catkin_package(CATKIN_DEPENDS roscpp sensor_msgs geometry_msgs tf2_ros)

set(UKF_SOURCE_FILES     
		src/quadrotor_ukf.cpp
		src/vio_utils.cpp
)

set(UKF_HEADER_FILES  
		include/quadrotor_ukf/quadrotor_ukf.h
		include/quadrotor_ukf/vio_utils.h 
)


###########
## Build ##
###########
add_library(ukf_filter_lib ${UKF_SOURCE_FILES})
target_link_libraries(ukf_filter_lib ${catkin_LIBRARIES} ${EIGEN_LIBS})
add_executable(quadrotor_ukf  src/quadrotor_ukf_node.cpp)
add_dependencies(quadrotor_ukf ${catkin_EXPORTED_TARGETS} )
target_link_libraries(quadrotor_ukf ukf_filter_lib ${catkin_LIBRARIES} ${EIGEN_LIBS})

add_executable(quadrotor_ukf_redo  src/quadrotor_ukf_node_redo.cpp)
add_dependencies(quadrotor_ukf_redo ${catkin_EXPORTED_TARGETS} )
target_link_libraries(quadrotor_ukf_redo ukf_filter_lib ${catkin_LIBRARIES} ${EIGEN_LIBS})


install(TARGETS ukf_filter_lib quadrotor_ukf
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

