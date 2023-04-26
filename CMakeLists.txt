cmake_minimum_required(VERSION 3.0.2)
project(liosam_2gridmap)

# add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  pcl_ros
  std_msgs
  lio_sam
)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES liosam_2gridmap
 CATKIN_DEPENDS roscpp rospy std_msgs
#  DEPENDS system_lib
)

###########
## Build ##
###########


include_directories(
# include
  ${catkin_INCLUDE_DIRS}
)


add_executable(${PROJECT_NAME}_node src/liosam_2gridmap_node.cpp)

target_link_libraries(${PROJECT_NAME}_node
  ${catkin_LIBRARIES}
)


