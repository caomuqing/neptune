cmake_minimum_required(VERSION 2.8.3)
project(behavior_selector)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS
	rospy
	message_generation
	std_msgs
    rqt_gui
    rqt_gui_py
)

## Generate services in the 'srv' folder
add_service_files(
  FILES
  MissionModeChange.srv
)

generate_messages(
	DEPENDENCIES
	std_msgs
)

catkin_package()


## Your package locations should be listed before other locations
include_directories(
 ${catkin_INCLUDE_DIRS}
)

install(PROGRAMS scripts/rqt_pkg
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(DIRECTORY
  resource
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
)

install(FILES
  plugin.xml
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
)