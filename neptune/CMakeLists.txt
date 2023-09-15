cmake_minimum_required(VERSION 2.8.3)
project(neptune)

## Add support for C++11, supported in ROS Kinetic and newer
add_definitions(-std=c++11)
#add_definitions(-std=c99)

find_package(catkin REQUIRED COMPONENTS 
	roscpp
	rospy
	snapstack_msgs
	mader_msgs
	sensor_msgs
	decomp_util
	decomp_ros_utils
	rviz_visual_tools
	separator
)
find_package( Eigen3 REQUIRED )
find_package(decomp_util REQUIRED)
find_package(CGAL REQUIRED )

include_directories(${EIGEN3_INCLUDE_DIR}
	${catkin_INCLUDE_DIRS}
	include
	third_party)
include_directories(${DECOMP_UTIL_INCLUDE_DIRS})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_finders")




option(USE_GUROBI "Use GUROBI as the solver" ON) #if OFF, then it will use NLOPT as the solver


include(${CGAL_USE_FILE})

set(CMAKE_BUILD_TYPE Release) #From terminal: catkin config -DCMAKE_BUILD_TYPE=Release (or RelWithDebInfo) or Debug

catkin_package()


add_library(cgal_utils STATIC src/cgal_utils)
target_include_directories (cgal_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cgal_utils PUBLIC ${DECOMP_UTIL_LIBRARIES} ${catkin_LIBRARIES})

add_library(entangle_utils STATIC src/entangle_utils)
target_include_directories (entangle_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(entangle_utils PUBLIC ${DECOMP_UTIL_LIBRARIES} ${catkin_LIBRARIES})

add_library(utils STATIC src/utils.cpp)
target_include_directories (utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(utils PUBLIC ${catkin_LIBRARIES})


add_library(gjk STATIC src/gjk.cpp) 
target_include_directories (gjk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gjk PUBLIC ${catkin_LIBRARIES}) 

add_library(kinodynamic_search STATIC src/kinodynamic_search.cpp src/bspline_utils.cpp src/utils.cpp src/entangle_utils.cpp) 
target_include_directories (kinodynamic_search PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(kinodynamic_search PUBLIC ${catkin_LIBRARIES} gjk) 



find_package(GUROBI REQUIRED)

if(GUROBI_FOUND)
	message(STATUS "GUROBI FOUND")
else(GUROBI_FOUND)
	message(FATAL_ERROR "GUROBI NOT FOUND!")
endif(GUROBI_FOUND)

include_directories(${GUROBI_INCLUDE_DIRS})
file(GLOB GurobiSOFiles $ENV{GUROBI_HOME}/lib/libgurobi*[0-9].so) #files that are start with libgurobi and end with number.so
set(GUROBI_LIBRARIES "$ENV{GUROBI_HOME}/lib/libgurobi_c++.a;${GurobiSOFiles};$ENV{GUROBI_HOME}/lib/" )

add_library(poly_solver STATIC src/solver_gurobi_poly.cpp)
target_compile_definitions(poly_solver PUBLIC USE_GUROBI_FLAG=1) 
target_include_directories (poly_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${DECOMP_UTIL_INCLUDE_DIRS})
target_link_libraries(poly_solver PUBLIC ${GUROBI_LIBRARIES} ${DECOMP_UTIL_LIBRARIES} utils ${catkin_LIBRARIES} cgal_utils)




add_executable(neptune_test_node src/node.cpp src/neptune_ros.cpp src/neptune.cpp src/Astar_searcher.cpp 
			   src/solveMultiLinearPath.cpp src/bezier_base.cpp)
target_link_libraries(neptune_test_node PUBLIC ${catkin_LIBRARIES} kinodynamic_search poly_solver utils gjk)  #${CGAL_LIBS}
add_dependencies(neptune_test_node ${catkin_EXPORTED_TARGETS} )


add_executable(mtlp_node src/mtlp_node.cpp)
add_dependencies(mtlp_node ${catkin_EXPORTED_TARGETS})
target_link_libraries(mtlp_node ${catkin_LIBRARIES} utils cgal_utils)




