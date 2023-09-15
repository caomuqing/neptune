#!/bin/bash
# Author: Cao Muqing

#This file should be called from the root of the workspace 

#Example: If the structure is
# ws --> src --> neptune
# This file should be called from the directory ws

path_to_ws=$(pwd)


xterm -e "source ~/.bashrc && source $path_to_ws/devel/setup.bash && roscore" & sleep 2
xterm -e "source ~/.bashrc && source $path_to_ws/devel/setup.bash && rosrun rviz rviz -d $path_to_ws/src/neptune/neptune/rviz_cfgs/benchmark_single.rviz" & sleep 3


xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune many_drones_mq.launch action:=start" & sleep 3

xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune neptune_single_benchmark.launch quad:=firefly1" & sleep 1



