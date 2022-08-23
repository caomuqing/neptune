#!/bin/bash
# Author: Cao Muqing

#This file should be called from the root of the workspace 

#Example: If the structure is
# ws --> src --> neptune
# This file should be called from the directory ws

path_to_ws=$(pwd)

xterm -e "source ~/.bashrc && source $path_to_ws/devel/setup.bash && roscore" & sleep 2
xterm -e "source ~/.bashrc && source $path_to_ws/devel/setup.bash && roslaunch ros_tcp_endpoint endpoint.launch" & sleep 2

xterm -e "source ~/.bashrc && source $path_to_ws/devel/setup.bash && roslaunch tcc unity.launch" & sleep 1

xterm -e "source ~/.bashrc && source $path_to_ws/devel/setup.bash && rosrun rviz rviz -d $path_to_ws/src/neptune/neptune/rviz_cfgs/neptune.rviz" & sleep 3


xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune many_drones_mq.launch action:=start" & sleep 3

xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune kino_test.launch mode:=multi_obstacle quad:=firefly1" & sleep 1
xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune kino_test.launch mode:=multi_obstacle quad:=firefly2" & sleep 1
xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune kino_test.launch mode:=multi_obstacle quad:=firefly3" & sleep 1
xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune kino_test.launch mode:=multi_obstacle quad:=firefly4" & sleep 1

xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune kino_test.launch mode:=multi_obstacle quad:=firefly5" & sleep 1
sleep 3
xterm -hold -e "source $path_to_ws/devel/setup.bash && rosrun neptune auto_commands_together.py" & sleep 3