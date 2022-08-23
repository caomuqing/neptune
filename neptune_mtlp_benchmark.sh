#!/bin/bash
# Author: Cao Muqing

#This file should be called from the root of the workspace 

#Example: If the structure is
# ws --> src --> neptune
# This file should be called from the directory ws

path_to_ws=$(pwd)


xterm -e "source ~/.bashrc && source $path_to_ws/devel/setup.bash && roscore" & sleep 2
xterm -hold -e "source ~/.bashrc && source $path_to_ws/devel/setup.bash && rosrun rviz rviz -d $path_to_ws/src/neptune/neptune/rviz_cfgs/benchmark_mtlp.rviz" & sleep 3



xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune many_drones_mq.launch action:=start" & sleep 3

xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune neptune_mtlp_benchmark.launch" & sleep 3
xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune kino_test.launch quad:=firefly1" & sleep 1
xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune kino_test.launch quad:=firefly2" & sleep 1
xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune kino_test.launch quad:=firefly3" & sleep 1
xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune kino_test.launch quad:=firefly4" & sleep 1

xterm -hold -e "source $path_to_ws/devel/setup.bash && roslaunch neptune kino_test.launch quad:=firefly5" & sleep 1

sleep 3
#xterm -e "source ~/.bashrc && source $path_to_ws/devel/setup.bash && cd /home/iot/catkin_ws_neptune/bags/ && rosbag record /firefly/record /mtlp_log /firefly/log_for_plot" & sleep 3

xterm -hold -e "source $path_to_ws/devel/setup.bash && rosrun neptune benchmark_mtlp.py" & sleep 3


