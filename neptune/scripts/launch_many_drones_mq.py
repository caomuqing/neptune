#!/usr/bin/env python
# coding=utf-8

# /* ----------------------------------------------------------------------------
#  * Copyright 2022, Cao Muqing
#  * Nanyang Technological University
#  * All Rights Reserved
#  * Authors: Cao Muqing, et al.
#  * Acknowledgement: Jesus Tordesillas
#  * See LICENSE file for the license information
#  * -------------------------------------------------------------------------- */


import math
import os
import sys
import time
from random import *
# import numpy as np
# from pyquaternion import Quaternion
from tf.transformations import quaternion_from_euler, euler_from_quaternion

num_of_agents=5; 
radius=10;
circle_start = True;
init_height = 0;

def create_session(session_name, commands):

    os.system("tmux new -d -s "+str(session_name)+" -x 300 -y 300")

    for i in range(len(commands)):
        print('splitting ',i)
        os.system('tmux split-window ; tmux select-layout tiled')
   
    for i in range(len(commands)):
        os.system('tmux send-keys -t '+str(session_name)+':0.'+str(i) +' "'+ commands[i]+'" '+' C-m')
    print("Commands sent")


def convertToStringCommand(action,quad,x,y,z,goal_x,goal_y,goal_z, yaw):
    if(action=="start"):
        return "source ~/catkin_ws_neptune/devel/setup.bash && roslaunch neptune perfect_tracker_and_sim.launch gazebo:=false quad:="+quad+" x:="+str(x)+" y:="+str(y)+" z:="+str(z)+" yaw:="+str(yaw);
    if(action=="send_goal"):
        return "source ~/catkin_ws_neptune/devel/setup.bash && rostopic pub /"+quad+"/term_goal geometry_msgs/PoseStamped '{header: {stamp: now, frame_id: 'world'}, pose: {position: {x: "+str(goal_x)+", y: "+str(goal_y)+", z: "+str(goal_z)+"}, orientation: {x: 0.0, y: 0.0, z: 0.0, w: 0.0}}}'"
    if(action=="neptune"):
        if(quad=="SQ01s"):
            return "source ~/catkin_ws_neptune/devel/setup.bash && roslaunch neptune kino_test.launch quad:="+quad
        else:
            return "source ~/catkin_ws_neptune/devel/setup.bash && roslaunch neptune kino_test.launch quad:="+quad
        #+ " >> "+quad+".txt"
        # return "script -q -c 'roslaunch neptune neptune.launch quad:="+quad + "' "+quad+".txt"
        

if __name__ == '__main__':
    # formation="sphere", "square" "circle"
    formation="sphere"
    commands = []


    if(formation=="sphere"):
        # num_mer=int(math.sqrt(num_of_agents)); #Num of meridians
        # num_of_agents_per_mer=int(math.sqrt(num_of_agents));    #Num of agents per meridian
        if(num_of_agents%3==0):
            num_mer=max(int(num_of_agents/4.0),3); #Num of meridians
        else: #either divisible by 4 or will force it by changing num of agents
            num_mer=max(int(num_of_agents/4.0),4); #Num of meridians
        num_of_agents_per_mer=int(num_of_agents/num_mer);    #Num of agents per meridian

    if(formation=="circle" or formation=="square"):
        num_mer=num_of_agents
        num_of_agents_per_mer=1

    print("num_mer= ", num_mer)
    print("num_of_agents_per_mer= ", num_of_agents_per_mer)

    id_number=1;
    shift_z=radius;
    shift_z=1.0
    one_slice = 2*math.pi/num_of_agents;
    #TODO: Implement the square as well for other number_of_agents
    # square_starts=[[-3.0, -2.4, 1.0],   #for small scene like iot 
    #                 [3.8, -2.4, 1.0], 
    #                 [0.0, 2.5, 1.0], 
    #                 [4.0, 4.0, 1.0],
    #                 [-4.0, 0.0, 1.0],
    #                 [-4.0, -4.0, 1.0],
    #                 [0.0, -4.0, 1.0],
    #                 [4.0, -4.0, 1.0] ]

    square_starts=[[-10.0, 0.4, 1.0],   #for square start
                    [0.8, -10.4, 1.0], 
                    [10.0, -0.5, 1.0], 
                    [0.0, 10.0, 1.0],
                    [-8.0, -8.0, 1.0],
                    [8.0, -8.0, 1.0],
                    [8.0, 8.0, 1.0],
                    [-8.0, 8.0, 1.0] ]

    # square_starts=[[-24,-5, 1.0],   #for 8 drones hexadecagon
    #                 [-20,13, 1.0], 
    #                 [-5,24, 1.0], 
    #                 [13,20, 1.0],
    #                 [24,5, 1.0],
    #                 [20,-13, 1.0],
    #                 [5,-24, 1.0],
    #                 [-13,-20, 1.0] ]

    square_goals=  [[-4.0, 0.0, 1.0],
                    [-4.0, -4.0, 1.0],
                    [0.0, -4.0, 1.0],
                    [4.0, -4.0, 1.0],
                    [4.0, 0.0, 1.0],
                    [4.0, 4.0, 1.0],
                    [0.0, 4.0, 1.0],
                    [-4.0, 4.0, 1.0]];

    square_yaws_deg=  [-180.0, -135.0, -90.0, -45.0, 0.0, 45.0, 90.0, 135.0];

    if circle_start:
        for x in range(0, num_of_agents):
            square_starts[x][0] = -radius*math.cos(one_slice*x);
            square_starts[x][1] = -radius*math.sin(one_slice*x);

    for i in range(1, num_of_agents+1):
        x= square_starts[id_number-1][0]
        y=square_starts[id_number-1][1]
        z=init_height

        pitch=0.0;
        roll=0.0;
        yaw= 0#+math.pi  

        goal_x=square_starts[id_number-1][0]
        goal_y=0
        goal_z=0
            
        quad="firefly" + str(id_number);
        id_number=id_number+1;

        commands.append(convertToStringCommand(sys.argv[1],quad,x,y,z,goal_x,goal_y,goal_z, yaw));

        x_tmp="{:5.3f}".format(x);
        y_tmp="{:5.3f}".format(y);
        z_tmp="{:5.3f}".format(z);

        goal_x_tmp="{:5.3f}".format(goal_x);
        goal_y_tmp="{:5.3f}".format(goal_y);
        goal_z_tmp="{:5.3f}".format(goal_z);

        print (' "start": [',x_tmp,', ',y_tmp,', ',z_tmp,'], "goal": [',goal_x_tmp,', ',goal_y_tmp,', ',goal_z_tmp,']  ')


    print("len(commands)= " , len(commands))
    session_name=sys.argv[1] + "_session"
    os.system("tmux kill-session -t" + session_name)
    create_session(session_name, commands)
    if(sys.argv[1]!="send_goal"):
        os.system("tmux attach") #comment if you don't want to visualize all the terminals
    else: ##if send_goal, kill after some time
        time.sleep(num_of_agents); #The more agents, the more I've to wait to make sure the goal is sent correctly
        os.system("tmux kill-session -t" + session_name)

    # half_of_agents=num_of_agents/2.0
    # dist_bet_groups=6.0
    # random_01=randint(0, 1)
    # print("random_01= ", random_01)

    # positions=[];
    # thetas=[];

    # z_value=0.0;


    # for i in range(1,num_of_agents+1):
    #     theta=(2*math.pi)*i/(1.0*num_of_agents)
    #     thetas.append(theta)

    # for i in range(1,num_of_agents+1):

    #     # group = (i>half_of_agents)

    #     # x= i if group == 0 else (i-half_of_agents)
    #     # y= dist_bet_groups*group   
    #     # z=0

    #     # goal_x=half_of_agents-x
    #     # goal_y=random_01*(dist_bet_groups-y) + (1-random_01)*y 
    #     # goal_z=0

    #     if(sphere):
    #         x=radius*math.cos(theta)*math.sin(phi)
    #         y=radius*math.sin(theta)*math.sin(phi)
    #         z=radius*cos(phi)

    #     theta=thetas[i-1];

    #     x=radius*math.cos(theta)
    #     y=radius*math.sin(theta)
    #     z=1.0 #z_value
    #     #z_value=1.0 #From now on, stay on z=1 meter

    #     pitch=0.0;
    #     roll=0.0;
    #     yaw= theta+math.pi  

    #     theta=theta+math.pi

    #     goal_x=radius*math.cos(theta)
    #     goal_y=radius*math.sin(theta)
    #     goal_z=z

    #     thetas[i-1]=theta;

      
    #     # quat = quaternion_from_euler(yaw, pitch, roll, 'szyx')
    #     # print (quat)






#import libtmux
    # panes=win.list_panes()
    # print("panes.size=", len(panes))
    # for i in range(len(panes)):
    #     panes[i].send_keys(commands[i])


    # logging.info(panes)
    # logging.info(win)

            #win.select_layout(layout='tiled')
        #win.split_window()
        #win.cmd('select-layout tiled')  
        #win.select_layout(layout='tiled')
        #win.cmd('split-window', '-h')    

                #win.cmd('select-layout tiled')
        #os.system("tmux attach")
        #os.system("tmux select-layout tiled")

            #os.system("tmux attach")

    # win = session.new_window(attach=False, window_name="win")
    # win.select_layout(layout='tiled')
    #logging.info(commands)
    # pane_NUM = 3
    # WINDOW_NUM = int(math.ceil(len(commands)/4.0))  # in python3, we can use 4 also

    # server = libtmux.Server()
    # session = server.new_session(session_name)
    # panes = []
