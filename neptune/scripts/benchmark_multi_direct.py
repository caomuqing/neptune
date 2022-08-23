#!/usr/bin/env python

# /* ----------------------------------------------------------------------------
#  * Copyright 2022, Cao Muqing
#  * Nanyang Technological University
#  * All Rights Reserved
#  * Authors: Cao Muqing, et al.
#  * -------------------------------------------------------------------------- */

import rospy
from mader_msgs.msg import Mode
from snapstack_msgs.msg import Goal, State
from geometry_msgs.msg import Pose, PoseStamped
from snapstack_msgs.msg import QuadFlightMode
from nav_msgs.msg import Odometry
from std_msgs.msg import Float32MultiArray

import math
import random
import numpy as np
from numpy import linalg as LA
import sys

number_of_robots = 12;
x_min = -10.0;
x_max = 10.0;
y_min = -10.0;
y_max = 10.0;
tether_length = 15; #previous was 13, this value limites the horizontal distance between base and goal
goal_height = 5.0;
close_range = 3.0;
circle_init_bases = True;
base_positions = [[-10, 5, 0.0],  #unity sim scene1
                  [-5.0, -10, 0.0],
                  [10.0, -5.0, 0.0],
                  [5.0, 10.0, 0.0],
                  [-9.5, -5.0, 0.0],
                  [5.0, -9.5, 0.0],
                  [9.5, 5.0, 0.0],
                  [-5.0, 9.5, 0.0]];  
odom_topic_name = 'unity/odom';
dc = 0.05;
quad_name = "firefly"; #firefly

class auto_commands:

    def __init__(self):
        radius = 10.0;
        one_slice = 2*math.pi/number_of_robots;
        init_pos_as_circle = False;

        self.initialized=False;
        self.start_pub_time = rospy.Time.now();
        self.get_goal_time = [rospy.Time.now() for i in range(number_of_robots)];
        self.complete_time = [100.0 for i in range(number_of_robots)];
        self.gotten_goal = [0 for i in range(number_of_robots)];

        self.dist_travelled = [0.0 for i in range(number_of_robots)];
        self.completed_current = [0 for i in range(number_of_robots)];
        self.gotten_odom = [0 for i in range(number_of_robots)];

        self.pos= [[0.0, 0.0, 0.0] for i in range(number_of_robots)];

        self.pos_prev=[[0.0, 0.0, 0.0] for i in range(number_of_robots)];
        self.goals=[[0.0, 0.0, 0.0] for i in range(number_of_robots)];   
        self.jerkSum = [0.0 for i in range(number_of_robots)];

    


        # random.seed(10);
        self.pose = Pose();
        self.round =1
        self.currentrun = 0
        self.pubRecord = rospy.Publisher('/firefly/record',PoseStamped,queue_size=3,latch=True)
        self.pubGoal_multi = rospy.Publisher('/start_and_goals',Float32MultiArray,queue_size=1,latch=True)

    #In rospy, the callbacks are all of them in separate threads
    def stateCB(self, data, args):
        idx = args
        self.pos[idx][0] = data.pos.x
        self.pos[idx][1] = data.pos.y
        self.pos[idx][2] = data.pos.z
        self.pose.orientation = data.quat
        if np.prod(self.gotten_odom) !=1:
            print("INIT: "+str(idx)+" robot has x = "+str(self.pos[idx][0]) + " y = "+str(self.pos[idx][1]))
        self.checkAndPublish(idx)
        self.pos_prev[idx][0] = data.pos.x
        self.pos_prev[idx][1] = data.pos.y
        self.pos_prev[idx][2] = data.pos.z
        self.gotten_odom[idx] = 1;

    def goalOneCB(self, data, args):
        idx = args
        self.goals[idx][0] = data.pose.position.x;
        self.goals[idx][1] = data.pose.position.y;
        self.goals[idx][2] = data.pose.position.z;         
        self.get_goal_time[idx] = rospy.Time.now();   
        self.gotten_goal[idx] = 1;                    
        # self.publishgoals();

    #In rospy, the callbacks are all of them in separate threads
    def odomCB(self, data, args):
        idx = args        
        self.pos[idx][0] = data.pose.pose.position.x
        self.pos[idx][1] = data.pose.pose.position.y
        self.pos[idx][2] = data.pose.pose.position.z
        self.pose.orientation = data.pose.pose.orientation
        if np.prod(self.gotten_odom) !=1:
            print("INIT: "+str(idx)+" robot has x = "+str(self.pos[idx][0]) + " y = "+str(self.pos[idx][1]))
        self.checkAndPublish(idx)
        self.pos_prev[idx][0] = data.pose.pose.position.x
        self.pos_prev[idx][1] = data.pose.pose.position.y
        self.pos_prev[idx][2] = data.pose.pose.position.z
        self.gotten_odom[idx] = 1;

    def cmdCB(self, data, args):
        idx = args        
        if self.completed_current[idx] != 1 and self.gotten_goal[idx]:
            self.jerkSum[idx] += LA.norm(np.array([data.j.x, data.j.y, data.j.z]))*0.05;

    def checkAndPublish(self, idx):
        if self.completed_current[idx] != 1 and self.gotten_goal[idx]:
            self.dist_travelled[idx] += \
                    math.sqrt(sum((self.pos[idx][i]-self.pos_prev[idx][i])**2 for i in range(0,3)));
        
        if not self.completed_current[idx] and self.gotten_goal[idx] and \
            LA.norm(np.array(self.pos[idx][0:3])-np.array(self.goals[idx][0:3]))<0.30:
            self.completed_current[idx] = 1;
            self.complete_time[idx] = (rospy.Time.now()-self.get_goal_time[idx]).to_sec();
            print("Robot "+str(idx)+" has completed current run !")
            if np.prod(self.completed_current) ==1:
                rospy.logwarn("Robots have completed the journey !")
                rospy.logwarn("The longest time taken is %f", max(self.complete_time));
                rospy.logwarn("The average distance taken is %f", sum(self.dist_travelled)/number_of_robots);
                rospy.logwarn("The average jerk is %f", sum(self.jerkSum)/number_of_robots);
                rospy.signal_shutdown("shutting down")
                sys.exit(1) 
                # print("The "+str(self.pos[idx][0]) + " y = "+str(self.pos[idx][1]))

    def setupAndPublishRecord(self, success):
        msgRecord=PoseStamped()  
        msgRecord.header.frame_id="world"
        msgRecord.header.stamp = rospy.get_rostime()
        msgRecord.header.seq = self.currentrun-1;
        msgRecord.pose.position.x = (rospy.Time.now()-self.start_pub_time).to_sec();
        msgRecord.pose.position.y = sum(self.dist_travelled)/number_of_robots;
        msgRecord.pose.position.z = success;            
        self.pubRecord.publish(msgRecord);




                  
def startNode():
    c = auto_commands()
    #s = rospy.Service("/change_mode",MissionModeChange,c.srvCB)
    for i in range(0,number_of_robots):
        rospy.Subscriber('/'+quad_name+str(i+1)+"/state", State, c.stateCB, i)
        rospy.Subscriber('/'+quad_name+str(i+1)+"/term_goal", PoseStamped, c.goalOneCB, i)
        rospy.Subscriber('/'+quad_name+str(i+1)+"/"+odom_topic_name, Odometry, c.odomCB, i)
        rospy.Subscriber('/'+quad_name+str(i+1)+"/goal", Goal, c.cmdCB, i)

    rospy.spin()

if __name__ == '__main__':
    rospy.init_node('auto_commands')  
    startNode()
    print "auto_commands started" 