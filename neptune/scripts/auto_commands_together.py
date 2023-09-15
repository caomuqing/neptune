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
import math
import random
import numpy as np

odom_topic_name = 'vins_estimator/odometry';
#odom_topic_name = 'unity/odom';
number_of_robots = 5;

class auto_commands:

    def __init__(self):
        radius = 10.0;
        one_slice = 2*math.pi/number_of_robots;
        init_pos_as_circle = True;

        self.initialized=False;
        self.start_pub_time = None;
        self.dist_travelled = [0.0] * 10;
        self.completed_current = [0] * number_of_robots;

        self.pos=[[4.0, 0.0, 1.0], 
                [4.0, 4.0, 1.0], 
                [0.0, 4.0, 1.0], 
                [-4.0, 4.0, 1.0],
                [-4.0, 0.0, 1.0],
                [-4.0, -4.0, 1.0],
                [0.0, -4.0, 1.0],
                [4.0, -4.0, 1.0] ]
        self.pos_prev=[[4.0, 0.0, 1.0], 
                [4.0, 4.0, 1.0], 
                [0.0, 4.0, 1.0], 
                [-4.0, 4.0, 1.0],
                [-4.0, 0.0, 1.0],
                [-4.0, -4.0, 1.0],
                [0.0, -4.0, 1.0],
                [4.0, -4.0, 1.0] ]
        # self.goals=  [[0.0, 10.0, 2.0], #for square init config
        #             [-10.0, 0.0, 2.0],
        #             [0.0, -10.0, 2.0],
        #             [10.0, 0.0, 2.0],
        #             [4.0, 0.0, 2.0],
        #             [4.0, 4.0, 2.0],
        #             [0.0, 4.0, 2.0],
        #             [-4.0, 4.0, 2.0]]  
        self.goals0=  [[-1.5, -2.0, 1.0], #for unity scene config
                    [-1.5, 2.0, 1.0],
                    [-2.0, 0.0, 1.0],
                    [0.0, 10.0, 1.0],
                    [-8.0, -8.0, 1.0],
                    [8.0, -8.0, 2.0],
                    [8.0, 8.0, 2.0],
                    [-8.0, 8.0, 2.0] ]          
        self.goals1=  [[1.5, 2.0, 1.0], #for unity scene config
                    [1.5, -2.0, 1.0],
                    [2.0, 0.0, 1.0],
                    [0.0, -10.0, 1.0],
                    [4.0, 4.0, 1.0],
                    [-4.0, 4.0, 2.0],
                    [-4.0, -4.0, 2.0],
                    [4.0, -4.0, 2.0] ]        
        # self.goals=  [[14, 0, 2.0], #for 8 drones hexadecagon
        #             [10, -10, 2.0],
        #             [0, -14, 2.0],
        #             [-10, -10, 2.0],
        #             [-14, 0, 2.0],
        #             [-10, 10, 2.0],
        #             [0, 14, 2.0],
        #             [10, 10, 2.0] ]    
        # self.goals2=  [[-24,-5, 1.0],   #for 8 drones hexadecagon
        #             [-20,13, 1.0], 
        #             [-5,24, 1.0], 
        #             [13,20, 1.0],
        #             [24,5, 1.0],
        #             [20,-13, 1.0],
        #             [5,-24, 1.0],
        #             [-13,-20, 1.0] ]                                                            
        # self.goals=  [[2.7, 0.5, 1.5], # for unity iot scene setup
        #             [-3.0, 0.0, 1.5],
        #             [0.0, -2.5, 1.5],
        #             [10.0, 0.0, 2.0],
        #             [4.0, 0.0, 2.0],
        #             [4.0, 4.0, 2.0],
        #             [0.0, 4.0, 2.0],
        #             [-4.0, 4.0, 2.0]]    
        if init_pos_as_circle:
            for x in range(0, number_of_robots):
                _radius = radius;
                if number_of_robots==8 and x%2==1:
                    _radius = _radius +0.7;
                self.goals0[x][0] = -_radius*math.cos(one_slice*x);
                self.goals0[x][1] = -_radius*math.sin(one_slice*x);
                self.goals1[x][0] = -_radius*math.cos(one_slice*x + math.pi);
                self.goals1[x][1] = -_radius*math.sin(one_slice*x + math.pi);

            # self.goals1[x][0] = -_radius*math.cos(one_slice*x + math.pi) + 0.5*math.cos(-1.571-one_slice*x);
            # self.goals1[x][1] = -_radius*math.sin(one_slice*x + math.pi) - 0.5*math.sin(-1.571-one_slice*x);          

        self.goals = [self.goals0, self.goals1]

        random.seed(10);
        self.pose = Pose();
        self.round =1
        pubGoal1 = rospy.Publisher('/firefly1/term_goal',PoseStamped,queue_size=1,latch=True)
        pubGoal2 = rospy.Publisher('/firefly2/term_goal',PoseStamped,queue_size=1,latch=True)
        pubGoal3 = rospy.Publisher('/firefly3/term_goal',PoseStamped,queue_size=1,latch=True)
        pubGoal4 = rospy.Publisher('/firefly4/term_goal',PoseStamped,queue_size=1,latch=True)
        pubGoal5 = rospy.Publisher('/firefly5/term_goal',PoseStamped,queue_size=1,latch=True)
        pubGoal6 = rospy.Publisher('/firefly6/term_goal',PoseStamped,queue_size=1,latch=True)
        pubGoal7 = rospy.Publisher('/firefly7/term_goal',PoseStamped,queue_size=1,latch=True)
        pubGoal8 = rospy.Publisher('/firefly8/term_goal',PoseStamped,queue_size=1,latch=True)
        self.currentrun = 0
        self.pubGoal = [pubGoal1, pubGoal2, pubGoal3, pubGoal4, pubGoal5, pubGoal6, pubGoal7, pubGoal8]
        self.pubRecord = rospy.Publisher('/firefly/record',PoseStamped,queue_size=3,latch=True)

    #In rospy, the callbacks are all of them in separate threads
    def stateCB(self, data, args):
        idx = args
        self.pos[idx][0] = data.pos.x
        self.pos[idx][1] = data.pos.y
        self.pos[idx][2] = data.pos.z
        self.pose.orientation = data.quat
        self.checkAndPublish(idx)
        self.pos_prev[idx][0] = data.pos.x
        self.pos_prev[idx][1] = data.pos.y
        self.pos_prev[idx][2] = data.pos.z

    #In rospy, the callbacks are all of them in separate threads
    def odomCB(self, data, args):
        idx = args        
        self.pos[idx][0] = data.pose.pose.position.x
        self.pos[idx][1] = data.pose.pose.position.y
        self.pos[idx][2] = data.pose.pose.position.z
        self.pose.orientation = data.pose.pose.orientation
        self.checkAndPublish(idx)
        self.pos_prev[idx][0] = data.pose.pose.position.x
        self.pos_prev[idx][1] = data.pose.pose.position.y
        self.pos_prev[idx][2] = data.pose.pose.position.z

    def checkAndPublish(self, idx):
        if self.completed_current[idx] != 1:
            self.dist_travelled[idx] += \
                    math.sqrt(sum((self.pos[idx][i]-self.pos_prev[idx][i])**2 for i in range(0,3)));

        if not self.initialized and idx==0:
            # xx = raw_input("type any to proceed: ");
            self.currentrun = self.currentrun+1
            self.publishgoals();
            self.initialized = True;
            print("INIT: current number of mission completed is "+str(self.currentrun-1))
        
        currentno = self.currentrun%2            
        if (self.pos[idx][0]-self.goals[currentno][idx][0])**2<0.5 and \
        (self.pos[idx][1]-self.goals[currentno][idx][1])**2<0.5:
            self.completed_current[idx] = 1;
        if np.prod(self.completed_current) == 1 and idx==0:
            # xx = raw_input("type any to proceed: ");
            self.currentrun = self.currentrun+1            
            msgRecord=PoseStamped()  
            msgRecord.header.frame_id="world"
            msgRecord.header.stamp = rospy.get_rostime()
            msgRecord.header.seq = self.currentrun-1;
            msgRecord.pose.position.x = (rospy.Time.now()-self.start_pub_time).to_sec();
            msgRecord.pose.position.y = sum(self.dist_travelled)/number_of_robots;
            self.pubRecord.publish(msgRecord);
            self.publishgoals();
            print("current number of mission completed is "+str(self.currentrun-1))

    
    def publishgoals(self):
        runno = self.currentrun%2;
        for i in range(0, number_of_robots):
            msg=PoseStamped()  
            msg.pose.orientation = self.pose.orientation
            msg.header.frame_id="world"
            msg.header.stamp = rospy.get_rostime()            
            msg.pose.position.x=self.goals[runno][i][0]
            msg.pose.position.y=self.goals[runno][i][1]
            msg.pose.position.z=self.goals[runno][i][2]   
            # print("goals1 x is %f, y is %f ", msg.pose.position.x, msg.pose.position.y)   
            if self.currentrun<3:                  
                self.pubGoal[i].publish(msg)
            self.dist_travelled[i] = 0.0;  
            self.completed_current[i] = 0;                     
        self.start_pub_time = rospy.Time.now();            



                  
def startNode():
    c = auto_commands()
    #s = rospy.Service("/change_mode",MissionModeChange,c.srvCB)
    rospy.Subscriber("/firefly1/state", State, c.stateCB, 0)
    rospy.Subscriber("/firefly2/state", State, c.stateCB, 1)
    rospy.Subscriber("/firefly3/state", State, c.stateCB, 2)
    rospy.Subscriber("/firefly4/state", State, c.stateCB, 3)
    rospy.Subscriber("/firefly5/state", State, c.stateCB, 4)
    rospy.Subscriber("/firefly6/state", State, c.stateCB, 5)
    rospy.Subscriber("/firefly7/state", State, c.stateCB, 6)
    rospy.Subscriber("/firefly8/state", State, c.stateCB, 7)

    rospy.Subscriber('/firefly1/' + odom_topic_name, Odometry, c.odomCB, 0)
    rospy.Subscriber('/firefly2/' + odom_topic_name, Odometry, c.odomCB, 1)
    rospy.Subscriber('/firefly3/' + odom_topic_name, Odometry, c.odomCB, 2)
    rospy.Subscriber('/firefly4/' + odom_topic_name, Odometry, c.odomCB, 3)
    rospy.Subscriber('/firefly5/' + odom_topic_name, Odometry, c.odomCB, 4)
    rospy.Subscriber('/firefly6/' + odom_topic_name, Odometry, c.odomCB, 5)
    rospy.Subscriber('/firefly7/' + odom_topic_name, Odometry, c.odomCB, 6)
    rospy.Subscriber('/firefly8/' + odom_topic_name, Odometry, c.odomCB, 7)

    rospy.spin()

if __name__ == '__main__':
    rospy.init_node('auto_commands')  
    startNode()
