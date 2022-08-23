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

number_of_robots = 5;
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

class auto_commands:

    def __init__(self):
        radius = 10.0;
        one_slice = 2*math.pi/number_of_robots;
        init_pos_as_circle = False;

        self.initialized=False;
        self.start_pub_time = rospy.Time.now();
        self.dist_travelled = [0.0 for i in range(number_of_robots)];
        self.completed_current = [0 for i in range(number_of_robots)];
        self.gotten_odom = [0 for i in range(number_of_robots)];

        self.pos= [[0.0, 0.0, 0.0] for i in range(number_of_robots)];

        self.pos_prev=[[0.0, 0.0, 0.0] for i in range(number_of_robots)];

        self.goals0=  [[-1.5, -2.0, goal_height], #for unity scene config
                    [-1.5, 2.0, goal_height],
                    [-2.0, 0.0, goal_height],
                    [0.0, 10.0, goal_height],
                    [-8.0, -8.0, goal_height],
                    [8.0, -8.0, goal_height],
                    [8.0, 8.0, goal_height],
                    [-8.0, 8.0, goal_height] ]          
        self.goals1=  [[1.5, 2.0, 1.0], #for unity scene config
                    [1.5, -2.0, 1.0],
                    [2.0, 0.0, 1.0],
                    [0.0, -10.0, 1.0],
                    [4.0, 4.0, 1.0],
                    [-4.0, 4.0, 2.0],
                    [-4.0, -4.0, 2.0],
                    [4.0, -4.0, 2.0] ]        
        self.heartbeattime = [rospy.Time.now() for i in range(number_of_robots)];

        if circle_init_bases:
            one_slice = 3.1415927*2/number_of_robots;
            circle_init_dist_to_agent = 2.5;
            circle_init_radius = 10;            
            for i in range(0, number_of_robots):
                theta = one_slice*i;
                base_positions[i][0] = -circle_init_radius*math.cos(theta)- \
                                circle_init_dist_to_agent*math.cos(1.571-theta);
                base_positions[i][1] = -circle_init_radius*math.sin(theta)+ \
                                circle_init_dist_to_agent*math.sin(1.571-theta);   

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

        self.goals = self.goals0
        # self.seed_number = 19; #previous was 9
        # np.random.seed(self.seed_number);

        # random.seed(10);
        self.pose = Pose();
        self.round =1
        self.currentrun = 0
        self.failingcase = 0
        self.pubGoal = [None for i in range(0,number_of_robots)]
        for i in range(0,number_of_robots):
            self.pubGoal[i]=rospy.Publisher('/firefly'+str(i+1)+'/term_goal',PoseStamped,queue_size=1,latch=True)
        self.pubRecord = rospy.Publisher('/firefly/record',PoseStamped,queue_size=3,latch=True)
        self.pubGoal_multi = rospy.Publisher('/start_and_goals',Float32MultiArray,queue_size=1,latch=True)
        self.mtlp_log_received = False;

    #In rospy, the callbacks are all of them in separate threads
    def stateCB(self, data, args):
        idx = args
        # if  self.initialized and (rospy.Time.now()-self.heartbeattime[idx]).to_sec()>15:
        #     print("HEARTBEAT: robot "+str(idx))
        #     self.heartbeattime[idx] = rospy.Time.now();

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
        for i in range(number_of_robots):
            if i!=idx:
                self.goals[i][0] = self.pos[i][0];
                self.goals[i][1] = self.pos[i][1];
                self.goals[i][2] = goal_height;      
            else:
                self.goals[i][0] = data.pose.position.x;
                self.goals[i][1] = data.pose.position.y;
                self.goals[i][2] = goal_height;                               
        self.publishgoals();

    def mtlp_cb(self, data):
        self.mtlp_log_received = True;

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

    def checkAndPublish(self, idx):
        if self.completed_current[idx] != 1:
            self.dist_travelled[idx] += \
                    math.sqrt(sum((self.pos[idx][i]-self.pos_prev[idx][i])**2 for i in range(0,3)));

        if not self.initialized and idx==0 and np.prod(self.gotten_odom) ==1:
            # xx = raw_input("type any to proceed: ");
            self.currentrun = self.currentrun+1
            self.updateGoalsInitial();
            self.publishgoals();
            self.initialized = True;
            return
            # print("INIT: current number of mission completed is "+str(self.currentrun-1))
        
        if not self.completed_current[idx] and \
            LA.norm(np.array(self.pos[idx][0:3])-np.array(self.goals[idx][0:3]))<0.50:
            self.completed_current[idx] = 1;
            # print("Robot "+str(idx)+" has completed current run !")
        elif self.completed_current[idx] and \
            LA.norm(np.array(self.pos[idx][0:3])-np.array(self.goals[idx][0:3]))>0.50:
            self.completed_current[idx] = 0;

        if self.initialized and np.prod(self.completed_current) == 1 and idx==0\
        and self.mtlp_log_received:
            # xx = raw_input("type any to proceed: ");
            # if xx!="s":
            #     return;
            self.updateGoalsRandom();
            self.currentrun = self.currentrun+1            
            self.setupAndPublishRecord(1);
            self.publishgoals();

        if  idx==0 and self.initialized and (rospy.Time.now()-self.start_pub_time).to_sec()>40 \
        and self.currentrun<101 and self.mtlp_log_received:
            # xx = raw_input("type any to proceed: ");
            # if xx!="s":
            #     return;
            # self.seed_number = self.seed_number +1;
            # np.random.seed(self.seed_number);
            self.updateGoalsRandom();
            # self.currentrun = 1           
            # print("CHANGE: resetting seed number to "+str(self.seed_number)+" !!!")
            self.currentrun = self.currentrun+1            
            print("CHANGE: resetting goal points !!!")
            self.failingcase = self.failingcase+1            
            self.setupAndPublishRecord(0);
            self.publishgoals();

    def setupAndPublishRecord(self, success):
        msgRecord=PoseStamped()  
        msgRecord.header.frame_id="world"
        msgRecord.header.stamp = rospy.get_rostime()
        msgRecord.header.seq = self.currentrun-1;
        msgRecord.pose.position.x = (rospy.Time.now()-self.start_pub_time).to_sec();
        msgRecord.pose.position.y = sum(self.dist_travelled)/number_of_robots;
        msgRecord.pose.position.z = success;            
        self.pubRecord.publish(msgRecord);

    def updateGoalsRandom(self):
        for i in range(0, number_of_robots):  
            condition_failed = True;
            # print("updating goals ...")
            while condition_failed:
                condition_failed = False;
                self.goals[i][0] = np.random.uniform(x_min + 4.0, x_max-4.0);
                self.goals[i][1] = np.random.uniform(y_min + 4.0, y_max-4.0);
                if LA.norm(np.array(self.goals[i][0:2])-np.array(base_positions[i][0:2]))>tether_length:
                    condition_failed = True;
                    continue;
                for j in range(0, number_of_robots):
                    if LA.norm(np.array(self.goals[i])-np.array(self.pos[j]))<close_range/4:
                        condition_failed = True;
                for j in range(0, i):
                    if LA.norm(np.array(self.goals[i])-np.array(self.goals[j]))<close_range:
                        condition_failed = True;                
        # print("goals updated!")

    def updateGoalsInitial(self):
        for i in range(0, number_of_robots):  
            self.goals[i][0] = self.pos[i][0];
            self.goals[i][1] = self.pos[i][1];

    def publishgoals(self):
        runno = self.currentrun%2;
        msg_mtlp = Float32MultiArray(); 
        for i in range(0, number_of_robots):
            msg=PoseStamped()  
            msg.pose.orientation = self.pose.orientation
            msg.header.frame_id="world"
            msg.header.stamp = rospy.get_rostime()            
            msg.pose.position.x=self.goals[i][0]
            msg.pose.position.y=self.goals[i][1]
            msg.pose.position.z=self.goals[i][2]   
            # print("goals1 x is %f, y is %f ", msg.pose.position.x, msg.pose.position.y)   

            msg_mtlp.data.append(self.pos[i][0]); #start positions
            msg_mtlp.data.append(self.pos[i][1]);
            msg_mtlp.data.append(goal_height);
            msg_mtlp.data.append(self.goals[i][0]); #start positions
            msg_mtlp.data.append(self.goals[i][1]);
            msg_mtlp.data.append(goal_height);

            if self.currentrun<101:                  
                self.pubGoal[i].publish(msg)
                self.pubGoal_multi.publish(msg_mtlp);
                self.mtlp_log_received = False
            else:
                print("finishing 100 runs: num of failed runs is "+str(self.failingcase))
                quit()
                sys.exit()

            self.dist_travelled[i] = 0.0;  
            self.completed_current[i] = 0;   
        self.start_pub_time = rospy.Time.now();            
        print("current number of mission completed is "+str(self.currentrun-1))



                  
def startNode():
    c = auto_commands()
    #s = rospy.Service("/change_mode",MissionModeChange,c.srvCB)
    for i in range(0,number_of_robots):
        rospy.Subscriber("/firefly"+str(i+1)+"/state", State, c.stateCB, i)
        rospy.Subscriber("/firefly"+str(i+1)+"/goalind", PoseStamped, c.goalOneCB, i)
        rospy.Subscriber('/firefly'+str(i+1)+"/"+odom_topic_name, Odometry, c.odomCB, i)
    rospy.Subscriber("/mtlp_log", PoseStamped, c.mtlp_cb)

    rospy.spin()

if __name__ == '__main__':
    rospy.init_node('auto_commands')  
    startNode()
    print "auto_commands started" 