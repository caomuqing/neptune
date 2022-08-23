#!/usr/bin/env python
import rospy
import numpy as np

from snapstack_msgs.msg import QuadFlightMode
from behavior_selector.srv import MissionModeChange

NOT_FLYING = 0
FLYING = 1

class Behavior_Selector:

    def __init__(self):
        self.status = NOT_FLYING
        self.pubEvent    = rospy.Publisher("globalflightmode", QuadFlightMode, queue_size=1, latch=True)
        self.flightevent = QuadFlightMode()
        
    def sendEvent(self):
        self.flightevent.header.stamp = rospy.get_rostime()
        self.pubEvent.publish(self.flightevent)

    def change_mode(self,req):
        if req.mode == req.KILL:
            self.status = NOT_FLYING
            self.flightevent.mode = QuadFlightMode.KILL
            self.sendEvent()
        
        if req.mode == req.END and self.status == FLYING:
            self.flightevent.mode = QuadFlightMode.LAND
            self.sendEvent()

        if req.mode == req.START:
            self.status = FLYING
            self.flightevent.mode = QuadFlightMode.GO
            self.sendEvent()

    def srvCB(self,req):
        self.change_mode(req)
        return True
                  
def startNode():
    c = Behavior_Selector()
    s = rospy.Service("change_mode", MissionModeChange, c.srvCB)
    rospy.spin()

if __name__ == '__main__':
    rospy.init_node('behavior_selector')
    print "Starting behavior selector"   
    startNode()
