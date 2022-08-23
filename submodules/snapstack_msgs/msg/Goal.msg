# Use this message to command the outer loop to track
# a trajectory generated from a high-level trajectory planner.

Header header

# Current time-slice of desired trajectory
geometry_msgs/Vector3 p # position
geometry_msgs/Vector3 v # velocity
geometry_msgs/Vector3 a # acceleration
geometry_msgs/Vector3 j # jerk

float64 yaw # heading / yaw angle
float64 dyaw # d/dt{unrolled, unpitched body heading w.r.t world}
# n.b., recall that dyaw = d/dt{psi} != r. Angular heading rate r is defined in
# the body frame, but yaw is the heading of the local level frame w.r.t world.
# For slow, nearly-level flight, dyaw ~= r. For more agile flight, it will
# be useful to make sure you are commanding the correct quantity.
# See, e.g., eq (7) and (8) in 
# https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?article=2324&context=facpub

bool power # true if motors should be able to spin

# Trajectory tracking mode constants
uint8 MODE_POSITION_CONTROL     = 0
uint8 MODE_VELOCITY_CONTROL     = 1
uint8 MODE_ACCELERATION_CONTROL = 2

# Trajectory tracking mode for x/y and z components.
# The default is POSITION control, which uses position and velocity error
# to calculate the control effort. VELOCITY control only uses vel error.
# ACCELERATION mode does not use tracking error and could be used to provide
# a control signal computed from something other than the default PID cntrl.
uint8 mode_xy
uint8 mode_z
