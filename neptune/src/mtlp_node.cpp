/* ----------------------------------------------------------------------------
 * Copyright 2022, Cao Muqing
 * Nanyang Technological University
 * All Rights Reserved
 * Authors: Cao Muqing, et al.
 * -------------------------------------------------------------------------- */

#include <ros/ros.h>
#include <math.h>
#include <std_msgs/Float32MultiArray.h>
#include <nav_msgs/Odometry.h>
#include "solveMultiLinearPath.cpp"
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include "utils.hpp"
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>

ros::Publisher solveDataPub, pathVisualPub_, dataForLogPub_, pub_traj_of_r1_;
std::vector<Eigen::Vector3d> basepoints_;
visualization_msgs::MarkerArray markerarray_just_published_;
int mid_ = 1;
int num_of_robot_ = 3;
mtlp* mtlp_pointer_;
double x_min, x_max, y_min, y_max, z_min, z_max;
ros::Publisher pathIndividualPub[10];
visualization_msgs::MarkerArray markerarray_ind_[10];

void poseCallback(const nav_msgs::Odometry::ConstPtr &msg) {


}

void addPointToMarker(visualization_msgs::Marker &m, const Eigen::Vector3d& xyz)
{
  geometry_msgs::Point p;
  p.x = xyz.x();
  p.y = xyz.y();
  p.z = xyz.z();
  m.points.push_back(p);  
}

void clearMarkerArray(visualization_msgs::MarkerArray* tmp, ros::Publisher* publisher)
{
  if ((*tmp).markers.size() == 0)
  {
    return;
  }

  for (int i = 0; i < (*tmp).markers.size(); i++)
  {
    visualization_msgs::Marker m;
    m.header.frame_id = "world";
    m.header.stamp = ros::Time::now();
    m.type = visualization_msgs::Marker::LINE_LIST;
    m.action = visualization_msgs::Marker::DELETEALL;
    m.id = mid_;
    m.scale.x = 0.15;
    m.scale.y = 0;
    m.scale.z = 0;
    (*tmp).markers[i] = m;
    mid_++;
  }

  (*publisher).publish(*tmp);
  (*tmp).markers.clear();

	for (int i=0; i<num_of_robot_; i++)
	{
	  for (int j = 0; j < markerarray_ind_[i].markers.size(); j++)
	  {
	    visualization_msgs::Marker m;
	    m.header.frame_id = "world";
	    m.header.stamp = ros::Time::now();
	    m.type = visualization_msgs::Marker::LINE_LIST;
	    m.action = visualization_msgs::Marker::DELETEALL;
	    m.id = mid_;
	    m.scale.x = 0.15;
	    m.scale.y = 0;
	    m.scale.z = 0;
	    markerarray_ind_[i].markers[j] = m;
	    mid_++;
	  }

	  pathIndividualPub[i].publish(markerarray_ind_[i]);
	  markerarray_ind_[i].markers.clear();		
	}  
}

void startGoalCallback(const std_msgs::Float32MultiArray::ConstPtr &msg) {
	if (msg->data.size()%6 !=0) ROS_WARN("Size of input message not correct!");

	clearMarkerArray(&markerarray_just_published_, &pathVisualPub_);
	
	int num_of_robot = msg->data.size()/6;
	if (num_of_robot!=num_of_robot_)
	{
		 ROS_WARN("Size of input message not consistent with the number of robots!");
		 return;
	}
	std::vector<Eigen::Matrix<double, 3, 2>> startEndPts;
	std::vector<Eigen::Vector3d> basepoints;
	for (int i=0; i<num_of_robot; i++)
	{
		Eigen::Matrix<double, 3, 2> startEndForRobot;
		for (int k=0; k<2; k++)
		{
			startEndForRobot(0,k) =  msg->data[i*6+k*3+0];
			startEndForRobot(1,k) =  msg->data[i*6+k*3+1];
			startEndForRobot(2,k) =  msg->data[i*6+k*3+2];
		}
		startEndPts.push_back(startEndForRobot);
		basepoints.push_back(basepoints_[i]);
	}
	ROS_INFO("start solving for solveMultiTetherLinearPath...");
	ros::Time time_start = ros::Time::now();
	std::vector<std::vector<Eigen::Vector3d>> outputPath;
	bool result = mtlp_pointer_->solveMultiTetherLinearPath(basepoints, startEndPts, outputPath);
	if (result)
	{
		ROS_INFO("Solving Success!");
		double solve_time = (ros::Time::now()-time_start).toSec();
		ROS_INFO("Solving time is %f ms", solve_time*1000);
		double path_total_length = 0.0;
		for (int i=0; i<num_of_robot; i++)
		{
			for (int j=0; j<outputPath[i].size()-1; j++)
			{
				path_total_length += (outputPath[i][j]-outputPath[i][j+1]).norm();
			}
		}
		double avg_path_length = path_total_length/num_of_robot;
		geometry_msgs::PoseStamped msg;
		msg.header.stamp = ros::Time::now();
		msg.pose.position.x = solve_time*1000;
		msg.pose.position.y = avg_path_length;
		dataForLogPub_.publish(msg);

		// publish path for visualization
	   for (int i = 0; i < outputPath.size(); i++)
	   {
			visualization_msgs::Marker m;
			m.type = visualization_msgs::Marker::LINE_LIST;
			// m.action = visualization_msgs::Marker::DELETEALL;
			m.action = visualization_msgs::Marker::ADD;
			m.id = mid_;
			m.scale.x = 0.10;
			m.scale.y = 0;
			m.scale.z = 0;
			m.ns = "Tether_" + std::to_string(i);
			m.color = mu::getColorJet(i, 0, num_of_robot-1);  // note that par_.v_max is per axis!
    		m.header.stamp = ros::Time::now();
   		m.header.frame_id = "world";

			m.pose.position.x = 0.0;
			m.pose.position.y = 0.0;
			m.pose.position.z = 0.0;
			m.pose.orientation.x = 0.0;
			m.pose.orientation.y = 0.0;
			m.pose.orientation.z = 0.0;
			m.pose.orientation.w = 1.0;	

			visualization_msgs::Marker m_tri = m;
			m_tri.ns = "Tether_tri" + std::to_string(i);
			m_tri.type = visualization_msgs::Marker::TRIANGLE_LIST;
			m_tri.scale.x = 1;
			m_tri.scale.y = 1;
			m_tri.scale.z = 1;			
			std_msgs::ColorRGBA co;
			co.r = 0;
			co.g = 102/255;
			co.b = 1.0;
			co.a = 1;			
			m_tri.color = co;
			// addPointToMarker(m, basepoints_[i]);
			// addPointToMarker(m, outputPath[i][0]);
			addPointToMarker(m, outputPath[i][0]);
			for (int j=1; j<outputPath[i].size()-1; j++)
			{
				addPointToMarker(m, outputPath[i][j]);
				addPointToMarker(m, outputPath[i][j]);
				// addPointToMarker(m, basepoints_[i]);
				// addPointToMarker(m, outputPath[i][j]);
			}			
			addPointToMarker(m, outputPath[i].back());
			// addPointToMarker(m, outputPath[i].back());
			// addPointToMarker(m, basepoints_[i]);

			markerarray_just_published_.markers.push_back(m);
			markerarray_ind_[i].markers.push_back(m);

			mid_++;

			visualization_msgs::Marker m2;
			m2.type = visualization_msgs::Marker::LINE_LIST;
			m2.action = visualization_msgs::Marker::ADD;
			m2.id = mid_;
			m2.scale.x = 0.05;
			m2.scale.y = 0;
			m2.scale.z = 0;
			m2.ns = "Tether_" + std::to_string(i);
			std_msgs::ColorRGBA c;
			c.r = 128;
			c.g = 128;
			c.b = 128;
			c.a = 0.8;			
			m2.color = c;
    		m2.header.stamp = ros::Time::now();
   		m2.header.frame_id = "world";

			m2.pose.orientation.x = 0.0;
			m2.pose.orientation.y = 0.0;
			m2.pose.orientation.z = 0.0;
			m2.pose.orientation.w = 1.0;	

			// addPointToMarker(m2, basepoints_[i]);
			// addPointToMarker(m2, outputPath[i][0]);
			// for (int j=1; j<outputPath[i].size()-1; j++)
			// {
			// 	addPointToMarker(m2, outputPath[i][j]);
			// 	addPointToMarker(m2, basepoints_[i]);
			// }
			addPointToMarker(m2, outputPath[i].back());
			addPointToMarker(m2, basepoints_[i]);

			for (int j=0; j<outputPath[i].size()-1; j++)
			{
				addPointToMarker(m_tri, outputPath[i][j]);
				addPointToMarker(m_tri, outputPath[i][j+1]);
				addPointToMarker(m_tri, basepoints_[i]);

			}
			// markerarray_just_published_.markers.push_back(m2); //uncomment to show cable
			// markerarray_ind_[i].markers.push_back(m2);		//uncomment to show cable
			// if (!m_tri.points.empty()) markerarray_ind_[i].markers.push_back(m_tri);

			pathIndividualPub[i].publish(markerarray_ind_[i]);
			mid_++;			
	   }	
	   pathVisualPub_.publish(markerarray_just_published_);	

		nav_msgs::Path traj1;
		traj1.header.stamp = ros::Time::now();
		traj1.header.frame_id = "world";
		geometry_msgs::PoseStamped temp_path;
		for (int i=0; i<outputPath[0].size(); i++)
		{
			temp_path.pose.position.x = outputPath[0][i](0);
			temp_path.pose.position.y = outputPath[0][i](1);
			temp_path.pose.position.z = outputPath[0][i](2);
			temp_path.pose.orientation.w = 1;
			temp_path.pose.orientation.x = 0;
			temp_path.pose.orientation.y = 0;
			temp_path.pose.orientation.z = 0;
			traj1.poses.push_back(temp_path);
		}
		pub_traj_of_r1_.publish(traj1);
	}
	else
	{
		ROS_WARN("Solving Failed!");
	}
}

int main(int argc, char *argv[]) {
	// initialize ROS
	ros::init(argc, argv, "solve_multi_linear");
	ros::NodeHandle nh("~");


	 std::vector<double> base_pos;
	 if (!nh.getParam("base_positions", base_pos))
	 {
	 	ROS_WARN("not getting param! exitting");
	 	exit(-1);
	 }

	bool circle_init_bases;
  	mu::safeGetParam(nh, "circle_init_bases", circle_init_bases);
  	mu::safeGetParam(nh, "num_of_agents", num_of_robot_);
	mu::safeGetParam(nh, "x_min", x_min);
	mu::safeGetParam(nh, "x_max", x_max);

	mu::safeGetParam(nh, "y_min", y_min);
	mu::safeGetParam(nh, "y_max", y_max);

	mu::safeGetParam(nh, "z_min", z_min);
	mu::safeGetParam(nh, "z_max", z_max);
	double radius = 0.0;
	mu::safeGetParam(nh, "drone_radius", radius);
	double voxel_size;
  	mu::safeGetParam(nh, "a_star_fraction_voxel_size", voxel_size);


	double circle_init_dist_to_agent = 2.5;
  	double circle_init_radius = 10;

	basepoints_.clear();

	if (circle_init_bases)
	{
  		double one_slice = 3.1415927*2/num_of_robot_;
		for (int i=0; i<num_of_robot_; i++)
		{
			double theta = one_slice*i;
			Eigen::Vector3d tmpv(0, 0, 0);
			tmpv(0) = -circle_init_radius*cosf(theta)- 
			                circle_init_dist_to_agent*cosf(1.571-theta);
			tmpv(1) = -circle_init_radius*sinf(theta)+ 
			                circle_init_dist_to_agent*sinf(1.571-theta);   
			basepoints_.push_back(tmpv);                   
		}
	}
	else
	{
		for (int i=0; i<base_pos.size()/2; i++)
		{
		    basepoints_.push_back(Eigen::Vector3d(base_pos[i*2], base_pos[i*2+1], 0));
		}		
	}	
	// std::cout<<"base positions are "<<std::endl;
	// std::cout<<basepoints_<<std::endl;
	mtlp_pointer_ = new mtlp(voxel_size, x_min, x_max, y_min, y_max, z_min, z_max, radius);
	ros::Subscriber sub_state_est = nh.subscribe("flight_pilot/state_estimate", 1, &poseCallback);
	ros::Subscriber startgoals_sub = nh.subscribe("/start_and_goals", 1, &startGoalCallback);
  	pathVisualPub_ = nh.advertise<visualization_msgs::MarkerArray>("multi_linear_path", 1);
  	dataForLogPub_ = nh.advertise<geometry_msgs::PoseStamped>("/mtlp_log", 1);
  	pub_traj_of_r1_ = nh.advertise<nav_msgs::Path>("/path_mtlp_r1", 1);
  	// ros::Timer timer_main_loop_ = pnh.createTimer(ros::Duration(0.10), mainLoopCallback);
	for (int i=0; i<num_of_robot_;i++){
		pathIndividualPub[i]=nh.advertise<visualization_msgs::MarkerArray>("/mtlp/firefly"+ std::to_string(i+1), 50);
	}
  	ros::spin();
	ros::waitForShutdown();
 
		// ros::spinOnce();
	// }
	return 0;
}