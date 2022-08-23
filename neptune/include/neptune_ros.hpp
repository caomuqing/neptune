/* ----------------------------------------------------------------------------
 * Nanyang Technological University
 * Authors: Cao Muqing, et al.
 * Acknowledgement: Jesus Tordesillas
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once

#include <Eigen/Dense>

#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/MarkerArray.h>
#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <rviz_visual_tools/rviz_visual_tools.h>
#include <snapstack_msgs/State.h>
#include <snapstack_msgs/Goal.h>
#include <trajectory_msgs/MultiDOFJointTrajectory.h>
#include <nav_msgs/Odometry.h>

#include <mader_msgs/Mode.h>
#include <mader_msgs/DynTraj.h>

#include "utils.hpp"
#include "neptune.hpp"
#include "mader_types.hpp"
#include "kinodynamic_search.hpp"

#include "timer.hpp"
#include "Astar_searcher.h"

namespace rvt = rviz_visual_tools;

class NeptuneRos
{
public:
  NeptuneRos(ros::NodeHandle nh1, ros::NodeHandle nh2, ros::NodeHandle nh3);
  ~NeptuneRos();

private:
  std::unique_ptr<Neptune> neptune_ptr_;
  KinodynamicSearch* kino_ptr_;

  void publishOwnTraj(const mt::PieceWisePol& pwp);
  // class methods
  void pubTraj(const std::vector<mt::state>& data);
  void terminalGoalCB(const geometry_msgs::PoseStamped& msg);
  void pubState(const mt::state& msg, const ros::Publisher pub);
  void stateCB(const snapstack_msgs::State& msg);
  void odomCB(const nav_msgs::Odometry& msg);  
  void modeCB(const mader_msgs::Mode& msg);
  void pubCB(const ros::TimerEvent& e);
  void replanCB(const ros::TimerEvent& e);
  void trajCB(const mader_msgs::DynTraj& msg);
  void autoCMD(const ros::TimerEvent& e);

  void clearMarkerActualTraj();
  void clearMarkerColoredTraj();

  void pubActualTraj();
  void pubSampledAgentsTraj(mt::SampledPointsofCurves input);
  void pubSampledTraj(const std::vector<Eigen::Vector3d>& path, int index=1);
  void pubTetherShapeRough();

  void clearMarkerArray(visualization_msgs::MarkerArray* tmp, ros::Publisher* publisher);

  void visualizeStaticObs(std::vector<mt::Polygon_Std>& convexHullOfStaticObs);

  void setUpCheckingPosAndStaticObs(Eigen::Vector2d pos);

  void updateEntStateStaticObs(Eigen::Vector2d pkplus1);
  void addPointToMarker(visualization_msgs::Marker &m, const Eigen::Vector2d& xy, double z);
  void reactiveController(bool& controller_enabled, Eigen::Vector2d& tarPos, int& itr);

  mt::state state_;

  std::string world_name_ = "world";

  rvt::RvizVisualToolsPtr visual_tools_;

  visualization_msgs::Marker setpoint_;

  ros::NodeHandle nh1_;
  ros::NodeHandle nh2_;
  ros::NodeHandle nh3_;

  ros::Publisher pub_point_G_;
  ros::Publisher pub_point_G_term_;
  ros::Publisher pub_goal_;
  ros::Publisher pub_traj_safe_;
  ros::Publisher pub_setpoint_;
  ros::Publisher pub_actual_traj_;
  ros::Publisher pub_sampled_traj_;
  ros::Publisher pub_sampled_traj2_, pub_sampled_traj3_, pub_sampled_traj4_;
  ros::Publisher pub_sampled_agents_traj_;
  ros::Publisher pub_tether_line_, pub_path_mesh_;
  ros::Publisher pub_log_;

  ros::Publisher pub_traj_safe_colored_;

  ros::Publisher pub_traj_;
  ros::Publisher pub_cmd_;

  ros::Publisher poly_safe_pub_;

  ros::Subscriber sub_goal_;
  ros::Subscriber sub_mode_;
  ros::Subscriber sub_state_;
  ros::Subscriber sub_traj_;
  ros::Subscriber sub_odom_;

  ros::Timer pubCBTimer_;
  ros::Timer replanCBTimer_;
  ros::Time odomCBTime_;
  ros::Timer autoCMDTimer_;

  mt::parameters par_;  // where all the parameters are

  std::string name_drone_;

  std::vector<std::string> traj_;  // trajectory of the dynamic obstacle

  visualization_msgs::MarkerArray traj_safe_colored_;

  int actual_trajID_ = 0;

  int id_;  // id of the drone

  bool published_initial_position_ = false;

  mt::PieceWisePol pwp_last_;

  MADER_timers::Timer timer_stop_;
  MADER_timers::Timer timer_pub_;
  MADER_timers::Timer timer_auto_cmd_;

  double T_span_=0.5;
  mt::state G_;
  mt::state G_prev_;  
  bool gotTerminalGoal_=false;
  int search_process_=0;
  Eigen::VectorXi beta_, beta2_, beta3_;
  std::vector<Eigen::Vector2d> previousCheckingPos_;
  std::vector<Eigen::Vector2d> previousCheckingPosAgent_;
  std::vector<Eigen::Vector2d> latestCheckingPosAgent_;

  int num_of_agents_;
  std::vector<Eigen::Vector2d> pb_;
  std::vector<mt::Polygon_Std> convexHullOfStaticObs_;

  std::vector<Eigen::Matrix<double, 2, 2>> staticObsRep_;

  bool staticObsRepInitialized = false;

  std::vector<Eigen::Vector2i> alphas_;
  std::vector<Eigen::Vector2d> betas_; 
  std::vector<int> active_cases_;

  eu::ent_state entangle_state_;
  eu::ent_state_orig entangle_state_orig_;

  std::vector<std::vector<Eigen::Vector2d>> bendPtsForAgents_, bendPtsForAgents_prev_;  
  std::vector<MADER_timers::Timer> previousCheckingAgentTimer_;

  AstarPathFinder * astar_path_finder_;
  std::vector<mt::Polygon_Std> inflatedStaticObs_;
  std::vector<mt::Polygon_Std> inflatedStaticObsLarger_;
  bool autoCMD_initialized_ = false;
  bool auto_command_ = false;
  unsigned seed_;
  std::vector<Eigen::Vector2d> staticObsLongestDist_;

  double runtime_total_astar_ = 0.0, runtime_length_astar_ = 0.0;
  double runtime_total_kino_ = 0.0, runtime_length_kino_ = 0.0;
  double node_expanded_astar_ = 0.0, node_expanded_kino_ = 0.0;
  int num_successful_runs_ = 0, num_successful_astar_= 0, num_successful_kino_=0;
  int num_targets_sent_ = 0;
  double astar_voxel_scale_ = 2.0;
  std::vector<Eigen::Vector3d> path_taken_;
  double length_of_path_taken_ = 0.0;
  nav_msgs::Odometry log_msg_;
  bool single_benchmark_ = false;
  bool use_vicon_odom_ = false;
  bool just_got_new_goal_ = false;
  bool on_way_back_ = false;
  bool use_moveback_ = false;
  MADER_timers::Timer timer_start_moveback_;

};
