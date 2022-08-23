/* ----------------------------------------------------------------------------
 * Nanyang Technological University
 * Authors: Cao Muqing, et al.
 * Acknowledgement: Jesus Tordesillas
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "neptune_ros.hpp"

#include <sensor_msgs/point_cloud_conversion.h>
#include <sensor_msgs/point_cloud_conversion.h>
#include <nav_msgs/Path.h>

#include <decomp_ros_msgs/PolyhedronArray.h>
#include <decomp_ros_utils/data_ros_utils.h>  //For DecompROS::polyhedron_array_to_ros
#include <decomp_geometry/polyhedron.h>       //For hyperplane
#include <Eigen/Geometry>

#include <stdio.h>
#include <random>
#include <shape_msgs/Mesh.h>

typedef MADER_timers::Timer MyTimer;

// this object is created in the mader_ros_node
NeptuneRos::NeptuneRos(ros::NodeHandle nh1, ros::NodeHandle nh2, ros::NodeHandle nh3) : nh1_(nh1), nh2_(nh2), nh3_(nh3)
{
  mu::safeGetParam(nh1_, "visual", par_.visual);
  mu::safeGetParam(nh1_, "color_type", par_.color_type);

  mu::safeGetParam(nh1_, "dc", par_.dc);
  mu::safeGetParam(nh1_, "goal_radius", par_.goal_radius);
  mu::safeGetParam(nh1_, "drone_radius", par_.drone_radius);
  mu::safeGetParam(nh1_, "rad_term_pos", par_.rad_term_pos);

  mu::safeGetParam(nh1_, "w_max", par_.w_max);
  mu::safeGetParam(nh1_, "alpha_filter_dyaw", par_.alpha_filter_dyaw);

  mu::safeGetParam(nh1_, "x_min", par_.x_min);
  mu::safeGetParam(nh1_, "x_max", par_.x_max);

  mu::safeGetParam(nh1_, "y_min", par_.y_min);
  mu::safeGetParam(nh1_, "y_max", par_.y_max);

  mu::safeGetParam(nh1_, "z_min", par_.z_min);
  mu::safeGetParam(nh1_, "z_max", par_.z_max);

  std::vector<double> v_max_tmp;
  std::vector<double> a_max_tmp;

  mu::safeGetParam(nh1_, "v_max", v_max_tmp);
  mu::safeGetParam(nh1_, "a_max", a_max_tmp);
  mu::safeGetParam(nh1_, "j_max", par_.j_max);

  par_.v_max << v_max_tmp[0], v_max_tmp[1], v_max_tmp[2];
  par_.a_max << a_max_tmp[0], a_max_tmp[1], a_max_tmp[2];

  mu::safeGetParam(nh1_, "factor_alpha", par_.factor_alpha);

  mu::safeGetParam(nh1_, "num_pol", par_.num_pol);
  mu::safeGetParam(nh1_, "deg_pol", par_.deg_pol);
  mu::safeGetParam(nh1_, "weight", par_.weight);

  mu::safeGetParam(nh1_, "upper_bound_runtime_snlopt", par_.upper_bound_runtime_snlopt);
  mu::safeGetParam(nh1_, "lower_bound_runtime_snlopt", par_.lower_bound_runtime_snlopt);
  mu::safeGetParam(nh1_, "runtime_opt", par_.runtime_opt);

  mu::safeGetParam(nh1_, "a_star_samp_x", par_.a_star_samp_x);

  mu::safeGetParam(nh1_, "a_star_fraction_voxel_size", par_.a_star_fraction_voxel_size);

  mu::safeGetParam(nh1_, "a_star_bias", par_.a_star_bias);

  mu::safeGetParam(nh1_, "basis", par_.basis);

  mu::safeGetParam(nh1_, "res_plot_traj", par_.res_plot_traj);

  mu::safeGetParam(nh1_, "num_of_agents", par_.num_of_agents);
  mu::safeGetParam(nh1_, "tetherLength", par_.tetherLength);
  mu::safeGetParam(nh1_, "num_sample_per_interval", par_.num_sample_per_interval);
  mu::safeGetParam(nh1_, "T_span", par_.T_span);
  mu::safeGetParam(nh1_, "safe_factor", par_.safe_factor);
  mu::safeGetParam(nh1_, "use_linear_constraints", par_.use_linear_constraints);

  std::vector<double> base_pos;
  mu::safeGetParam(nh1_, "base_positions", base_pos);
  mu::safeGetParam(nh1_, "num_of_static_obst", par_.num_of_static_obst);
  std::vector<int> numOfVertEachObst;
  mu::safeGetParam(nh1_, "num_of_vertices_each_obst", numOfVertEachObst);
  std::vector<double> vertice_obst_x, vertice_obst_y;
  mu::safeGetParam(nh1_, "vertice_obst_x", vertice_obst_x);
  mu::safeGetParam(nh1_, "vertice_obst_y", vertice_obst_y);
  bool random_obst;
  mu::safeGetParam(nh1_, "random_obst", random_obst);
  mu::safeGetParam(nh1_, "auto_cmd", auto_command_);
  mu::safeGetParam(nh1_, "use_not_reaching_soln", par_.use_not_reaching_soln);
  bool circle_init_bases;
  mu::safeGetParam(nh1_, "circle_init_bases", circle_init_bases);
  mu::safeGetParam(nh1_, "enable_entangle_check", par_.enable_entangle_check);
  mu::safeGetParam(nh1_, "enable_single_benchmark", single_benchmark_);
  nh1_.param("use_vicon_odom", use_vicon_odom_, false);
  nh1_.param("use_DIRECT_backend", par_.use_direct, false);
  nh1_.param("use_moveback_strategy", use_moveback_, false);
  nh1_.param("move_when_goal_occupied", par_.move_when_goal_occupied, false);  

  double replan_dt;
  nh1_.param("replan_dt", replan_dt, par_.dc);

  std::cout << "Parameters obtained" << std::endl;

  // CHECK parameters
  std::cout << bold << "Parameters obtained, checking them..." << reset << std::endl;
  std::cout << bold << "Parameters checked" << reset << std::endl;
  std::cout << bold << "Running Neptune Ros!!!" << reset << std::endl;

  /////////////////////
  name_drone_ = ros::this_node::getNamespace();  // Return also the slashes (2 in Kinetic, 1 in Melodic)
  name_drone_.erase(std::remove(name_drone_.begin(), name_drone_.end(), '/'), name_drone_.end());  // Remove the slashes

  std::string id = name_drone_;
  id.erase(0, 7);  // Erase firefly
  id_ = std::stoi(id);
  // id_ = (int)name_drone_.back() - '0';
  std::cout<<bold<<green<<"id is "<<id_<<std::endl;
  if (id_ > par_.num_of_agents || id_<=0)
  {
    std::cout << bold <<red<< "Drone index exceeding the total number of robots, "<<
                              "exitting!!!" << reset << std::endl;
    exit(-1);
  }

  par_.id = id_;
  num_of_agents_ = par_.num_of_agents;
  for (int i=0; i<num_of_agents_; i++)
  {
    par_.pb.push_back(Eigen::Vector2d(base_pos[i*2], base_pos[i*2+1]));
  }
  double one_slice = 3.1415927*2/num_of_agents_;
  double circle_init_radius = 10;
  double circle_init_dist_to_agent = 2.5;
  if (circle_init_bases)
  {
    for (int i=0; i<num_of_agents_; i++)
    {
      double theta = one_slice*i;
      par_.pb[i](0) = -circle_init_radius*cosf(theta)- 
                      circle_init_dist_to_agent*cosf(1.571-theta);
      par_.pb[i](1) = -circle_init_radius*sinf(theta)+ 
                      circle_init_dist_to_agent*sinf(1.571-theta);                      
    }
  }
  neptune_ptr_ = std::unique_ptr<Neptune>(new Neptune(par_));

  // Publishers
  pub_goal_ = nh1_.advertise<snapstack_msgs::Goal>("goal", 1);
  pub_cmd_ = nh1_.advertise<trajectory_msgs::MultiDOFJointTrajectory>("command/trajectory", 1);

  pub_setpoint_ = nh1_.advertise<visualization_msgs::Marker>("setpoint", 1);
  pub_point_G_ = nh1_.advertise<geometry_msgs::PointStamped>("point_G", 1);
  pub_point_G_term_ = nh1_.advertise<geometry_msgs::PointStamped>("point_G_term", 1);
  pub_actual_traj_ = nh1_.advertise<visualization_msgs::Marker>("actual_traj", 1);
  pub_sampled_agents_traj_ = nh1_.advertise<visualization_msgs::Marker>("sampled_agents_traj", 1);
  pub_tether_line_ = nh1_.advertise<visualization_msgs::Marker>("tether_rough", 1);
  pub_path_mesh_ = nh1_.advertise<visualization_msgs::Marker>("path_mesh", 1);
  pub_sampled_traj_ = nh1_.advertise<nav_msgs::Path>("sampled_traj", 1);
  pub_sampled_traj2_ = nh1_.advertise<nav_msgs::Path>("sampled_traj2", 1);
  pub_sampled_traj3_ = nh1_.advertise<nav_msgs::Path>("sampled_traj3", 1);
  pub_sampled_traj4_ = nh1_.advertise<nav_msgs::Path>("sampled_traj4", 1);

  poly_safe_pub_ = nh1_.advertise<decomp_ros_msgs::PolyhedronArray>("poly_safe", 1, true);
  pub_traj_safe_colored_ = nh1_.advertise<visualization_msgs::MarkerArray>("traj_safe_colored", 1);
  pub_traj_ = nh1_.advertise<mader_msgs::DynTraj>("/trajs", 1, true);  // The last boolean is latched or not
  pub_log_ = nh1_.advertise<nav_msgs::Odometry>("/firefly/log_for_plot", 1);

  // Subscribers
  sub_goal_ = nh1_.subscribe("term_goal", 1, &NeptuneRos::terminalGoalCB, this);
  sub_mode_ = nh1_.subscribe("mode", 1, &NeptuneRos::modeCB, this);
  sub_state_ = nh1_.subscribe("state", 1, &NeptuneRos::stateCB, this);
  sub_traj_ = nh1_.subscribe("/trajs", 20, &NeptuneRos::trajCB, this);  // The number is the queue size
  sub_odom_ = nh1_.subscribe("ground_truth/odometry", 1, &NeptuneRos::odomCB, this);

  // Timers
  pubCBTimer_ = nh2_.createTimer(ros::Duration(par_.dc), &NeptuneRos::pubCB, this);
  autoCMDTimer_ = nh2_.createTimer(ros::Duration(0.5), &NeptuneRos::autoCMD, this);  
  replanCBTimer_ = nh3_.createTimer(ros::Duration(replan_dt), &NeptuneRos::replanCB, this);

  // For now stop all these subscribers/timers until we receive GO
  // sub_state_.shutdown();
  pubCBTimer_.stop();
  replanCBTimer_.stop();
  autoCMDTimer_.stop();

  // Rviz_Visual_Tools
  visual_tools_.reset(new rvt::RvizVisualTools("world", "/rviz_visual_tools"));
  visual_tools_->loadMarkerPub();  // create publisher before waitin
  ros::Duration(0.5).sleep();
  visual_tools_->deleteAllMarkers();
  visual_tools_->enableBatchPublishing();

  // Markers
  setpoint_ = mu::getMarkerSphere(0.35, mu::orange_trans);

  timer_stop_.Reset();
  timer_pub_.Reset();
  timer_auto_cmd_.Reset();

  clearMarkerActualTraj();

  convexHullOfStaticObs_.clear();
  seed_ = std::chrono::system_clock::now().time_since_epoch().count();

  if (random_obst)
  {
    std::vector<Eigen::Vector2d> addedpoints;
    std::uniform_real_distribution<double> unif_x(par_.x_min, par_.x_max);
    std::uniform_real_distribution<double> unif_y(par_.y_min, par_.y_max);

    std::default_random_engine re = std::default_random_engine(seed_);

    for (int i=0; i<par_.num_of_static_obst; i++)
    {
      bool valid_pt = false;
      Eigen::Matrix<double, 2, 4> vertices;      

      while (!valid_pt)
      {
        double ran_x = unif_x(re);    
        double ran_y = unif_y(re);    
        for (int j=0; j<addedpoints.size(); j++)
        {
          double dist_to_an_obst = (Eigen::Vector2d(ran_x, ran_y) - addedpoints[j]).norm();
          if (dist_to_an_obst < 3 * par_.a_star_fraction_voxel_size ||
            dist_to_an_obst > 2 * par_.drone_radius &&
            dist_to_an_obst<par_.a_star_fraction_voxel_size*2.83 + 8 * par_.drone_radius)
          {
            goto trynext;
          }
        }
        if ((Eigen::Vector2d(ran_x, ran_y) - par_.pb[id_-1]).norm()<6.0)
          goto trynext;
        valid_pt = true;
        vertices << ran_x + 0.25, ran_x + 0.25, ran_x - 0.25, ran_x - 0.25,
                    ran_y + 0.25, ran_y - 0.25, ran_y - 0.25, ran_y + 0.25;   
        convexHullOfStaticObs_.push_back(vertices);                  
        addedpoints.push_back(Eigen::Vector2d(ran_x, ran_y));
        trynext:;                 
      }

    }    
  }
  else
  {
    int _idx = 0;
    for (int i=0; i<par_.num_of_static_obst; i++)
    {
      Eigen::MatrixXd vertices(2, numOfVertEachObst[i]);
      for (int j=0; j<numOfVertEachObst[i]; j++)
      {
        vertices.col(j) = Eigen::Vector2d(vertice_obst_x[_idx], vertice_obst_y[_idx]);
        _idx++;
      }
      convexHullOfStaticObs_.push_back(vertices);
    }    
  }
  neptune_ptr_->setStaticObst(convexHullOfStaticObs_);
  visualizeStaticObs(convexHullOfStaticObs_);

  /******* kinodynamic searcher setip ************/
  T_span_ = par_.T_span;
  pb_ = par_.pb;

  beta_ = Eigen::MatrixXi::Zero(num_of_agents_ + par_.num_of_static_obst, 1);
  beta2_ = Eigen::MatrixXi::Zero(num_of_agents_ + par_.num_of_static_obst, 1);
  beta3_ = Eigen::MatrixXi::Zero(num_of_agents_ + par_.num_of_static_obst, 1);

  Eigen::Vector2d _pos(-1000, -1000);
  MADER_timers::Timer _timer(true);

  for (int i=0; i<num_of_agents_; i++)
  {
    previousCheckingPosAgent_.push_back(_pos);
    previousCheckingAgentTimer_.push_back(_timer);
    latestCheckingPosAgent_.push_back(_pos);
  }

  for (int i=0; i<num_of_agents_ + par_.num_of_static_obst; i++)
  {
    entangle_state_.active_cases.push_back(0);
  }  

  for (int i=0; i<num_of_agents_; i++)
  {
    std::vector<Eigen::Vector2d> bendPts;
    bendPtsForAgents_.push_back(bendPts);
    bendPtsForAgents_prev_.push_back(bendPts);
  }  

  /******* end of kinodynamic searcher setip ******/  
  
  /**************initialize kino_ptr for testing************/
  kino_ptr_ = new KinodynamicSearch(par_.num_pol, par_.deg_pol, id_, par_.safe_factor, T_span_, 
                                    par_.num_sample_per_interval, pb_, par_.use_not_reaching_soln, 
                                    par_.enable_entangle_check);
  kino_ptr_->setTetherLength(par_.tetherLength);
  kino_ptr_->setMaxValuesAndSamples(v_max_tmp[0], a_max_tmp[0], par_.j_max, 
                                    par_.a_star_samp_x);
  kino_ptr_->setXYZMinMaxAndRa(par_.x_min, par_.x_max, par_.y_min, par_.y_max, par_.z_min, 
                               par_.z_max, 5.0, par_.a_star_fraction_voxel_size);
  kino_ptr_->setRunTime(15.0 * 1000.0);
  kino_ptr_->setBias(1.1);
  kino_ptr_->setGoalSize(par_.goal_radius);

  inflatedStaticObs_.clear();
  inflatedStaticObsLarger_.clear();
  double safe_dist = 2 * par_.drone_radius;

  for (int i=0; i<convexHullOfStaticObs_.size(); i++)
  {
    std::vector<Point_2> points_cgal;
    std::vector<Point_2> points_cgal_larger;

    for (int j=0; j<convexHullOfStaticObs_[i].cols(); j++)
    {
      points_cgal.push_back(Point_2(convexHullOfStaticObs_[i](0,j) + safe_dist, 
                                    convexHullOfStaticObs_[i](1,j) + safe_dist));
      points_cgal.push_back(Point_2(convexHullOfStaticObs_[i](0,j) + safe_dist, 
                                    convexHullOfStaticObs_[i](1,j) - safe_dist));
      points_cgal.push_back(Point_2(convexHullOfStaticObs_[i](0,j) - safe_dist, 
                                    convexHullOfStaticObs_[i](1,j) + safe_dist));
      points_cgal.push_back(Point_2(convexHullOfStaticObs_[i](0,j) - safe_dist,
                                    convexHullOfStaticObs_[i](1,j) - safe_dist));
      points_cgal_larger.push_back(Point_2(convexHullOfStaticObs_[i](0,j) + 2*safe_dist, 
                                           convexHullOfStaticObs_[i](1,j) + 2*safe_dist));
      points_cgal_larger.push_back(Point_2(convexHullOfStaticObs_[i](0,j) + 2*safe_dist, 
                                           convexHullOfStaticObs_[i](1,j) - 2*safe_dist));
      points_cgal_larger.push_back(Point_2(convexHullOfStaticObs_[i](0,j) - 2*safe_dist, 
                                           convexHullOfStaticObs_[i](1,j) + 2*safe_dist));
      points_cgal_larger.push_back(Point_2(convexHullOfStaticObs_[i](0,j) - 2*safe_dist,
                                           convexHullOfStaticObs_[i](1,j) - 2*safe_dist));
    }
    inflatedStaticObs_.push_back(cu::convexHullOfPoints2d(points_cgal));
    inflatedStaticObsLarger_.push_back(cu::convexHullOfPoints2d(points_cgal_larger));
  }
  
  kino_ptr_->setStaticObstVert(inflatedStaticObs_);  

  astar_path_finder_  = new AstarPathFinder();
  astar_path_finder_  -> init(par_.a_star_fraction_voxel_size*astar_voxel_scale_, 
                              par_.x_min, par_.x_max, 
                              par_.y_min, par_.y_max, par_.z_min, par_.z_max, 
                              par_.tetherLength);
  astar_path_finder_ -> setStaticObstVert(inflatedStaticObs_, convexHullOfStaticObs_);
  entangle_state_orig_.contPts.push_back(pb_[id_-1]);

  /**************end for initialize kino_ptr for testing************/


  ////// to avoid having to click on the GUI (TODO)
  mader_msgs::Mode tmp;
  tmp.mode = 1;
  modeCB(tmp);
  // ros::Duration(1.0).sleep();  // TODO
  // bool success_service_call = system("rosservice call /change_mode 'mode: 1'");
  ////
  odomCBTime_ = ros::Time::now();
  ROS_INFO("Planner initialized");
}

NeptuneRos::~NeptuneRos()
{
  sub_state_.shutdown();
  pubCBTimer_.stop();
  replanCBTimer_.stop();
  autoCMDTimer_.stop();

}


void NeptuneRos::trajCB(const mader_msgs::DynTraj& msg)
{
  if (msg.id == id_ || msg.id > num_of_agents_)
  {  // This is my own trajectory
    return;
  }
  // std::cout<<bold<<blue<<"receiving traj fro agent id "<<msg.id<<std::endl;

  Eigen::Vector3d W_pos(msg.pos.x, msg.pos.y, msg.pos.z);  // position in world frame
  double dist = (state_.pos - W_pos).norm();

  mt::dynTraj tmp;
  tmp.function.push_back(msg.function[0]);
  tmp.function.push_back(msg.function[1]);
  tmp.function.push_back(msg.function[2]);

  tmp.bbox << msg.bbox[0], msg.bbox[1], msg.bbox[2];

  tmp.id = msg.id;

  tmp.is_agent = msg.is_agent;

  if (msg.is_agent)
  {
    tmp.pwp = mu::pwpMsg2Pwp(msg.pwp);
  }

  tmp.time_received = ros::Time::now().toSec();

  neptune_ptr_->updateTrajObstacles(tmp);

  // std::cout<<bold<<blue<<"receiving traj is"<<msg.id<<std::endl;
  /**********checking for entangling beta ***************/

  latestCheckingPosAgent_[msg.id-1] = Eigen::Vector2d(msg.pos.x, msg.pos.y);

  std::vector<Eigen::Vector2d> bendptVec;

  for (int i=0; i<msg.bendpt.size(); i++)
  {
    Eigen::Vector2d bentpt(msg.bendpt[i].x, msg.bendpt[i].y);
    bendptVec.push_back(bentpt);
  }
  bendPtsForAgents_[msg.id-1] = bendptVec; // save in a vector
  if (bendPtsForAgents_prev_[msg.id-1].empty()) bendPtsForAgents_prev_[msg.id-1] = bendptVec; 

  if (previousCheckingPosAgent_[msg.id-1](0)<-900) //first time setup
  {
    previousCheckingPosAgent_[msg.id-1] = Eigen::Vector2d(msg.pos.x, msg.pos.y);
    neptune_ptr_->setBendPtsForAgent(msg.id, bendptVec); 
  }  
}

// This trajectory contains all the future trajectory (current_pos --> A --> final_point_of_traj), because it's the
// composition of pwp
void NeptuneRos::publishOwnTraj(const mt::PieceWisePol& pwp)
{
  std::vector<std::string> s;  // mu::pieceWisePol2String(pwp); The rest of the agents will use the pwp field, not the
                               // string
  s.push_back("");
  s.push_back("");
  s.push_back("");

  mader_msgs::DynTraj msg;
  msg.function = s;
  msg.bbox.push_back(2 * par_.drone_radius);
  msg.bbox.push_back(2 * par_.drone_radius);
  msg.bbox.push_back(2 * par_.drone_radius);
  msg.pos.x = state_.pos.x();
  msg.pos.y = state_.pos.y();
  msg.pos.z = state_.pos.z();
  msg.id = id_;

  msg.is_agent = true;

  msg.pwp = mu::pwp2PwpMsg(pwp);

  geometry_msgs::Vector3 v3;
  v3.x = pb_[id_-1].x();
  v3.y = pb_[id_-1].y();
  v3.z = 0.0;
  msg.bendpt.push_back(v3);

  for (int i=0; i<entangle_state_.bendPointsIdx.size(); i++)
  {
    Eigen::Vector2i _id = entangle_state_.alphas[entangle_state_.bendPointsIdx[i]];
    if (_id(0) < num_of_agents_) // if agent, only base is possible
    {
      v3.x = pb_[_id(0)-1].x();
      v3.y = pb_[_id(0)-1].y();
    }
    else //if static obst, depends on the case
    {
      v3.x = staticObsRep_[_id(0)-num_of_agents_-1].col(_id(1)).x();
      v3.y = staticObsRep_[_id(0)-num_of_agents_-1].col(_id(1)).y();      
    }
    msg.bendpt.push_back(v3);
  }
  // std::cout<<"msg.pwp.times[0]= "<<msg.pwp.times[0];

  pub_traj_.publish(msg);
}

void NeptuneRos::replanCB(const ros::TimerEvent& e)
{
  if (ros::ok() && published_initial_position_ == true)
  {
    // std::cout<<bold<<white<<"in replan cb!"<<std::endl;
    bool just_got_new_goal = just_got_new_goal_;
    if (just_got_new_goal_) just_got_new_goal_ = false;

    /*-------------------------single agent benchmarking---------------------------------*/
    if (single_benchmark_)
    {
      if (gotTerminalGoal_)//||getchar()=='s')
      {
        int search_result;
        std::vector<Eigen::Vector3d> currentSamplePath;
        astar_path_finder_ -> setHsig(entangle_state_orig_);
        astar_path_finder_-> AstarGraphSearch(G_prev_.pos.head<2>(), G_.pos.head<2>(), search_result);

        int search_result_kino;
        std::vector<Eigen::Vector3d> currentSamplePath_kino;
        kino_ptr_->run(currentSamplePath_kino, search_result_kino);

        if (search_result==1) 
        {
          // std::vector<Eigen::Vector2d> v = astar_path_finder_-> getVisitedNodes();
          // visual_tools_->deleteAllMarkers();

          // for (int i=0; i<v.size(); i++)
          // {
          //   auto color = visual_tools_->getRandColor();
          //   Eigen::Isometry3d pose;
          //   pose.translation() = Eigen::Vector3d(v[i](0), v[i](1), 0.0);

          //   // Calculate the rotation matrix from the original normal z_0 = (0,0,1) to new normal n = (A,B,C)
          //   Eigen::Vector3d z_0 = Eigen::Vector3d::UnitZ();
          //   Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(z_0, z_0);
          //   pose.linear() = q.toRotationMatrix();

          //   visual_tools_->publishCuboid(pose, par_.a_star_fraction_voxel_size*astar_voxel_scale_, 
                                            // par_.a_star_fraction_voxel_size*astar_voxel_scale_, 
          //                                par_.a_star_fraction_voxel_size, color);

          //   visual_tools_->trigger();          
          //   ros::Duration(0.005).sleep();

          // }        

          if (search_result_kino!=1)
          {
            int current_process=2;
            kino_ptr_->clearProcess();

            while (current_process==2)
            {
              std::vector<Eigen::Vector3d> cSamplePath;
              kino_ptr_->runonce(cSamplePath, current_process);
              pubSampledTraj(cSamplePath, 4);
              ros::Duration(0.005).sleep();
            }
            // getchar();
          }         
        }
        if (search_result_kino==1)
        {
          // the case when original astar is not feasible at start point
          // if (eu::tetherLength_orig(entangle_state_orig_.contPts) > par_.tetherLength)
          //   search_process_=1;

        } 

        if (search_result==1 && search_result_kino==1) 
        {
          num_successful_astar_++;
          double runtime_total, runtime_length_check;
          int node_used = 0;

          astar_path_finder_->getRuntime(runtime_total, runtime_length_check, node_used);
          // std::cout<<bold<<green<<"[homotopic Astar] runtime is "<<runtime_total
          // <<", time spent for length check is "<<runtime_length_check<<std::endl;
          runtime_total_astar_ = runtime_total_astar_*(num_successful_astar_-1)/num_successful_astar_ +
                                 runtime_total/num_successful_astar_;
          runtime_length_astar_ = runtime_length_astar_*(num_successful_astar_-1)/num_successful_astar_ +
                                 runtime_length_check/num_successful_astar_;                               
          node_expanded_astar_ = node_expanded_astar_*(num_successful_astar_-1)/num_successful_astar_ +
                                 (double)node_used/num_successful_astar_;   
          // astar_path_finder_->getHsigTerminal(entangle_state_orig_);
          // G_prev_ = G_;
          log_msg_.pose.pose.position.x = runtime_total;
          log_msg_.pose.pose.position.y = runtime_length_check;
          log_msg_.pose.pose.position.z = node_used;


          num_successful_kino_++;    
          kino_ptr_->getRuntime(runtime_total, runtime_length_check, node_used);        
          // std::cout<<bold<<green<<"[Kinodynamic Astar] runtime is "<<runtime_total
          // <<", time spent for length check is "<<runtime_length_check<<std::endl;

          runtime_total_kino_ = runtime_total_kino_*(num_successful_kino_-1)/num_successful_kino_ +
                                 runtime_total/num_successful_kino_;
          runtime_length_kino_ = runtime_length_kino_*(num_successful_kino_-1)/num_successful_kino_ +
                                 runtime_length_check/num_successful_kino_;                               
          node_expanded_kino_ = node_expanded_kino_*(num_successful_kino_-1)/num_successful_kino_ +
                                 (double)node_used/num_successful_kino_;     
          log_msg_.twist.twist.linear.x = runtime_total;
          log_msg_.twist.twist.linear.y = runtime_length_check;
          log_msg_.twist.twist.linear.z = node_used;

          search_process_ = 1; //if one of the search is not successful, then we will send a new target soon                               
          // if (node_expanded_kino_*3.5 < node_expanded_astar_)
          // {
          //    //this means the path computed are likely different, executing this kino path will lead to infeasible
          //    // starting condition for astar
          //   search_process_ = -1; 
          // }
          double length_astar, length_kino;
          currentSamplePath = astar_path_finder_-> getPath(length_astar);
          pubSampledTraj(currentSamplePath);    
          log_msg_.pose.pose.orientation.x = length_astar;

          kino_ptr_->getPath(currentSamplePath_kino, length_kino);
          pubSampledTraj(currentSamplePath_kino, 2);
          log_msg_.twist.twist.angular.x = length_kino;

          // getchar();

        }
        gotTerminalGoal_=false;

      }
      
      if (search_process_<0) return;
    }
    /*-------------------------end of single agent benchmarking---------------------------------*/

    //record path for visualization
    if (!path_taken_.empty() && (state_.pos - path_taken_.back()).norm()>0.05)
    {
      path_taken_.push_back(state_.pos);
      if ((state_.pos - G_.pos).norm()>par_.goal_radius)
      {
        length_of_path_taken_+= 
          (path_taken_[path_taken_.size()-1]-path_taken_[path_taken_.size()-2]).norm();      
      }      
    }
    else if (path_taken_.empty()) path_taken_.push_back(state_.pos);



    std::vector<mt::state> safe_traj;

    std::vector<Hyperplane3D> planes;
    mt::PieceWisePol pwp;

    bool replanned = neptune_ptr_->replanFull(safe_traj, pwp);

    if (par_.visual && replanned)
    {
      // Delete markers to publish stuff
      visual_tools_->deleteAllMarkers();
      visual_tools_->enableBatchPublishing();

      pubTraj(safe_traj);
    }

    if (replanned)
    {
      publishOwnTraj(pwp);
      pwp_last_ = pwp;
    }
    else
    {
      if (timer_stop_.ElapsedMs() > 50.0)  // publish every half a second. TODO set as param
      {
        // std::cout << bold << "trying to publish pwp_last_ " << reset << std::endl;
        publishOwnTraj(pwp_last_);  // This is needed because is drone DRONE1 stops, it needs to keep publishing his
                                    // last planned trajectory, so that other drones can avoid it (even if DRONE1 was
                                    // very far from the other drones with it last successfully planned a trajectory).
                                    // Note that these trajectories are time-indexed, and the last position is taken if
                                    // t>times.back(). See eval() function in the pwp struct
        timer_stop_.Reset();
      }
    }
    bool kino_success, opt_success;
    double kino_compute_time, opt_compute_time, objective_value;
    neptune_ptr_->getPlanningStats(kino_success, kino_compute_time, 
                                 opt_success, opt_compute_time, objective_value);

    if (!single_benchmark_ && (use_moveback_ && timer_start_moveback_.ElapsedMs()<4000))
    {
      log_msg_.header.stamp = ros::Time::now();
      // log_msg_.header.frame_id = world_name_;
      log_msg_.pose.pose.position.x = kino_success? 1.0:0.0;
      log_msg_.pose.pose.position.y = kino_compute_time>2000?2000:kino_compute_time; //prevent large number when undefined
      log_msg_.pose.pose.position.z = opt_success? 1.0:0.0;    
      log_msg_.pose.pose.orientation.x = opt_compute_time>2000?2000:opt_compute_time;    
      log_msg_.pose.pose.orientation.y = objective_value;  
      log_msg_.pose.pose.orientation.y = just_got_new_goal?1.0:0.0;    
      pub_log_.publish(log_msg_);        
    }

  }
    // std::cout << "[Callback] Leaving replanCB" << std::endl;

  /*************DEBUGGING RUN KINO ONCE*************/
  // for (int i=0; i<1; i++) //if you want to run multiple times per iteration
  // {
  //   if (gotTerminalGoal_||search_process_==2 &&getchar()=='s')
  //   {
  //     gotTerminalGoal_=false;
  //     int current_process;
  //     std::vector<Eigen::Vector3d> currentSamplePath;
  //     kino_ptr_->runonce(currentSamplePath, current_process);
  //     pubSampledTraj(currentSamplePath);
  //     search_process_ = current_process;
  //   } 
  // }
  /**************END DEBUGGING RUN KINO ONCE************/

  //   // view sampled traj for other agents
  //   double t_start = ros::Time::now().toSec();
  //   double t_final = t_start + T_span_ * par_.num_pol;  
  //   mt::SampledPointsofCurves SampledPointsForAll = neptune_ptr_->SamplePointsOfCurves(t_start, t_final);
  //   pubSampledAgentsTraj(SampledPointsForAll);

}


void NeptuneRos::stateCB(const snapstack_msgs::State& msg)
{
  mt::state state_tmp;
  state_tmp.setPos(msg.pos.x, msg.pos.y, msg.pos.z);
  state_tmp.setVel(msg.vel.x, msg.vel.y, msg.vel.z);
  state_tmp.setAccel(0.0, 0.0, 0.0);
  // std::cout << bold << red << "MSG_QUAT= " << msg.quat << reset << std::endl;
  double roll, pitch, yaw;
  mu::quaternion2Euler(msg.quat, roll, pitch, yaw);
  state_tmp.setYaw(yaw);
  state_ = state_tmp;
  // std::cout << bold << red << "STATE_YAW= " << state_.yaw << reset << std::endl;
  // neptune_ptr_->updateState(state_tmp);

  if (published_initial_position_ == false)
  {
    pwp_last_ = mu::createPwpFromStaticPosition(state_);
    publishOwnTraj(pwp_last_);
    published_initial_position_ = true;
  }
  if (neptune_ptr_->IsTranslating() == true && par_.visual)
  {
    pubActualTraj();
  }

  if (previousCheckingPos_.empty())
  {
    setUpCheckingPosAndStaticObs(Eigen::Vector2d(msg.pos.x, msg.pos.y));
  }
  else
  {
    updateEntStateStaticObs(Eigen::Vector2d(msg.pos.x, msg.pos.y));
  }  
  // std::cout << bold << white << "dist to goal is = " << (state_.pos-G_.pos).norm() << reset << std::endl;

}

void NeptuneRos::odomCB(const nav_msgs::Odometry& msg)
{
  if ((ros::Time::now()-odomCBTime_).toSec()<0.02)
  {
    return;
  }
  odomCBTime_ = ros::Time::now();

  mt::state state_tmp;
  state_tmp.setPos(msg.pose.pose.position.x, msg.pose.pose.position.y, msg.pose.pose.position.z);

  Eigen::Quaterniond orientation_W_B(msg.pose.pose.orientation.w, msg.pose.pose.orientation.x, 
                                     msg.pose.pose.orientation.y, msg.pose.pose.orientation.z);
  Eigen::Vector3d velocity_body(msg.twist.twist.linear.x, msg.twist.twist.linear.y, msg.twist.twist.linear.z);

  Eigen::Vector3d velocity_world = velocity_body;//orientation_W_B *velocity_body;    
  if (use_vicon_odom_) velocity_world = orientation_W_B *velocity_body; 

  state_tmp.setVel(velocity_world(0), velocity_world(1), velocity_world(2));
  state_tmp.setAccel(0.0, 0.0, 0.0);
  // std::cout << bold << red << "MSG_QUAT= " << msg.quat << reset << std::endl;
  double roll, pitch, yaw;
  mu::quaternion2Euler(msg.pose.pose.orientation, roll, pitch, yaw);
  state_tmp.setYaw(yaw);
  state_ = state_tmp;
  // std::cout << bold << red << "STATE_YAW= " << state_.yaw << reset << std::endl;
  // neptune_ptr_->updateState(state_tmp);


  if (published_initial_position_ == false)
  {
    pwp_last_ = mu::createPwpFromStaticPosition(state_);
    publishOwnTraj(pwp_last_);
    published_initial_position_ = true;
  }
  if (neptune_ptr_->IsTranslating() == true && par_.visual)
  {
    pubActualTraj();
  }

  if (previousCheckingPos_.empty())
  {
    setUpCheckingPosAndStaticObs(Eigen::Vector2d(msg.pose.pose.position.x, msg.pose.pose.position.y));
  }
  else
  {
    updateEntStateStaticObs(Eigen::Vector2d(msg.pose.pose.position.x, msg.pose.pose.position.y));
  }
  
}

// this function not only updates the static obstacles, naming should be changed
void NeptuneRos::updateEntStateStaticObs(Eigen::Vector2d pkplus1)
{
  std::vector<Eigen::Vector2i> alphasToAdd;
  if (previousCheckingPos_.empty()) return;
  if ((previousCheckingPos_[num_of_agents_]-pkplus1).norm()<0.05 && 
      previousCheckingAgentTimer_[0].ElapsedMs() < 100) return;

  previousCheckingAgentTimer_[0].Reset();
  for (int i=0; i<num_of_agents_; i++) //update individual agent
  {
    if (i == id_-1) continue;
    //check agent msg not received yet
    if (previousCheckingPosAgent_[i](0)<-900 || bendPtsForAgents_[i].empty()) continue; 

    eu::entangleHSigToAddAgentInd(alphasToAdd, previousCheckingPos_[i], pkplus1,  
                              previousCheckingPosAgent_[i], latestCheckingPosAgent_[i], 
                              pb_[id_-1], bendPtsForAgents_[i], bendPtsForAgents_prev_[i], i+1);
    previousCheckingPos_[i] = pkplus1;
    previousCheckingPosAgent_[i] = latestCheckingPosAgent_[i];
    neptune_ptr_->setBendPtsForAgent(i+1, bendPtsForAgents_[i]);  
  }
  bendPtsForAgents_prev_ = bendPtsForAgents_;
  eu::entangleHSigToAddStatic(alphasToAdd, previousCheckingPos_[num_of_agents_], pkplus1,  
                              staticObsRep_, num_of_agents_);     

  eu::addAlphaBetaToList(alphasToAdd, entangle_state_, previousCheckingPos_[num_of_agents_], 
                         pb_, pb_[id_-1], staticObsRep_, num_of_agents_, bendPtsForAgents_);
  eu::updateBendPts(entangle_state_, pkplus1, pb_, pb_[id_-1], staticObsRep_, num_of_agents_);

  /*----------------------- below for original hsig update, comment if not in use---*/
  // eu::updateHSig_orig(entangle_state_orig_, previousCheckingPos_[num_of_agents_], 
  //                     pkplus1, staticObsRep_);

  // entangle_state_orig_.contPts.push_back(pkplus1);
  // eu::updateContPts_orig(entangle_state_orig_.contPts, convexHullOfStaticObs_);
  /*----------------------- above for original hsig update, comment if not in use---*/

  previousCheckingPos_[num_of_agents_] = pkplus1;

  // eu::shrinkAlphaBetas(entangle_state_.alphas, entangle_state_.betas, entangle_state_.active_cases, num_of_agents_);
  neptune_ptr_->setAlphasBetas(entangle_state_, previousCheckingPos_, previousCheckingPosAgent_);
  neptune_ptr_->updateState(state_);
  if (par_.enable_entangle_check) pubTetherShapeRough();

  for (int i=0; i<num_of_agents_; i++)
  {
    if (entangle_state_.active_cases[i]>2)
    {
      std::cout<<bold<<red<<"entangled with agent "<<i<<" !!!"<<std::endl;
      // exit(-1);
    }
  }
}

void NeptuneRos::setUpCheckingPosAndStaticObs(Eigen::Vector2d pos)
{
  for (int i=0; i<num_of_agents_ + 1; i++)
  {
    previousCheckingPos_.push_back(pos);
  }

  std::vector<Eigen::Vector2d> centers;
  for (int i=0; i<par_.num_of_static_obst; i++)
  {
    int num_of_vertices = convexHullOfStaticObs_[i].cols();
    Eigen::Vector2d center(0.0, 0.0);
    for (int j=0; j<num_of_vertices; j++)
    {
      center += convexHullOfStaticObs_[i].col(j);
    }
    center = center/num_of_vertices;
    centers.push_back(center);
  }

  double c11 = pb_[id_-1].x();
  double c22 = pb_[id_-1].y();
  double d11 = pos.x();
  double d22 = pos.y();  

  double theta = atan2((d22-c22), (d11-c11));

  bool setting_correct = false;
  while (!setting_correct && par_.enable_entangle_check)
  {
    if (theta > atan2((d22-c22), (d11-c11)) + 3.14)
    {
      std::cout<<bold<<red<<"cannot find a feasible vertex representation for all obsts!"
                          <<" exitting!"<<std::endl;
      exit(-1);  
    }
    // std::cout<<bold<<blue<<"theta for this round is "<<theta/3.14159*180<<std::endl;
    staticObsRep_.clear();
    staticObsLongestDist_.clear();

    double m = tan(theta);
    for (int i=0; i<par_.num_of_static_obst; i++)
    {
      double e1 = centers[i](1) - centers[i](0)*m; //y = mx +e1
      for (int j = i+1; j<par_.num_of_static_obst; j++)
      {
        double e2 = centers[j](1) - centers[j](0)*m; //y = mx +e2
        double distance = fabs(e1-e2)/sqrt(1+m*m);
        // std::cout<<bold<<blue<<"distance is "<<distance<<std::endl;
        if (distance < par_.a_star_fraction_voxel_size*1.6) goto nextround;
      }

      double e2 = pb_[id_-1].y() - pb_[id_-1].x()*m; //check also with the basepoint
      double distance = fabs(e1-e2)/sqrt(1+m*m);
      if (distance < par_.a_star_fraction_voxel_size*1.6) goto nextround;

      double a1 = centers[i](0);
      double a2 = centers[i](1);
      double b1 = centers[i](0) + 10;
      double b2 = centers[i](1) + 10*m;
      double c1 = pb_[id_-1].x();
      double c2 = pb_[id_-1].y();
      double d1 = pos.x();
      double d2 = pos.y();
      double _y = ((c2-a2)*(b1-a1) - (b2-a2)*(c1-a1)) / ((b2-a2)*(d1-c1)-(b1-a1)*(d2-c2));

      //check for interaction with the cable line
      if (fabs((b2-a2)*(d1-c1)-(b1-a1)*(d2-c2))<0.01 || _y<0 || _y>1.0) //no intersection
      {
        //do nothing
      }
      else //fail the check
      {
        goto nextround;
      }      
    }

    for (int i=0; i<par_.num_of_static_obst; i++)
    {
      int num_of_vertices = convexHullOfStaticObs_[i].cols();
      double _x_max = -1e-5;
      double _x_min = 1e-5;
      Eigen::Vector2d max_vert, min_vert;      
      for (int j=0; j<num_of_vertices; j++)
      {
        Eigen::Vector2d vertex1 = convexHullOfStaticObs_[i].col(j);
        Eigen::Vector2d vertex2 = convexHullOfStaticObs_[i].col((j+1)%num_of_vertices);
        double a1 = centers[i](0);
        double a2 = centers[i](1);
        double b1 = centers[i](0) + 10;
        double b2 = centers[i](1) + 10*m;
        double c1 = vertex1.x();
        double c2 = vertex1.y();
        double d1 = vertex2.x();
        double d2 = vertex2.y();
        double _y = ((c2-a2)*(b1-a1) - (b2-a2)*(c1-a1)) / ((b2-a2)*(d1-c1)-(b1-a1)*(d2-c2));   
        double _x; 


        //check for interaction with the cable line
        if (fabs((b2-a2)*(d1-c1)-(b1-a1)*(d2-c2))<0.01 || _y<0 || _y>1.0) 
        {
          //do nothing
        }
        else //fail the check
        {
          if (fabs(b1-a1)<0.01)  _x = (c2-a2)/(b2-a2) + _y * (d2-c2)/(b2-a2);
          else _x = (c1-a1)/(b1-a1) + _y * (d1-c1)/(b1-a1);          
          if (_x>_x_max)
          {
            _x_max = _x;
            max_vert << c1+_y*(d1-c1), c2+_y*(d2-c2);
          }
          else if (_x<_x_min)
          {
            _x_min = _x;
            min_vert << c1+_y*(d1-c1), c2+_y*(d2-c2);
          }
        }
        // if (i==5 && j==3)  //just for debugging 
        // {
        //   std::cout<<bold<<blue<<"_y is "<<_y<<std::endl;
        //   std::cout<<bold<<blue<<"_x is "<<_x<<std::endl;
        //   std::cout<<bold<<blue<<"_x_min is "<<_x_min<<std::endl;
        //   std::cout<<bold<<blue<<"_x_max is "<<_x_max<<std::endl;
        //   std::cout<<bold<<blue<<"min_vert is "<<min_vert<<std::endl;
        //   std::cout<<bold<<blue<<"max_vert is "<<max_vert<<std::endl;
        // } 
      }
      if (_x_max<0 ||_x_min>0) goto nextround;
      Eigen::Matrix<double, 2, 2> pts;
      pts << min_vert, max_vert;       
      staticObsRep_.push_back(pts);     

      double maxl1 = 0, maxl2 = 0;
      int maxidx1 = 0, maxidx2 = 0;
      for (int j=0; j<num_of_vertices; j++)
      {
        if ((convexHullOfStaticObs_[i].col(j)-min_vert).norm() > maxl1)
        {
          maxl1 = (convexHullOfStaticObs_[i].col(j)-min_vert).norm();
          maxidx1 = j;
        }
        if ((convexHullOfStaticObs_[i].col(j)-max_vert).norm() > maxl2)
        {
          maxl2 = (convexHullOfStaticObs_[i].col(j)-max_vert).norm();
          maxidx2 = j;
        }        
      }
      Eigen::Vector2d longestdist(maxl1, maxl2); 
      staticObsLongestDist_.push_back(longestdist);
    }
    setting_correct = true;

    nextround:
    {
      theta = theta + 1.0/180.0 * 3.1415927;
    }
  }
  neptune_ptr_->updateState(state_);
  neptune_ptr_ -> setStaticObstRep(staticObsRep_, staticObsLongestDist_);
  G_prev_ = state_; //for astar planning
  
  astar_path_finder_ -> setStaticObstRep(staticObsRep_);
  kino_ptr_ -> setStaticObstRep(staticObsRep_, staticObsLongestDist_); //comment out when not using

  std::cout<<bold<<green<<"Finished setting up static obstacles representation! "<<std::endl;      
}



void NeptuneRos::modeCB(const mader_msgs::Mode& msg)
{
  // neptune_ptr_->changeMode(msg.mode);

  if (msg.mode != msg.GO)
  {  // MADER DOES NOTHING
    // sub_state_.shutdown();
    pubCBTimer_.stop();
    replanCBTimer_.stop();
    // std::cout << on_blue << "**************stopping replanCBTimer" << reset << std::endl;
    autoCMDTimer_.stop();

    neptune_ptr_->resetInitialization();
  }
  else
  {  // The mode changed to GO (the mode changes to go when takeoff is finished)
    // sub_state_ = nh_.subscribe("state", 1, &NeptuneRos::stateCB, this);  // TODO duplicated from above
    // std::cout << on_blue << "**************starting replanCBTimer" << reset << std::endl;
    pubCBTimer_.start();
    replanCBTimer_.start();
    autoCMDTimer_.start();
  }
}

void NeptuneRos::autoCMD(const ros::TimerEvent& e)
{
  if (!auto_command_ || timer_auto_cmd_.ElapsedMs()<5000) return;
  if (num_successful_kino_ == 100) return;
  geometry_msgs::PoseStamped msg;
  msg.header.stamp = ros::Time::now();
  msg.pose.position.z = 0.05;

  if ((state_.pos - G_.pos).norm()<par_.goal_radius || !autoCMD_initialized_ || 
       timer_auto_cmd_.ElapsedMs()>45000 ||
       (timer_auto_cmd_.ElapsedMs()>3000 && gotTerminalGoal_==false && search_process_<0))
  {
    if (timer_auto_cmd_.ElapsedMs()<45000 &&
        state_.vel.head<2>().norm()>0.1 || state_.accel.head<2>().norm()>0.1) return;

    autoCMD_initialized_ = true;
    double safe_di = 3* par_.drone_radius;
    std::uniform_real_distribution<double> unif_x(par_.x_min+safe_di, par_.x_max-safe_di);
    std::uniform_real_distribution<double> unif_y(par_.y_min+safe_di, par_.y_max-safe_di);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re = std::default_random_engine(seed);
    bool goal_valid = false;
    while (!goal_valid)
    {
      double ran_x = unif_x(re);
      double ran_y = unif_y(re);
      // std::cout<<bold<<red<<"ranx is "<<ran_x<<" , rany is "<<ran_y<<std::endl;
      if ((Eigen::Vector2d(ran_x, ran_y)-state_.pos.head(2)).norm()<5.0)
        goto nextround;
      if ((Eigen::Vector2d(ran_x, ran_y)-pb_[id_-1]).norm()>par_.tetherLength*0.85)
        goto nextround;
      for (int i=0; i<par_.num_of_static_obst; i++)
      {
        Point_2 points[inflatedStaticObsLarger_[i].cols()];
        for (int j=0; j<inflatedStaticObsLarger_[i].cols(); j++)
        {
          points[j] = Point_2(inflatedStaticObsLarger_[i](0,j), inflatedStaticObsLarger_[i](1,j));
        }

        if (cu::check_inside(Point_2(ran_x, ran_y), points, 
                                     points+inflatedStaticObsLarger_[i].cols(), K()))
        {
          goal_valid = false;
          goto nextround;
        }        
      }
      msg.pose.position.x = ran_x;
      msg.pose.position.y = ran_y;

      goal_valid = true;
      nextround:;
    }
    terminalGoalCB(msg);
    timer_auto_cmd_.Reset();
  }
}

void NeptuneRos::reactiveController(bool& controller_enabled, Eigen::Vector2d& tarPos, int& itr)
{
  if (itr>100)
  {
    controller_enabled = false;
    return;
  }
  for (int i=0; i<num_of_agents_; i++)
  {
    if (i==id_-1) continue;
    if (latestCheckingPosAgent_[i](0)<par_.x_min) continue;

    if ((tarPos-latestCheckingPosAgent_[i]).norm()<2*par_.drone_radius-0.01)
    {
      controller_enabled = true;
      itr ++;
      // std::cout<<bold<<red<<"itr is" <<itr<<std::endl;
      // std::cout<<bold<<red<<"tarPos original is" <<tarPos<<std::endl; 
      tarPos = latestCheckingPosAgent_[i] + (tarPos-latestCheckingPosAgent_[i]) *
              2.2*par_.drone_radius/(tarPos-latestCheckingPosAgent_[i]).norm();

      // std::cout<<bold<<red<<"latestCheckingPosAgent_[i] is" <<latestCheckingPosAgent_[i]<<std::endl;  
      // std::cout<<bold<<red<<"tarPos new is" <<tarPos<<std::endl;        
      reactiveController(controller_enabled, tarPos, itr);
      return;
    }
  }  
}


void NeptuneRos::pubCB(const ros::TimerEvent& e)
{
  mt::state next_goal;
  bool last_point;

  if (neptune_ptr_->getNextGoal(next_goal, last_point))
  {

    snapstack_msgs::Goal quadGoal;

    quadGoal.p = mu::eigen2rosvector(next_goal.pos);
    quadGoal.v = mu::eigen2rosvector(next_goal.vel);
    quadGoal.a = mu::eigen2rosvector(next_goal.accel);
    quadGoal.j = mu::eigen2rosvector(next_goal.jerk);
    quadGoal.dyaw = next_goal.dyaw;
    quadGoal.yaw = next_goal.yaw;
    quadGoal.header.stamp = ros::Time::now();
    quadGoal.header.frame_id = world_name_;

    pub_goal_.publish(quadGoal);

    setpoint_.header.stamp = ros::Time::now();
    setpoint_.pose.position.x = quadGoal.p.x;
    setpoint_.pose.position.y = quadGoal.p.y;
    setpoint_.pose.position.z = quadGoal.p.z;

    pub_setpoint_.publish(setpoint_);

    bool contEnabled = false;
    Eigen::Vector2d tarPos = state_.pos.head(2);
    int _itr = 0;
    reactiveController(contEnabled, tarPos, _itr);

    if (contEnabled)
    {
      trajectory_msgs::MultiDOFJointTrajectory trajset_msg;
      trajset_msg.header.stamp = ros::Time::now();
      trajectory_msgs::MultiDOFJointTrajectoryPoint trajpt_msg;
      geometry_msgs::Transform transform_msg;
      geometry_msgs::Twist accel_msg, vel_msg;
      transform_msg.translation.x = tarPos(0);
      transform_msg.translation.y = tarPos(1);
      transform_msg.translation.z = next_goal.pos(2);
      transform_msg.rotation.x = 0;
      transform_msg.rotation.y = 0;
      transform_msg.rotation.z = sinf(next_goal.yaw*0.5);
      transform_msg.rotation.w = cosf(next_goal.yaw*0.5);
      trajpt_msg.transforms.push_back(transform_msg);
      vel_msg.linear.x = 0.0;
      vel_msg.linear.y = 0.0;
      vel_msg.linear.z = 0.0;
      accel_msg.linear.x = 0.0;
      accel_msg.linear.y = 0.0;
      accel_msg.linear.z = 0.0;
      trajpt_msg.velocities.push_back(vel_msg);
      trajpt_msg.accelerations.push_back(accel_msg);
      trajset_msg.points.push_back(trajpt_msg);    
      pub_cmd_.publish(trajset_msg);   
    }
    else if (last_point && 
            (next_goal.pos.head(2)-G_.pos.head(2)).norm()<3.0*par_.goal_radius)
    {
      //publish for rotors simulation
      trajectory_msgs::MultiDOFJointTrajectory trajset_msg;
      trajset_msg.header.stamp = ros::Time::now();
      trajectory_msgs::MultiDOFJointTrajectoryPoint trajpt_msg;
      geometry_msgs::Transform transform_msg;
      geometry_msgs::Twist accel_msg, vel_msg;
      transform_msg.translation.x = G_.pos(0);
      transform_msg.translation.y = G_.pos(1);
      transform_msg.translation.z = next_goal.pos(2);
      transform_msg.rotation.x = 0;
      transform_msg.rotation.y = 0;
      transform_msg.rotation.z = sinf(next_goal.yaw*0.5);
      transform_msg.rotation.w = cosf(next_goal.yaw*0.5);
      trajpt_msg.transforms.push_back(transform_msg);
      vel_msg.linear.x = 0.0;
      vel_msg.linear.y = 0.0;
      vel_msg.linear.z = 0.0;
      accel_msg.linear.x = 0.0;
      accel_msg.linear.y = 0.0;
      accel_msg.linear.z = 0.0;
      trajpt_msg.velocities.push_back(vel_msg);
      trajpt_msg.accelerations.push_back(accel_msg);
      trajset_msg.points.push_back(trajpt_msg);    
      pub_cmd_.publish(trajset_msg);   
    }
    else
    {
      //publish for rotors simulation
      trajectory_msgs::MultiDOFJointTrajectory trajset_msg;
      trajset_msg.header.stamp = ros::Time::now();
      trajectory_msgs::MultiDOFJointTrajectoryPoint trajpt_msg;
      geometry_msgs::Transform transform_msg;
      geometry_msgs::Twist accel_msg, vel_msg;
      transform_msg.translation = quadGoal.p;
      transform_msg.rotation.x = 0;
      transform_msg.rotation.y = 0;
      transform_msg.rotation.z = sinf(next_goal.yaw*0.5);
      transform_msg.rotation.w = cosf(next_goal.yaw*0.5);
      trajpt_msg.transforms.push_back(transform_msg);
      vel_msg.linear = quadGoal.v;
      accel_msg.linear = quadGoal.a;
      trajpt_msg.velocities.push_back(vel_msg);
      trajpt_msg.accelerations.push_back(accel_msg);
      trajset_msg.points.push_back(trajpt_msg);    
      pub_cmd_.publish(trajset_msg);        
    }
  }

  if (timer_pub_.ElapsedMs() > 100.0)
  {
    pubSampledTraj(path_taken_, 3);

    if (par_.enable_entangle_check)
    {
      for (int i=0; i<entangle_state_.alphas.size(); i++)
      {
        std::cout<<bold<<blue<<"alphas "<<i<<"th entry is" <<entangle_state_.alphas[i]<<std::endl;
      }
      for (int i=0; i<entangle_state_.active_cases.size(); i++)
      {
        // std::cout<<bold<<blue<<"active cases for "<<i<<"th agent is" <<entangle_state_.active_cases[i]<<std::endl;
      }
      // for (int i=0; i<entangle_state_.bendPointsIdx.size(); i++)
      // {
      //   std::cout<<bold<<blue<<i<<"th bendPointsIdx is " <<
      //             entangle_state_.alphas[entangle_state_.bendPointsIdx[i]](0)<<std::endl;
      // }

      Eigen::Vector2d pkplus1(state_.pos(0), state_.pos(1));
      double length = eu::getTetherLength(entangle_state_, pb_, pb_[id_-1], pkplus1, 
                                          staticObsRep_, staticObsLongestDist_, num_of_agents_);
      std::cout<<bold<<blue<<"tether length is "<<length<<std::endl;
      // std::cout<<bold<<blue<<"timer_auto_cmd_ is "<<timer_auto_cmd_.ElapsedMs()<<std::endl;

      for (int i=0; i<entangle_state_orig_.hsig.size(); i++)
      {
        // std::cout<<bold<<blue<<"[hsig] "<<i<<"th entry is" <<entangle_state_orig_.hsig[i]<<std::endl;
      }      
    }
  
    // std::cout<<bold<<green<<"[homotopic Astar] Avg runtime is "<<runtime_total_astar_
    // <<", Avg time spent for length check is "<<runtime_length_astar_
    // <<", avg node expanded is "<<node_expanded_astar_<<std::endl;
    // std::cout<<bold<<green<<"[kino-homo Astar] Avg runtime is "<<runtime_total_kino_
    // <<", Avg time spent for length check is "<<runtime_length_kino_
    // <<", avg node expanded is "<<node_expanded_kino_<<std::endl;
    // std::cout<<bold<<green<<"[homotopic Astar] successful run is "<<num_successful_astar_<<
    //                         "/"<<num_targets_sent_<<std::endl;
    // std::cout<<bold<<green<<"[kino-homo Astar] successful run is "<<num_successful_kino_<<
    //                         "/"<<num_targets_sent_<<std::endl; 

    // pubPathMeshPlane();
    timer_pub_.Reset();
  }

  if (on_way_back_)
  {
    if ((state_.pos.head(2)-pb_[id_-1]).norm()<6.33-(par_.num_of_agents-4)*0.33 
      || timer_start_moveback_.ElapsedMs()>10000)
    {
      neptune_ptr_->setTerminalGoal(G_);  
      mt::state G;  // projected goal
      neptune_ptr_->getG(G);

      pubState(G_, pub_point_G_term_);
      pubState(G, pub_point_G_);
      on_way_back_ = false;

    }
  }
}

void NeptuneRos::clearMarkerArray(visualization_msgs::MarkerArray* tmp, ros::Publisher* publisher)
{
  if ((*tmp).markers.size() == 0)
  {
    return;
  }
  int id_begin = (*tmp).markers[0].id;

  for (int i = 0; i < (*tmp).markers.size(); i++)
  {
    visualization_msgs::Marker m;
    m.header.frame_id = "world";
    m.header.stamp = ros::Time::now();
    m.type = visualization_msgs::Marker::ARROW;
    m.action = visualization_msgs::Marker::DELETE;
    m.id = i + id_begin;
    m.scale.x = 0.15;
    m.scale.y = 0;
    m.scale.z = 0;
    (*tmp).markers[i] = m;
  }

  (*publisher).publish(*tmp);
  (*tmp).markers.clear();
}

void NeptuneRos::pubTraj(const std::vector<mt::state>& data)
{
  // Trajectory
  nav_msgs::Path traj;
  traj.poses.clear();
  traj.header.stamp = ros::Time::now();
  traj.header.frame_id = world_name_;

  geometry_msgs::PoseStamped temp_path;

  int increm = (int)std::max(data.size() / par_.res_plot_traj, 1.0);  // this is to speed up rviz

  for (int i = 0; i < data.size(); i = i + increm)
  {
    temp_path.pose.position.x = data[i].pos(0);
    temp_path.pose.position.y = data[i].pos(1);
    temp_path.pose.position.z = data[i].pos(2);
    temp_path.pose.orientation.w = 1;
    temp_path.pose.orientation.x = 0;
    temp_path.pose.orientation.y = 0;
    temp_path.pose.orientation.z = 0;
    traj.poses.push_back(temp_path);
  }

  pub_traj_safe_.publish(traj);
  clearMarkerArray(&traj_safe_colored_, &pub_traj_safe_colored_);

  double scale = 0.15;

  traj_safe_colored_ = mu::trajectory2ColoredMarkerArray(data, par_.v_max.maxCoeff(), increm, name_drone_, scale,
                                                         par_.color_type, id_, par_.num_of_agents);
  pub_traj_safe_colored_.publish(traj_safe_colored_);
}

void NeptuneRos::pubSampledTraj(const std::vector<Eigen::Vector3d>& path, int index)
{
  nav_msgs::Path traj;
  traj.poses.clear();
  traj.header.stamp = ros::Time::now();
  traj.header.frame_id = world_name_;

  geometry_msgs::PoseStamped temp_path;

  for (int i=0; i<path.size(); i++)
  {
    temp_path.pose.position.x = path[i](0);
    temp_path.pose.position.y = path[i](1);
    temp_path.pose.position.z = path[i](2);
    temp_path.pose.orientation.w = 1;
    temp_path.pose.orientation.x = 0;
    temp_path.pose.orientation.y = 0;
    temp_path.pose.orientation.z = 0;
    traj.poses.push_back(temp_path);
  }

  if (index==1) pub_sampled_traj_.publish(traj);
  else if (index==2) pub_sampled_traj2_.publish(traj);
  else if (index==3) pub_sampled_traj3_.publish(traj);
  else if (index==4) pub_sampled_traj4_.publish(traj);
}

void NeptuneRos::pubActualTraj()
{
  static geometry_msgs::Point p_last = mu::pointOrigin();

  mt::state current_state;
  neptune_ptr_->getState(current_state);
  Eigen::Vector3d act_pos = current_state.pos;

  visualization_msgs::Marker m;
  m.type = visualization_msgs::Marker::ARROW;
  m.action = visualization_msgs::Marker::ADD;
  m.id = actual_trajID_;  // % 3000;  // Start the id again after ___ points published (if not RVIZ goes very slow)
  m.ns = "ActualTraj_" + name_drone_;
  actual_trajID_++;
  // m.color = mu::getColorJet(current_state.vel.norm(), 0, par_.v_max.maxCoeff());  // mu::color(mu::red_normal);

  if (par_.color_type == "vel")
  {
    m.color = mu::getColorJet(current_state.vel.norm(), 0, par_.v_max.maxCoeff());  // note that par_.v_max is per axis!
  }
  else
  {
    m.color = mu::getColorJet(id_, 0, par_.num_of_agents);  // note that par_.v_max is per axis!
  }

  m.scale.x = 0.15;
  m.scale.y = 0.0001;
  m.scale.z = 0.0001;
  m.header.stamp = ros::Time::now();
  m.header.frame_id = world_name_;

  // pose is actually not used in the marker, but if not RVIZ complains about the quaternion
  m.pose.position = mu::pointOrigin();
  m.pose.orientation.x = 0.0;
  m.pose.orientation.y = 0.0;
  m.pose.orientation.z = 0.0;
  m.pose.orientation.w = 1.0;

  geometry_msgs::Point p;
  p = mu::eigen2point(act_pos);
  m.points.push_back(p_last);
  m.points.push_back(p);
  p_last = p;

  if (m.id == 0)
  {
    return;
  }

  pub_actual_traj_.publish(m);
}

void NeptuneRos::addPointToMarker(visualization_msgs::Marker &m, const Eigen::Vector2d& xy, double z)
{
  geometry_msgs::Point p;
  p.x = xy.x();
  p.y = xy.y();
  p.z = z;
  m.points.push_back(p);  
}


void NeptuneRos::pubTetherShapeRough()
{

  visualization_msgs::Marker m;
  m.type = visualization_msgs::Marker::LINE_LIST;
  m.action = visualization_msgs::Marker::DELETEALL;
  m.id = 1;
  m.scale.x = 0.15;
  m.scale.y = 0;
  m.scale.z = 0;
  m.ns = "Tether_" + name_drone_;

  pub_tether_line_.publish(m);

  m.action = visualization_msgs::Marker::ADD;
  // m.color = mu::getColorJet(current_state.vel.norm(), 0, par_.v_max.maxCoeff());  // mu::color(mu::red_normal);

  m.color = mu::getColorJetInt(id_, 0, num_of_agents_);  // note that par_.v_max is per axis!
  m.header.stamp = ros::Time::now();
  m.header.frame_id = world_name_;

  // pose is actually not used in the marker, but if not RVIZ complains about the quaternion
  m.pose.position = mu::pointOrigin();
  m.pose.orientation.x = 0.0;
  m.pose.orientation.y = 0.0;
  m.pose.orientation.z = 0.0;
  m.pose.orientation.w = 1.0;

  Eigen::Vector2d pkplus1(state_.pos(0), state_.pos(1));

  addPointToMarker(m, pb_[id_-1], 0.0); 
  for (int i=0; i<entangle_state_.bendPointsIdx.size(); i++)
  {
    Eigen::Vector2d bendpt_inter;
    double not_used;
    eu::getBendPt2dwIdx(bendpt_inter, entangle_state_, pb_, staticObsRep_, staticObsLongestDist_,
                        num_of_agents_, entangle_state_.bendPointsIdx[i], not_used);
    addPointToMarker(m, bendpt_inter, 0.0); 
    addPointToMarker(m, bendpt_inter, 0.0); 

  }
  // for (int i=1; i<entangle_state_orig_.contPts.size()-1; i++)
  // {
  //   Eigen::Vector2d bendpt_inter = entangle_state_orig_.contPts[i];
  //   addPointToMarker(m, bendpt_inter, 0.0); 
  //   addPointToMarker(m, bendpt_inter, 0.0); 

  // }  
  addPointToMarker(m, pkplus1, state_.pos(2));  

  if (m.id == 0) return;
  // std::cout<<bold<<green<<"PUBLISHING "<<std::endl;
  // visualization_msgs::Marker m_copy = m;
  // m_copy.action = visualization_msgs::Marker::DELETEALL;
  // clearMarkerArray(&m_copy, &pub_tether_line_);
  pub_tether_line_.publish(m);  
}


void NeptuneRos::pubSampledAgentsTraj(mt::SampledPointsofCurves input)
{

  visualization_msgs::Marker m;
  m.type = visualization_msgs::Marker::LINE_LIST;
  m.action = visualization_msgs::Marker::ADD;
  m.id = actual_trajID_;  // % 3000;  // Start the id again after ___ points published (if not RVIZ goes very slow)
  m.ns = "ActualTraj_" + name_drone_;
  actual_trajID_++;
  // m.color = mu::getColorJet(current_state.vel.norm(), 0, par_.v_max.maxCoeff());  // mu::color(mu::red_normal);

  m.color = mu::getColorJet(id_, 0, par_.num_of_agents);  // note that par_.v_max is per axis!

  m.scale.x = 0.15;
  m.scale.y = 0.0001;
  m.scale.z = 0.0001;
  m.header.stamp = ros::Time::now();
  m.header.frame_id = world_name_;

  // pose is actually not used in the marker, but if not RVIZ complains about the quaternion
  m.pose.position = mu::pointOrigin();
  m.pose.orientation.x = 0.0;
  m.pose.orientation.y = 0.0;
  m.pose.orientation.z = 0.0;
  m.pose.orientation.w = 1.0;

  for (int i=0; i<input.size(); i++)
  {
    for (int j=0; j<input[i].size(); j++)
    {
      for (int k=0; k<input[i][j].cols()-1; k++)
      {
        // std::cout<<bold<<green<<"i j k : "<<i<<j<<k<<std::endl;
        // std::cout<<bold<<green<<"column k is : "<<input[i][j].col(k)<<std::endl;
        geometry_msgs::Point p;
        p.x = input[i][j](0,k);
        p.y = input[i][j](1,k);
        p.z = 0.0;
        // std::cout<<bold<<green<<"pushing back 1 "<<i<<j<<k<<std::endl;

        m.points.push_back(p);
        p.x = input[i][j](0,k+1);
        p.y = input[i][j](1,k+1);
        // std::cout<<bold<<green<<"pushing back 2 "<<i<<j<<k<<std::endl;

        m.points.push_back(p);       
        // std::cout<<bold<<green<<"after pushing back "<<std::endl;

      }
    }
  }

  if (m.id == 0)
  {
    return;
  }
  // std::cout<<bold<<green<<"PUBLISHING "<<std::endl;

  pub_sampled_agents_traj_.publish(m);
}

void NeptuneRos::clearMarkerActualTraj()
{
  visualization_msgs::Marker m;
  m.type = visualization_msgs::Marker::ARROW;
  m.action = visualization_msgs::Marker::DELETEALL;
  m.id = 0;
  m.scale.x = 0.02;
  m.scale.y = 0.04;
  m.scale.z = 1;
  pub_actual_traj_.publish(m);
  actual_trajID_ = 0;
}

void NeptuneRos::clearMarkerColoredTraj()
{
  visualization_msgs::Marker m;
  m.type = visualization_msgs::Marker::ARROW;
  m.action = visualization_msgs::Marker::DELETEALL;
  m.id = 0;
  m.scale.x = 1;
  m.scale.y = 1;
  m.scale.z = 1;
  pub_actual_traj_.publish(m);
}

void NeptuneRos::pubState(const mt::state& data, const ros::Publisher pub)
{
  geometry_msgs::PointStamped p;
  p.header.stamp = ros::Time::now();
  p.header.frame_id = world_name_;
  p.point = mu::eigen2point(data.pos);
  pub.publish(p);
}

void NeptuneRos::terminalGoalCB(const geometry_msgs::PoseStamped& msg)
{
  mt::state G_term;
  double z;
  if (fabs(msg.pose.position.z) < 1e-5)  // This happens when you click in RVIZ (msg.z is 0.0)
  {
    z = 1.0;
  }
  else  // This happens when you publish by yourself the goal (should always be above the ground)
  {
    z = msg.pose.position.z;
  }
  G_term.setPos(msg.pose.position.x, msg.pose.position.y, z);
  std::cout<<bold<<red<<"getting g_term pos as "<< G_term.pos<<std::endl;

  if (!single_benchmark_ && use_moveback_)
  {
    mt::state near_base;
    near_base.setPos(pb_[id_-1].x(), pb_[id_-1].y(), z);    
    neptune_ptr_->setTerminalGoal(near_base);   
    on_way_back_ = true; 
    pubState(G_term, pub_point_G_term_);
    pubState(G_term, pub_point_G_);
    timer_start_moveback_.Reset();

  }
  else
  {
    neptune_ptr_->setTerminalGoal(G_term);    
    mt::state G;  // projected goal
    neptune_ptr_->getG(G);

    pubState(G_term, pub_point_G_term_);
    pubState(G, pub_point_G_);
  }



  /*****************set up kino_ptr for testing***************/
  if (single_benchmark_)
  {
    double t_start = ros::Time::now().toSec();
    double t_final = t_start + T_span_ * par_.num_pol;
    mt::SampledPointsofCurves SampledPointsForAll = neptune_ptr_->SamplePointsOfCurves(t_start, t_final);
    mt::ConvexHullsOfCurves_Std2d tmp_hull;
    mt::ConvexHullsOfCurves_Std2d hulls2d = neptune_ptr_->convexHullsOfCurves2d( t_start, t_final, tmp_hull);
    Eigen::Vector3d goalpt(msg.pose.position.x, msg.pose.position.y, z);

    std::vector<Eigen::Matrix<double, 4, 1>> coeff_z_init;
    neptune_ptr_->getInitialZPwp(state_, z, coeff_z_init);
    kino_ptr_->setInitZCoeffs(coeff_z_init);

    kino_ptr_->setUp(state_, goalpt, hulls2d, SampledPointsForAll, entangle_state_, bendPtsForAgents_);
    kino_ptr_->clearProcess();
    pubSampledAgentsTraj(SampledPointsForAll);

    if (length_of_path_taken_>0.1 && search_process_==1)
    {
      astar_path_finder_->getHsigTerminal(entangle_state_orig_);
      G_prev_ = G_;
      log_msg_.twist.twist.angular.y = length_of_path_taken_;
      log_msg_.header.stamp = ros::Time::now();
      log_msg_.header.frame_id = world_name_;
      pub_log_.publish(log_msg_);
    }    
  }

  /*****************end for set up kino_ptr for testing***************/

  G_ = G_term;
  gotTerminalGoal_ = true;
  just_got_new_goal_ = true; //for publish compute data for first compute
  search_process_ = -1;
  num_targets_sent_++;

  length_of_path_taken_ = 0.0;
  path_taken_.clear();
  clearMarkerActualTraj();
}


void NeptuneRos::visualizeStaticObs(std::vector<mt::Polygon_Std>& convexHullOfStaticObs)
{
  decomp_ros_msgs::PolyhedronArray polyarray_msg;
  polyarray_msg.header.stamp = ros::Time::now();
  polyarray_msg.header.frame_id = world_name_;

  int num_of_static_obst = convexHullOfStaticObs.size();

  for (int i=0; i<num_of_static_obst; i++)
  {
    int num_of_vertices = convexHullOfStaticObs[i].cols();
    decomp_ros_msgs::Polyhedron poly;
    Eigen::Vector2d center(0.0, 0.0);
    for (int j=0; j<num_of_vertices; j++)
    {
      center += convexHullOfStaticObs[i].col(j);
    }
    center = center/num_of_vertices;

    for (int j=0; j<num_of_vertices; j++)
    {    
      Eigen::Vector2d plane_dir = convexHullOfStaticObs[i].col((j+1)%num_of_vertices) - 
                                  convexHullOfStaticObs[i].col(j);
      Eigen::Vector2d normal;
      if (fabs(plane_dir(0)) > 0.01)
      {
        normal(1) = 1.0;
        normal(0) = -plane_dir(1)/plane_dir(0);
      } 
      else if (fabs(plane_dir(1)) > 0.01)
      {
        normal(0) = 1.0;
        normal(1) = -plane_dir(0)/plane_dir(1);
      }
      else continue; //skip this is likely to be the same point

      if((center - convexHullOfStaticObs[i].col(j)).dot(normal) >= 0)
      {
          normal = -1.0 * normal;
      }
      double L = normal.norm();
      normal   = normal / L;


      geometry_msgs::Point pt, n;
      pt.x = convexHullOfStaticObs[i].col(j).x();
      pt.y = convexHullOfStaticObs[i].col(j).y();
      pt.z = 0.0;
      n.x = normal(0);
      n.y = normal(1);
      n.z = 0.0;
      poly.points.push_back(pt);
      poly.normals.push_back(n);    
    }
    // add upper and lower cap
    geometry_msgs::Point pt;
    pt.x = convexHullOfStaticObs[i].col(0).x();
    pt.y = convexHullOfStaticObs[i].col(0).y();
    pt.z = 0.0;
    poly.points.push_back(pt);

    pt.z = 2.0;
    poly.points.push_back(pt);

    pt.x = 0.0;
    pt.y = 0.0;
    pt.z = -1.0;    
    poly.normals.push_back(pt);

    pt.z = 1.0;        
    poly.normals.push_back(pt);

    polyarray_msg.polyhedrons.push_back(poly);
  }

  poly_safe_pub_.publish(polyarray_msg);
  std::cout << bold << "static obst visualization published!" << reset << std::endl;
}

