/* ----------------------------------------------------------------------------
 * Nanyang Technological University
 * Authors: Cao Muqing, et al.
 * Acknowledgement: Jesus Tordesillas
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <Eigen/StdVector>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stdlib.h>

#include "neptune.hpp"
#include "timer.hpp"
#include "termcolor.hpp"

# include "gjk.hpp"
#include <ddp_optimizer.h>

using namespace termcolor;

// Uncomment the type of timer you want:
typedef MADER_timers::ROSTimer MyROSTimer;
// typedef ROSWallTimer MyTimer;
typedef MADER_timers::Timer MyTimer;

Neptune::Neptune(mt::parameters par) : par_(par)
{
  drone_status_ == DroneStatus::YAWING;
  G_.pos << 0, 0, 0;
  G_term_.pos << 0, 0, 0;

  changeDroneStatus(DroneStatus::GOAL_REACHED);
  resetInitialization();

  double sample_T=par_.T_span;

  mt::basisConverter basis_converter;

  if (par.basis == "MINVO")
  {
    A_rest_pos_basis_ = basis_converter.getArestMinvo();
  }
  else if (par.basis == "BEZIER")
  {
    A_rest_pos_basis_ = basis_converter.getArestBezier();
  }
  else if (par.basis == "B_SPLINE")
  {
    A_rest_pos_basis_ = basis_converter.getArestBSpline();
  }
  else
  {
    std::cout << red << "Basis " << par.basis << " not implemented yet" << reset << std::endl;
    std::cout << red << "============================================" << reset << std::endl;
    abort();
  }

  Eigen::DiagonalMatrix<double, 4> C_interval_convert;
  C_interval_convert.diagonal() << 1/std::pow(sample_T,3), 1/std::pow(sample_T,2), 1/sample_T, 1.0;

  A_rest_pos_basis_inverse_ = A_rest_pos_basis_.inverse();
  A_rest_pos_basis_t_inverse_ = (A_rest_pos_basis_*C_interval_convert).inverse();

  bspline_basis_ = basis_converter.getABSpline(par_.num_pol);
  Eigen::Matrix<double, 4, 4> M;
  M << 1, 4, 1, 0,   //////
      -3, 0, 3, 0,   //////
      3, -6, 3, 0,   //////
      -1, 3, -3, 1;  //////
  M = M / 6.0;       // *1/3!

  C_interval_convert.diagonal() << 1.0, 1/sample_T, 1/std::pow(sample_T,2), 1/std::pow(sample_T,3);

  for (int i=0; i<par_.num_pol; i++)
  {
    // bspline_basis_[i] = (bspline_basis_[i] /6 * C_interval_convert).transpose();       // *1/3!
    bspline_basis_[i] =  C_interval_convert * M;       // *1/3!

  }


  /******* kinodynamic searcher setip ******/
  num_agents_ = par_.num_of_agents;

  
  frontEndPlanner_ = std::unique_ptr<KinodynamicSearch>(new KinodynamicSearch(par_.num_pol, 
                      par_.deg_pol, par_.id, par_.safe_factor, par_.T_span, 
                      par_.num_sample_per_interval, par_.pb, par_.use_not_reaching_soln, par_.enable_entangle_check));
  frontEndPlanner_->setTetherLength(par_.tetherLength);
  frontEndPlanner_->setMaxValuesAndSamples(par_.v_max(0), par_.a_max(0), par_.j_max, 
                                    par_.a_star_samp_x);
  frontEndPlanner_->setXYZMinMaxAndRa(par_.x_min, par_.x_max, par_.y_min, par_.y_max, par_.z_min, par_.z_max, 5.0, par_.a_star_fraction_voxel_size);
  frontEndPlanner_->setBias(1.1);
  frontEndPlanner_->setGoalSize(par_.goal_radius);

  /******* end of kinodynamic searcher setip ******/

  /******* trajectory optimizer setip ******/
  backEndOptimizer_ = std::unique_ptr<PolySolverGurobi>(new PolySolverGurobi(par_.num_pol, par_.deg_pol, 
                      par_.id, par_.T_span, par_.pb, par_.weight, par_.rad_term_pos, par_.use_linear_constraints));
  backEndOptimizer_->setMaxValues(par_.x_min, par_.x_max, par_.y_min, par_.y_max, par_.z_min, par_.z_max,
                                  par_.v_max(0), par_.a_max(0), par_.j_max);
  backEndOptimizer_->setMaxRuntime(par_.runtime_opt);
  backEndOptimizer_->setTetherLength(par_.tetherLength);

  /******* end of trajectory optimizer setip ******/

  separator_solver_ = new separator::Separator();

  for (int i=0; i<num_agents_; i++)
  {
    std::vector<Eigen::Vector2d> bendPts;
    bendPtsForAgents_.push_back(bendPts);
  }
}

void Neptune::dynTraj2dynTrajCompiled(const mt::dynTraj& traj, mt::dynTrajCompiled& traj_compiled)
{
  mtx_t_.lock();
  for (auto function_i : traj.function)
  {
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    symbol_table_t symbol_table;
    symbol_table.add_variable("t", t_);
    symbol_table.add_constants();
    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    parser.compile(function_i, expression);

    traj_compiled.function.push_back(expression);
  }

  mtx_t_.unlock();

  traj_compiled.bbox = traj.bbox;
  traj_compiled.id = traj.id;
  traj_compiled.time_received = traj.time_received;  // ros::Time::now().toSec();

  traj_compiled.is_static =
      ((traj.is_agent == false) &&                           // is an obstacle and
       (traj.function[0].find("t") == std::string::npos) &&  // there is no dependence on t in the coordinate x
       (traj.function[1].find("t") == std::string::npos) &&  // there is no dependence on t in the coordinate y
       (traj.function[2].find("t") == std::string::npos))    // there is no dependence on t in the coordinate z
      ||                                                     // OR
      (traj.is_agent == true && fabs(traj.pwp.times.back() - traj.pwp.times.front()) < 1e-7);

  traj_compiled.pwp = traj.pwp;
}

void Neptune::updateTrajObstacles(mt::dynTraj traj)
{
  MyTimer tmp_t(true);

  if (started_check_ == true && traj.is_agent == true)
  {
    have_received_trajectories_while_checking_ = true;
  }

  mtx_trajs_.lock();

  std::vector<mt::dynTrajCompiled>::iterator obs_ptr =
      std::find_if(trajs_.begin(), trajs_.end(),
                   [=](const mt::dynTrajCompiled& traj_compiled) { return traj_compiled.id == traj.id; });

  bool exists_in_local_map = (obs_ptr != std::end(trajs_));

  mt::dynTrajCompiled traj_compiled;
  dynTraj2dynTrajCompiled(traj, traj_compiled);

  if (exists_in_local_map)
  {  // if that object already exists, substitute its trajectory
    *obs_ptr = traj_compiled;
  }
  else
  {  // if it doesn't exist, add it to the local map
    trajs_.push_back(traj_compiled);
    // ROS_WARN_STREAM("Adding " << traj_compiled.id);
  }

  // and now let's delete those trajectories of the obs/agents whose current positions are outside the local map
  // Note that these positions are obtained with the trajectory stored in the past in the local map
  std::vector<int> ids_to_remove;

  for (int index_traj = 0; index_traj < trajs_.size(); index_traj++)
  {
    bool traj_affects_me = false;

    mtx_t_.lock();
    t_ = ros::Time::now().toSec();

    Eigen::Vector3d center_obs;
    center_obs << trajs_[index_traj].function[0].value(),  ////////////////////
        trajs_[index_traj].function[1].value(),            ////////////////
        trajs_[index_traj].function[2].value();            /////////////////

    mtx_t_.unlock();

  }

  for (auto id : ids_to_remove)
  {
    // ROS_WARN_STREAM("Removing " << id);
    trajs_.erase(
        std::remove_if(trajs_.begin(), trajs_.end(), [&](mt::dynTrajCompiled const& traj) { return traj.id == id; }),
        trajs_.end());
  }

  mtx_trajs_.unlock();

  have_received_trajectories_while_checking_ = false;
  // std::cout << bold << blue << "updateTrajObstacles took " << tmp_t << reset << std::endl;
}



mt::ConvexHullsOfCurves_Std2d Neptune::convexHullsOfCurves2d(double t_start, double t_end, 
                                                  mt::ConvexHullsOfCurves_Std2d& result2)
{
  result2.clear();
  mt::ConvexHullsOfCurves_Std2d result;

  // for (auto traj : trajs_)
  // {
  //   result.push_back(convexHullsOfCurve2d(traj, t_start, t_end));
  // }

  for (int i=1; i<=num_agents_ ; i++) //id begin from 1 or 0? need to check
  {
    if (i==par_.id)
    {
      mt::ConvexHullsOfCurve_Std2d dummy; //push an empty vector
      result2.push_back(dummy);
    } 
    else
    {
      std::vector<mt::dynTrajCompiled>::iterator agent_ptr =
      std::find_if(trajs_.begin(), trajs_.end(),
                   [=](const mt::dynTrajCompiled& traj_compiled) { return traj_compiled.id == i; });    

      bool exists_in_traj = (agent_ptr != std::end(trajs_));  
      if (exists_in_traj)
      {
        // std::cout<<bold<<red<<"Exists in current traj!!"<<std::endl;
        mt::ConvexHullsOfCurve_Std2d tmp_hulls;
        result.push_back(convexHullsOfCurve2d(*agent_ptr, t_start, t_end, tmp_hulls));
        result2.push_back(tmp_hulls);
      } 
      else
      {
        // std::cout<<bold<<red<<"not exists in current traj!!"<<std::endl;

        mt::ConvexHullsOfCurve_Std2d dummy; //push an empty vector
        result2.push_back(dummy);      
      }
    }
  }

  return result;
}

mt::ConvexHullsOfCurve_Std2d Neptune::convexHullsOfCurve2d(mt::dynTrajCompiled& traj, double t_start, double t_end, 
                                                        mt::ConvexHullsOfCurve_Std2d& convexHulls2)
{
  mt::ConvexHullsOfCurve_Std2d convexHulls;
  double deltaT = (t_end - t_start) / (1.0 * par_.num_pol);  // num_pol is the number of intervals
  // std::cout<<bold<<red<<"adding obstacle for traj id "<<traj.id<<std::endl;

  for (int i = 0; i < par_.num_pol; i++)
  {
    mt::Polygon_Std polygon_tmp;
    convexHulls.push_back(convexHullOfInterval2d(traj, t_start + i * deltaT, 
                                  t_start + (i + 1) * deltaT, polygon_tmp));
    convexHulls2.push_back(polygon_tmp);
  }

  return convexHulls;
}

// this one includes finding convex hull without inflation
mt::Polygon_Std Neptune::convexHullOfInterval2d(mt::dynTrajCompiled& traj, double t_start, double t_end,
                                              mt::Polygon_Std& polygon2)
{
  std::vector<Eigen::Vector2d> points2;
  std::vector<Eigen::Vector2d> points = vertexesOfInterval2d(traj, t_start, t_end, points2);
  // std::cout<<bold<<green<<"vetex generated for this interval are "<<points<<std::endl;

  std::vector<Point_2> points_cgal;
  for (Eigen::Vector2d point_i : points)
  {
    points_cgal.push_back(Point_2(point_i.x(), point_i.y()));
  }

  std::vector<Point_2> points_cgal2;
  for (Eigen::Vector2d point_i : points2)
  {
    points_cgal2.push_back(Point_2(point_i.x(), point_i.y()));
  }

  polygon2 = cu::convexHullOfPoints2d(points_cgal2);
  return cu::convexHullOfPoints2d(points_cgal);
}

// this one do not include finding convex hull without inflation
mt::Polygon_Std Neptune::convexHullOfInterval2d(mt::dynTrajCompiled& traj, double t_start, double t_end)
{
  std::vector<Eigen::Vector2d> points2;
  std::vector<Eigen::Vector2d> points = vertexesOfInterval2d(traj, t_start, t_end, points2);
  // std::cout<<bold<<green<<"vetex generated for this interval are "<<points<<std::endl;

  std::vector<Point_2> points_cgal;
  for (Eigen::Vector2d point_i : points)
  {
    points_cgal.push_back(Point_2(point_i.x(), point_i.y()));
  }

  return cu::convexHullOfPoints2d(points_cgal);
}

// return a vector that contains all the vertexes of the polyhedral approx of an interval.
std::vector<Eigen::Vector2d> Neptune::vertexesOfInterval2d(mt::dynTrajCompiled& traj, double t_start, double t_end,
                                                        std::vector<Eigen::Vector2d>& points2)
{
  Eigen::Vector3d delta = Eigen::Vector3d::Zero();
  if (traj.is_agent == false) //not our concern atm
  {
    std::vector<Eigen::Vector2d> points;
    return points;
  }
  else
  {  // is an agent --> use the pwp field

    delta = traj.bbox / 2.0 + (par_.drone_radius) * Eigen::Vector3d::Ones();
    // std::cout << "****traj.bbox = " << traj.bbox << std::endl;
    // std::cout << "****par_.drone_radius = " << par_.drone_radius << std::endl;
    // std::cout << "****Inflation by delta= " << delta.transpose() << std::endl;

    return vertexesOfInterval2d(traj.pwp, t_start, t_end, delta, points2);
  }
}

std::vector<Eigen::Vector2d> Neptune::vertexesOfInterval2d(mt::PieceWisePol& pwp, double t_start, double t_end,
                                        const Eigen::Vector3d& delta, std::vector<Eigen::Vector2d>& points2)
{
  std::vector<Eigen::Vector2d> points;
  // std::cout<<bold<<blue<<"pwp.times are "<<std::endl;
  // std::cout<<bold<<blue<<pwp.times<<std::endl;

  // std::cout<<bold<<green<<"time interval we are looking is "<<std::endl;
  // std::cout<<bold<<green<<"t_start "<< t_start<<std::endl;
  // std::cout<<bold<<green<<"t_end "<< t_end<<std::endl;

  // if (t_start > pwp.times.back())
  // {
  //   Eigen::Matrix<double, 2, 4> P;
  //   P.row(0) = pwp.coeff_x.back();
  //   P.row(1) = pwp.coeff_y.back();
  //   Eigen::Matrix<double, 4, 1> T;
  //   T << par_.T_span*par_.T_span*par_.T_span, par_.T_span*par_.T_span, par_.T_span, 1;
  //   Vector2d x2d = P * T;

  //   points.push_back(Eigen::Vector2d(x2d(0) + delta.x(), x2d(1) + delta.y()));
  //   points.push_back(Eigen::Vector2d(x2d(0) + delta.x(), x2d(1) - delta.y()));
  //   points.push_back(Eigen::Vector2d(x2d(0) - delta.x(), x2d(1) - delta.y()));
  //   points.push_back(Eigen::Vector2d(x2d(0) - delta.x(), x2d(1) + delta.y()));

  //   points2.push_back(x2d);

  //   return points;
  // }

  std::vector<double>::iterator low = std::lower_bound(pwp.times.begin(), pwp.times.end(), t_start);
  std::vector<double>::iterator up = std::upper_bound(pwp.times.begin(), pwp.times.end(), t_end);

  int index_first_interval = low - pwp.times.begin() - 1;  // index of the interval [1,2]
  int index_last_interval = up - pwp.times.begin() - 1;    // index of the interval [5,6]

  // std::cout<<bold<<blue<<"index_first_interval is "<< index_first_interval<<std::endl;
  // std::cout<<bold<<blue<<"index_last_interval is "<< index_last_interval<<std::endl;

  mu::saturate(index_first_interval, 0, (int)(pwp.coeff_x.size() - 1));
  mu::saturate(index_last_interval, 0, (int)(pwp.coeff_x.size() - 1));

  Eigen::Matrix<double, 2, 4> P;
  Eigen::Matrix<double, 2, 4> V;

  // push all the complete intervals
  for (int i = index_first_interval; i <= index_last_interval; i++)
  {
    P.row(0) = pwp.coeff_x[i];
    P.row(1) = pwp.coeff_y[i];
    double _t;
    if (i!=index_last_interval)
    {
       _t = pwp.times[i+1] - pwp.times[i];
    }
    else if (t_end > pwp.times[i+1])
    {
      _t = pwp.times[i+1] - pwp.times[i];
    }
    else
    {
      _t = t_end - pwp.times[i];
    }

    if (_t>par_.T_span)
    {
      // std::cout << bold << red <<"DEBUG: _t is larger than T_span, not correct!!" << std::endl;
      // std::cout << bold << red <<"DEBUG: _t is "<<_t << std::endl;
      _t = par_.T_span;
    }
    else if (_t<0)
    {
      std::cout << bold << red <<"DEBUG: _t is smaller than 0, not correct!!" << std::endl;
      std::cout << bold << red <<"DEBUG: _t is "<<_t << std::endl;
      _t = 0;
    }

    Eigen::DiagonalMatrix<double, 4> C_interval_convert;
    C_interval_convert.diagonal() << _t*_t*_t, _t*_t, _t, 1.0;

    V = P * C_interval_convert * A_rest_pos_basis_inverse_; 

    for (int j = 0; j < V.cols(); j++)
    {
      double x = V(0, j);
      double y = V(1, j);

      if (delta.norm() < 1e-6)
      {  // no inflation
        points.push_back(Eigen::Vector2d(x, y));
      }
      else
      {
        points.push_back(Eigen::Vector2d(x + delta.x(), y + delta.y()));
        points.push_back(Eigen::Vector2d(x + delta.x(), y - delta.y()));
        points.push_back(Eigen::Vector2d(x - delta.x(), y - delta.y()));
        points.push_back(Eigen::Vector2d(x - delta.x(), y + delta.y()));
      }
      points2.push_back(V.col(j));
    }
  }

  return points;
}



bool Neptune::IsTranslating()
{
  return (drone_status_ == DroneStatus::GOAL_SEEN || drone_status_ == DroneStatus::TRAVELING);
}



mt::SampledPointsofCurves Neptune::SamplePointsOfCurves(double t_start, double t_end)
{
  mt::SampledPointsofCurves result;

  for (int i=1; i<=num_agents_ ; i++) //id begin from 1 or 0? need to check
  {
    if (i==par_.id)
    {
      mt::SampledPointsofIntervals dummy; //push an empty vector
      result.push_back(dummy);
    } 
    else
    {
      std::vector<mt::dynTrajCompiled>::iterator agent_ptr =
      std::find_if(trajs_.begin(), trajs_.end(),
                   [=](const mt::dynTrajCompiled& traj_compiled) { return traj_compiled.id == i; });    

      bool exists_in_traj = (agent_ptr != std::end(trajs_));  
      if (exists_in_traj)
      {
        // std::cout<<bold<<red<<"Exists in current traj!!"<<std::endl;

        result.push_back(SamplePointsOfIntervals(*agent_ptr, t_start, t_end));
      } 
      else
      {
        // std::cout<<bold<<red<<"not exists in current traj!!"<<std::endl;

        mt::SampledPointsofIntervals dummy; //push an empty vector
        result.push_back(dummy);      
      }
    }
  }

  return result;
}

mt::SampledPointsofIntervals Neptune::SamplePointsOfIntervals(mt::dynTrajCompiled& traj, double t_start, double t_end)
{
  mt::SampledPointsofIntervals SampledPointsInterval;
  double deltaT = (t_end - t_start) / (1.0 * par_.num_pol);  // num_pol is the number of intervals

  for (int i = 0; i < par_.num_pol; i++)
  { 
    mt::PointsofInterval ptsofInterval(2, par_.num_sample_per_interval+1);

    for (int j=0; j<=par_.num_sample_per_interval; j++)
    {
      double t_sample = t_start + deltaT * i + deltaT/par_.num_sample_per_interval * j;
      std::vector<double>::iterator low = std::upper_bound(traj.pwp.times.begin(), 
                                          traj.pwp.times.end(), t_sample);
      // std::cout << bold <<blue<<"DEBUG: t_sample is "<<t_sample << std::endl;
      // std::cout << bold <<red<<"DEBUG: index_interval is "<<low - traj.pwp.times.begin() - 1 << std::endl;

      if (low != std::end(traj.pwp.times)) //t_sample within time range
      {
        int index_interval = low - traj.pwp.times.begin() - 1; 
        mu::saturate(index_interval, 0, (int)(traj.pwp.coeff_x.size() - 1));
        double t_exceed = t_sample - traj.pwp.times[index_interval];
        // std::cout << bold <<"DEBUG: t_exceed is "<<t_exceed << std::endl;
        if (t_exceed<0)
        {
          // std::cout << bold << green <<"DEBUG: t_exceed is smaller than 0, not correct!!not correct!!" << std::endl;
          // std::cout << bold << red <<"DEBUG: t_exceed is "<<t_exceed << std::endl;
          t_exceed = 0;

        }
        else if (t_exceed>deltaT)
        {
          // std::cout << bold <<"DEBUG: t_exceed is larger than sample interval, not correct!!" << std::endl;
          // std::cout << bold << red <<"DEBUG: t_exceed is "<<t_exceed << std::endl;          
          t_exceed=deltaT;
        }
        Eigen::Matrix<double, 2, 4> P;
        Eigen::Matrix<double, 4, 1> T;
        T << t_exceed*t_exceed*t_exceed, t_exceed*t_exceed, t_exceed, 1;
        P.row(0) = traj.pwp.coeff_x[index_interval];
        P.row(1) = traj.pwp.coeff_y[index_interval];
        ptsofInterval.col(j) = P * T;
      }
      else
      {
        int index_interval = low - traj.pwp.times.begin() - 1; 
        //is this the last or second last time?? for now treat as second last        
        double t_exceed = traj.pwp.times[index_interval] - traj.pwp.times[index_interval-1];
        // std::cout << bold <<"TIME outside pwp: t_exceed is "<<t_exceed << std::endl;

        Eigen::Matrix<double, 2, 4> P;
        Eigen::Matrix<double, 4, 1> T;
        T << t_exceed*t_exceed*t_exceed, t_exceed*t_exceed, t_exceed, 1;
        P.row(0) = traj.pwp.coeff_x[index_interval-1];
        P.row(1) = traj.pwp.coeff_y[index_interval-1];
        ptsofInterval.col(j) = P * T;
      }

    }
    SampledPointsInterval.push_back(ptsofInterval);
  }

  // std::cout<<bold<<green<<"sampled points are "<<std::endl;
  // std::cout<<bold<<green<<SampledPointsInterval<<std::endl;

  return SampledPointsInterval;  
}

void Neptune::setTerminalGoal(mt::state& term_goal)
{
  mtx_G_term.lock();
  mtx_G.lock();
  mtx_state.lock();
  mtx_planner_status_.lock();

  G_term_.pos = term_goal.pos;
  Eigen::Vector3d temp = state_.pos;
  G_.pos = G_term_.pos;
  if (drone_status_ == DroneStatus::GOAL_REACHED)
  {
    changeDroneStatus(DroneStatus::YAWING);
  }
  if (drone_status_ == DroneStatus::GOAL_SEEN)
  {
    changeDroneStatus(DroneStatus::TRAVELING);
  }
  terminal_goal_initialized_ = true;

  // std::cout << bold << red << "[FA] Received Term Goal=" << G_term_.pos.transpose() << reset << std::endl;
  // std::cout << bold << red << "[FA] Received Proj Goal=" << G_.pos.transpose() << reset << std::endl;

  mtx_state.unlock();
  mtx_G.unlock();
  mtx_G_term.unlock();
  mtx_planner_status_.unlock();
}

void Neptune::setBetas(Eigen::VectorXi& beta, Eigen::VectorXi& beta2, Eigen::VectorXi& beta3,
    std::vector<Eigen::Vector2d>& previousCheckingPos,
    std::vector<Eigen::Vector2d>& previousCheckingPosAgent)
{
  mtx_beta_.lock();
  mtx_beta2_.lock();
  beta_ = beta;
  beta2_ = beta2;
  beta3_ = beta3;
  previousCheckingPos_ = previousCheckingPos;
  previousCheckingPosAgent_ = previousCheckingPosAgent;
  mtx_beta_.unlock();
  mtx_beta2_.unlock();  
}

void Neptune::setAlphasBetas(eu::ent_state& entangle_state, std::vector<Eigen::Vector2d>& previousCheckingPos,
                           std::vector<Eigen::Vector2d>& previousCheckingPosAgent)
{
  mtx_beta_.lock();

  entangle_state_ = entangle_state;
  previousCheckingPos_ = previousCheckingPos;
  previousCheckingPosAgent_ = previousCheckingPosAgent;
  
  mtx_beta_.unlock();
}

void Neptune::setBendPtsForAgent(int agent_id, std::vector<Eigen::Vector2d>& bendPts)
{
  if (agent_id > num_agents_ || agent_id < 1)
  {
    std::cout<<bold<<red<<"setBendPtsForAgent: agent ID invalid!"<<std::endl;
    return;
  }

  mtx_bendpts_.lock();

  bendPtsForAgents_[agent_id-1] = bendPts;
  
  mtx_bendpts_.unlock();
}

void Neptune::setStaticObst(std::vector<mt::Polygon_Std>& convexHullOfStaticObs)
{
  inflatedStaticObs_.clear();
  double safe_dist = 2 * par_.drone_radius+0.2;

  for (int i=0; i<convexHullOfStaticObs.size(); i++)
  {
    std::vector<Point_2> points_cgal;

    for (int j=0; j<convexHullOfStaticObs[i].cols(); j++)
    {
      points_cgal.push_back(Point_2(convexHullOfStaticObs[i](0,j) + safe_dist, 
                                    convexHullOfStaticObs[i](1,j) + safe_dist));
      points_cgal.push_back(Point_2(convexHullOfStaticObs[i](0,j) + safe_dist, 
                                    convexHullOfStaticObs[i](1,j) - safe_dist));
      points_cgal.push_back(Point_2(convexHullOfStaticObs[i](0,j) - safe_dist, 
                                    convexHullOfStaticObs[i](1,j) - safe_dist));
      points_cgal.push_back(Point_2(convexHullOfStaticObs[i](0,j) - safe_dist,
                                    convexHullOfStaticObs[i](1,j) + safe_dist));
    }
    inflatedStaticObs_.push_back(cu::convexHullOfPoints2d(points_cgal));
  }
  
  frontEndPlanner_->setStaticObstVert(inflatedStaticObs_);
  backEndOptimizer_->setStaticObstVert(inflatedStaticObs_);
}

void Neptune::setStaticObstRep(std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, 
                             std::vector<Eigen::Vector2d> staticObsLongestDist)
{
  staticObsRep_ = staticObsRep;
  frontEndPlanner_ -> setStaticObstRep(staticObsRep, staticObsLongestDist);
}

void Neptune::getG(mt::state& G)
{
  G = G_;
}

void Neptune::getState(mt::state& data)
{
  mtx_state.lock();
  data = state_;
  mtx_state.unlock();
}

void Neptune::updateState(mt::state data)
{
  mtx_state.lock();  
  state_ = data;
  mtx_state.unlock();

  if (!terminal_goal_initialized_ &&
      plan_.size()<=1)
  {
    if (plan_.size()==1) plan_.pop_front();
    mt::state tmp;
    tmp.pos = data.pos;
    tmp.yaw = data.yaw;
    plan_.push_back(tmp);
    previous_yaw_ = tmp.yaw;
  }

  state_initialized_ = true;
}


bool Neptune::initializedStateAndTermGoal()
{
  if (!state_initialized_ || !terminal_goal_initialized_)
  {
    return false;
  }
  return true;
}




// Checks that I have not received new trajectories that affect me while doing the optimization
bool Neptune::safetyCheckAfterReplan(mt::PieceWisePol& pwp_optimized, mt::SampledPointsofCurves& SampledPointsForAll,
                                   mt::state& current)
{
  started_check_ = true;

  bool result = true;
  bool need_to_rerun_entanglecheck = false;
  for (auto traj : trajs_)
  {
    if (traj.time_received > time_init_opt_ && traj.is_agent == true)
    {
      if (trajsAndPwpAreInCollision2d(traj, pwp_optimized, pwp_optimized.times.front(), pwp_optimized.times.back()))
      {
        ROS_ERROR_STREAM("Traj collides with " << traj.id);
        result = false;  // will have to redo the optimization
        break;
      }

      SampledPointsForAll[traj.id-1] = SamplePointsOfIntervals(traj, pwp_optimized.times.front(), 
                                                               pwp_optimized.times.back());

      frontEndPlanner_->updateSPocAndbendPtsForAgent(traj.id-1, SampledPointsForAll[traj.id-1], 
                                                     bendPtsForAgents_[traj.id-1]);
      need_to_rerun_entanglecheck = true;
    }
  }

  if (par_.enable_entangle_check && need_to_rerun_entanglecheck)
  {
    eu::ent_state ent_state_begin;
    PredictAlphasBetas(SampledPointsForAll, current, ent_state_begin);    
    if (frontEndPlanner_->entangleCheckGivenPwp(pwp_optimized, ent_state_begin))
    {
      result = false;
    }
  }

  // and now do another check in case I've received anything while I was checking. Note that mtx_trajs_ is locked!
  if (have_received_trajectories_while_checking_ == true)
  {
    ROS_ERROR_STREAM("Recvd traj while checking ");
    result = false;
  }
  started_check_ = false;

  return result;
}

bool Neptune::trajsAndPwpAreInCollision2d(mt::dynTrajCompiled traj, mt::PieceWisePol pwp_generated, double t_start,
                                      double t_end)
{

  double deltaT = (t_end - t_start) / pwp_generated.coeff_x.size();  // num_pol is the number of intervals
  if (abs(deltaT-par_.T_span)>0.1)
  {
    std::cout<<bold<<red<<"deltaT and Tspan not the same!"<<std::endl;
    std::cout<<bold<<red<<"deltaT: "<<deltaT<<std::endl;
    std::cout<<bold<<red<<"T_span: "<<par_.T_span<<std::endl;
    std::cout<<bold<<red<<"size of pwp: "<<pwp_generated.coeff_x.size()<<std::endl;

    return true;
  }
  for (int i = 0; i < pwp_generated.coeff_x.size(); i++)                     // for each interval
  {
    Eigen::Matrix<double, 2, 4> P;
    Eigen::Matrix<double, 2, 4> pointsA;    
    // This is my trajectory (no inflation)
    P.row(0) = pwp_generated.coeff_x[i];
    P.row(1) = pwp_generated.coeff_y[i];

    pointsA = P * A_rest_pos_basis_t_inverse_;    

    // This is the trajectory of the other agent/obstacle
    mt::Polygon_Std pointsB = convexHullOfInterval2d(traj, t_start + deltaT * i, t_start + deltaT * (i + 1));

    if (gjk::collision(pointsB, pointsA)) return true;

    // bool satisfies_LP = separator_solver_->solveModel(pointsB, pointsA);
    // if (satisfies_LP == false)
    // {
    //   return true;
    // }    

  }

  // if reached this point, they don't collide
  return false;
}

void Neptune::resetInitialization()
{
  planner_initialized_ = false;
  state_initialized_ = false;

  terminal_goal_initialized_ = false;
}

void Neptune::yaw(double diff, mt::state& next_goal)
{
  mu::saturate(diff, -par_.dc * par_.w_max, par_.dc * par_.w_max);
  double dyaw_not_filtered;

  dyaw_not_filtered = copysign(1, diff) * par_.w_max;

  dyaw_filtered_ = (1 - par_.alpha_filter_dyaw) * dyaw_not_filtered + par_.alpha_filter_dyaw * dyaw_filtered_;
  next_goal.dyaw = dyaw_filtered_;
  next_goal.yaw = previous_yaw_ + dyaw_filtered_ * par_.dc;
}

void Neptune::getDesiredYaw(mt::state& next_goal)
{
  double diff = 0.0;
  double desired_yaw = 0.0;

  switch (drone_status_)
  {
    case DroneStatus::YAWING:
      // desired_yaw = atan2(G_term_.pos[1] - next_goal.pos[1], G_term_.pos[0] - next_goal.pos[0]);
      desired_yaw = 0.0;
      diff = desired_yaw - state_.yaw;
      break;
    case DroneStatus::TRAVELING:
    case DroneStatus::GOAL_SEEN:
      // desired_yaw = atan2(M_.pos[1] - next_goal.pos[1], M_.pos[0] - next_goal.pos[0]);
      desired_yaw = 0.0;      
      diff = desired_yaw - state_.yaw;
      break;
    case DroneStatus::GOAL_REACHED:
      next_goal.dyaw = 0.0;
      next_goal.yaw = previous_yaw_;
      return;
  }

  mu::angle_wrap(diff);
  if (fabs(diff) < 0.04 && drone_status_ == DroneStatus::YAWING)
  {
    changeDroneStatus(DroneStatus::TRAVELING);
  }
  yaw(diff, next_goal);
}

bool Neptune::getNextGoal(mt::state& next_goal, bool& last_point)
{
  if (initializedStateAndTermGoal() == false)  // || (drone_status_ == DroneStatus::GOAL_REACHED && plan_.size() == 1))
  {
    // std::cout << "Not publishing new goal!!" << std::endl;
    return false;
  }

  mtx_goals.lock();
  mtx_plan_.lock();

  next_goal.setZero();
  next_goal = plan_.front();

  if (plan_.size() > 1)
  {
    plan_.pop_front();
    last_point = false;
  }
  else
  {
    last_point = true;
  }
  getDesiredYaw(next_goal);

  previous_yaw_ = next_goal.yaw;
  next_goal.yaw = 0.0;  //added mq

  mtx_goals.unlock();
  mtx_plan_.unlock();
  return true;
}

// Debugging functions
void Neptune::changeDroneStatus(int new_status)
{
  if (new_status == drone_status_)
  {
    return;
  }

  std::cout << "Changing DroneStatus from ";
  switch (drone_status_)
  {
    case DroneStatus::YAWING:
      std::cout << bold << "YAWING" << reset;
      break;
    case DroneStatus::TRAVELING:
      std::cout << bold << "TRAVELING" << reset;
      break;
    case DroneStatus::GOAL_SEEN:
      std::cout << bold << "GOAL_SEEN" << reset;
      break;
    case DroneStatus::GOAL_REACHED:
      std::cout << bold << "GOAL_REACHED" << reset;
      break;
  }
  std::cout << " to ";

  switch (new_status)
  {
    case DroneStatus::YAWING:
      std::cout << bold << "YAWING" << reset;
      break;
    case DroneStatus::TRAVELING:
      std::cout << bold << "TRAVELING" << reset;
      break;
    case DroneStatus::GOAL_SEEN:
      std::cout << bold << "GOAL_SEEN" << reset;
      break;
    case DroneStatus::GOAL_REACHED:
      std::cout << bold << "GOAL_REACHED" << reset;
      break;
  }

  std::cout << std::endl;

  drone_status_ = new_status;
}


void Neptune::PredictBetas(mt::SampledPointsofCurves& SampledPointsForAll, mt::state& current,
                         Eigen::VectorXi& beta_at_A, Eigen::VectorXi& beta2_at_A, Eigen::VectorXi& beta3_at_A)
{
  beta_at_A = beta_;
  beta2_at_A = beta2_;
  beta3_at_A = beta3_;

  Eigen::Vector2d pkplus1(current.pos(0), current.pos(1));

  for (int i=0; i<num_agents_; i++) //predict beta states
  {
    if (i==par_.id-1) continue;
    if (!SampledPointsForAll[i].empty())
    {
      Eigen::Vector2d pikplus1 = SampledPointsForAll[i][0].col(0);
      // eu::checkEntangleUpdate(beta_[i], previousCheckingPos_[i], pkplus1, 
      //                           previousCheckingPosAgent_[i], pikplus1, 
      //                           par_.pb[par_.id-1], par_.pb[i], par_.tetherLength, current.pos(2));  
      // eu::checkEntangleUpdate2(beta_at_A[i], beta2_at_A[i], previousCheckingPos_[i], pkplus1, 
      //                           previousCheckingPosAgent_[i], pikplus1, 
      //                           par_.pb[par_.id-1], par_.pb[i], par_.tetherLength, current.pos(2));  
      eu::checkEntangleUpdate3(beta_at_A[i], beta2_at_A[i], beta3_at_A[i], previousCheckingPos_[i], pkplus1, 
                                previousCheckingPosAgent_[i], pikplus1, 
                                par_.pb[par_.id-1], par_.pb[i], par_.tetherLength, current.pos(2));        
    }
  }

  for (int i=0; i<staticObsRep_.size(); i++)
  {
    eu::checkEntangleStaticUpdate(beta_at_A[num_agents_+i], beta2_at_A[num_agents_+i], beta3_at_A[num_agents_+i], 
                                  previousCheckingPos_[num_agents_+i], pkplus1, staticObsRep_[i].col(1), 
                                  par_.pb[par_.id-1], staticObsRep_[i].col(0), par_.tetherLength, current.pos(2));        
  }   
}

void Neptune::PredictAlphasBetas(mt::SampledPointsofCurves& SampledPointsForAll, mt::state& current,
                               eu::ent_state& entangle_state_A)
{
  entangle_state_A = entangle_state_;

  Eigen::Vector2d pkplus1(current.pos(0), current.pos(1));

  std::vector<Eigen::Vector2i> alphasToAdd;
  for (int i=0; i<num_agents_; i++) //predict entangle states
  {
    if (i==par_.id-1) continue;
    if (!SampledPointsForAll[i].empty())
    {
      Eigen::Vector2d pikplus1 = SampledPointsForAll[i][0].col(0);

      eu::entangleHSigToAddAgentInd(alphasToAdd, previousCheckingPos_[i], pkplus1,  
                                    previousCheckingPosAgent_[i], pikplus1, par_.pb[par_.id-1], 
                                    bendPtsForAgents_[i], i+1);      
      // previousCheckingPosAgent_[i] = pikplus1;
      // previousCheckingPos_[i] = pkplus1;      
    }
  }

  eu::entangleHSigToAddStatic(alphasToAdd, previousCheckingPos_[num_agents_], pkplus1,  
                              staticObsRep_, num_agents_);

  eu::addAlphaBetaToList(alphasToAdd, entangle_state_A, previousCheckingPos_[num_agents_], par_.pb, 
                         par_.pb[par_.id-1], staticObsRep_, num_agents_, bendPtsForAgents_);   
  eu::updateBendPts(entangle_state_A, pkplus1, par_.pb, 
                    par_.pb[par_.id-1], staticObsRep_, num_agents_);  
  // previousCheckingPos_[num_agents_] = pkplus1;

}

bool Neptune::replanKinodynamic(std::vector<mt::state>& X_safe_out, mt::PieceWisePol& pwp_out)
{
  MyTimer replanCB_t(true);
  frontEndPlanner_->clearProcess();

  if (initializedStateAndTermGoal() == false)
  {
    // std::cout << "Not Replanning" << std::endl;
    return false;
  }


  mtx_state.lock();
  mtx_G.lock();
  mtx_G_term.lock();

  mt::state state_local = state_;

  mt::state G_term = G_term_;  // Local copy of the terminal terminal goal

  mtx_G.unlock();
  mtx_G_term.unlock();
  mtx_state.unlock();

  // Check if we have reached the goal
  double dist_to_goal = (G_term.pos - state_local.pos).norm();
  if (dist_to_goal < par_.goal_radius)
  {
    changeDroneStatus(DroneStatus::GOAL_REACHED);
    exists_previous_pwp_ = false;
  }

  // Check if we have seen the goal in the last replan
  mtx_plan_.lock();
  double dist_last_plan_to_goal = (G_term.pos - plan_.back().pos).norm();
  mtx_plan_.unlock();
  if (dist_last_plan_to_goal < par_.goal_radius && drone_status_ == DroneStatus::TRAVELING)
  {
    changeDroneStatus(DroneStatus::GOAL_SEEN);
    std::cout << "Status changed to GOAL_SEEN!" << std::endl;
    exists_previous_pwp_ = false;
  }

  // Don't plan if drone is not traveling
  if (drone_status_ == DroneStatus::GOAL_REACHED) 
  //  || (drone_status_ == DroneStatus::GOAL_SEEN)) //for now, keep planning until reach goal, this helps curb infeasibility problem
    // || (drone_status_ == DroneStatus::YAWING) 
  {
    // std::cout << "No replanning needed because" << std::endl;
    // print_status();
    // return false;
  }

  std::cout << bold << on_white << "**********************IN REPLAN CB*******************" << reset << std::endl;

  //////////////////////////////////////////////////////////////////////////
  ///////////////////////// Select mt::state A /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  mt::state A;
  int k_index_end, k_index;

  // If k_index_end=0, then A = plan_.back() = plan_[plan_.size() - 1]

  mtx_plan_.lock();

  mu::saturate(deltaT_, par_.lower_bound_runtime_snlopt / par_.dc, par_.upper_bound_runtime_snlopt / par_.dc);

  int future_index = (int)(plan_.size() - deltaT_);
  k_index_end = std::max(future_index, 0);

  if (plan_.size() < ceil(par_.T_span/par_.dc))
  {
    k_index_end = 0;
  }

  k_index = plan_.size() - 1 - k_index_end;
  A = plan_.get(k_index);
  if (future_index<0)
  {
    A.vel.setZero();
    A.accel.setZero();
  }

  mtx_plan_.unlock();
  // std::cout << red << "A pos is:" << A.pos << reset << std::endl;

  // std::cout << blue << "k_index:" << k_index << reset << std::endl;
  // std::cout << blue << "k_index_end:" << k_index_end << reset << std::endl;
  // std::cout << blue << "plan_.size():" << plan_.size() << reset << std::endl;

  double runtime_snlopt;

  if (k_index_end != 0)
  {
    runtime_snlopt = k_index * par_.dc;  // std::min(, par_.upper_bound_runtime_snlopt);
  }
  else
  {
    runtime_snlopt = par_.upper_bound_runtime_snlopt;  // I'm stopped at the end of the trajectory --> take my
                                                       // time to replan
  }
  mu::saturate(runtime_snlopt, par_.lower_bound_runtime_snlopt, par_.upper_bound_runtime_snlopt);
  frontEndPlanner_->setRunTime(runtime_snlopt);

  // std::cout << green << "Runtime snlopt= " << runtime_snlopt << reset << std::endl;


  //////////////////////////////////////////////////////////////////////////
  ///////////////////////// Get point G ////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  // mt::state initial = A;
  // mt::state final = G_term;


  //////////////////////
  double time_now = ros::Time::now().toSec();  // TODO this ros dependency shouldn't be here
  double t_start = k_index * par_.dc + time_now;
  double t_final = t_start + par_.T_span * par_.num_pol;

  std::vector<Eigen::Matrix<double, 4, 1>> coeff_z_init;
  getInitialZPwp(A, G_term.pos[2], coeff_z_init);
  frontEndPlanner_->setInitZCoeffs(coeff_z_init);

  time_init_opt_ = ros::Time::now().toSec();
  mt::ConvexHullsOfCurves_Std2d tmp_hull;
  mtx_trajs_.lock();
  mt::SampledPointsofCurves SampledPointsForAll = SamplePointsOfCurves(t_start, t_final);
  mt::ConvexHullsOfCurves_Std2d hulls2d = convexHullsOfCurves2d( t_start, t_final, tmp_hull);
  mtx_trajs_.unlock();

  eu::ent_state entangle_state_A;

  mtx_beta_.lock();
  mtx_bendpts_.lock();
  // PredictBetas(SampledPointsForAll, A, beta_at_A, beta2_at_A, beta3_at_A);
  // frontEndPlanner_->setUp(A, G_term.pos, hulls2d, SampledPointsForAll, beta_at_A, beta2_at_A, beta3_at_A);
  PredictAlphasBetas(SampledPointsForAll, A, entangle_state_A);
  frontEndPlanner_->setUp(A, G_term.pos, hulls2d, SampledPointsForAll, entangle_state_A, bendPtsForAgents_);

  mtx_beta_.unlock();
  mtx_bendpts_.unlock();

  int search_result;
  std::vector<Eigen::Vector3d> currentSamplePath;
  bool got_solution = frontEndPlanner_->run(currentSamplePath, search_result);

  //////////////////////
  std::cout << on_cyan << bold << "Solved so far" << solutions_found_ << "/" << total_replannings_ << reset
            << std::endl;
  std::cout << on_cyan << bold <<  "Current dist to goal is "<<dist_to_goal << reset
            << std::endl;

  total_replannings_++;
  if (search_result==0) //runtime reached
  {
    int states_last_replan = ceil(replanCB_t.ElapsedMs() / (par_.dc * 1000));  // Number of states that
                                                                               // would have been needed for
                                                                               // the last replan
    deltaT_ = std::max(par_.factor_alpha * states_last_replan, 1.0);
    deltaT_ = std::min(1.0 * deltaT_, 2.0 / par_.dc);
    // return false;
  }

  if (!got_solution) //|| search_result!=1) //add search_result!=1 if we do not want to execute local safe path
  {
    std::cout << on_cyan << bold <<  "returning with no solution "<< reset
            << std::endl;
    return false;
  }

  // if (search_result!=1) //only local safe trajectory found, not goal-reaching
  // {
  //   // check if it is because goal is occupied by some other agent, if so, we can still use the local safe traj
  //   bool goal_occupied = false;
  //   Eigen::Matrix<double, 2, 4> goal_hull;
  //   goal_hull << G_term.pos(0) + par_.drone_radius, G_term.pos(0) + par_.drone_radius,
  //                G_term.pos(0) - par_.drone_radius, G_term.pos(0) - par_.drone_radius,
  //                G_term.pos(1) + par_.drone_radius, G_term.pos(1) - par_.drone_radius,
  //                G_term.pos(1) + par_.drone_radius, G_term.pos(1) - par_.drone_radius;

  //   for (int obs_index=0; obs_index<hulls2d.size(); obs_index++)
  //   {
  //     if (gjk::collision(hulls2d[obs_index].back(), goal_hull))
  //     {
  //       goal_occupied = true;
  //       break;
  //     }
  //   }

  //   if (!goal_occupied) return false;
  // }

  solutions_found_++;

  // get pwp out
  mt::PieceWisePol pwp_now;
  std::vector<mt::state> traj_out;
  frontEndPlanner_->generatePwpOut(pwp_now, traj_out, t_start, par_.dc);

  mtx_trajs_.lock();
  mtx_beta_.lock();
  mtx_bendpts_.lock();  
  bool is_safe_after_replan = safetyCheckAfterReplan(pwp_now, SampledPointsForAll, A);
  mtx_trajs_.unlock();
  mtx_beta_.unlock();
  mtx_bendpts_.unlock();  
  // std::cout << bold << "Check Timer=" << check_t << std::endl;

  if (is_safe_after_replan == false)
  {
    ROS_ERROR_STREAM("safetyCheckAfterOpt is not satisfied, returning");
    return false;
  }

  M_ = G_term;

  //////////////////////////////////////////////////////////////////////////
  ///////////////////////// Append to plan /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  mtx_plan_.lock();

  int plan_size = plan_.size();

  if ((plan_size - 1 - k_index_end) < 0)
  {
    std::cout << bold << red << "Already published the point A" << reset << std::endl;
    // std::cout << "plan_size= " << plan_size << std::endl;
    // std::cout << "k_index_end= " << k_index_end << std::endl;
    mtx_plan_.unlock();
    return false;
  }

  // std::cout << "Appending" << std::endl;
  std::cout << "before, plan_size=" << plan_.size() << std::endl;
  plan_.erase(plan_.end() - k_index_end - 1, plan_.end());  // this deletes also the initial condition...
  // std::cout << "middle, plan_size=" << plan_.size() << " sol.size()=" << (solver_->traj_solution_).size()
  // << std::endl;
  for (int i = 0; i < traj_out.size(); i++)  //... which is included in traj_solution_[0]
  {
    plan_.push_back(traj_out[i]);
    // std::cout << green << "traj_out pos is:" << traj_out[i].pos << reset << std::endl;

  }
  std::cout << "after, plan_size=" << plan_.size() << std::endl;
  
  mtx_plan_.unlock();

  ////////////////////
  ////////////////////

  if (exists_previous_pwp_ == true)
  {
    pwp_out = mu::composePieceWisePol(time_now, par_.dc, pwp_prev_, pwp_now);
    pwp_prev_ = pwp_out;
  }
  else
  {  //
    pwp_out = pwp_now;
    pwp_prev_ = pwp_now;
    exists_previous_pwp_ = true;
  }

  X_safe_out = plan_.toStdVector();

  ///////////////////////////////////////////////////////////
  ///////////////       OTHER STUFF    //////////////////////
  //////////////////////////////////////////////////////////

  // Check if we have planned until G_term
  mt::state F = plan_.back();  // Final point of the safe path (\equiv final point of the comitted path)
  double dist = (G_term_.pos - F.pos).norm();

  if (dist < par_.goal_radius)
  {
    changeDroneStatus(DroneStatus::GOAL_SEEN);
  }

  mtx_offsets.lock();

  int states_last_replan = ceil(replanCB_t.ElapsedMs() / (par_.dc * 1000));  // Number of states that
                                                                             // would have been needed for
                                                                             // the last replan
  deltaT_ = std::max(par_.factor_alpha * states_last_replan, 1.0);
  mtx_offsets.unlock();

  planner_initialized_ = true;

  return true;
}

bool Neptune::replanFull(std::vector<mt::state>& X_safe_out, mt::PieceWisePol& pwp_out)
{
  kino_success_ = false;
  opt_success_ = false;
  kino_compute_time_ = 100000000000000.0;
  opt_compute_time_ = 100000000000000.0;

  MyROSTimer replanCB_t(true);
  frontEndPlanner_->clearProcess();

  if (initializedStateAndTermGoal() == false)
  {
    // std::cout << "Not Replanning" << std::endl;
    return false;
  }


  mtx_state.lock();
  mtx_G.lock();
  mtx_G_term.lock();

  mt::state state_local = state_;

  mt::state G_term = G_term_;  // Local copy of the terminal terminal goal

  mtx_G.unlock();
  mtx_G_term.unlock();
  mtx_state.unlock();

  // Check if we have reached the goal
  double dist_to_goal = (G_term.pos - state_local.pos).norm();
  if (dist_to_goal < par_.goal_radius*2.0)
  {
    changeDroneStatus(DroneStatus::GOAL_REACHED);
    exists_previous_pwp_ = false;
  }

  // Check if we have seen the goal in the last replan
  mtx_plan_.lock();
  double dist_last_plan_to_goal = (G_term.pos - plan_.back().pos).norm();
  mtx_plan_.unlock();
  if (dist_last_plan_to_goal < par_.goal_radius && drone_status_ == DroneStatus::TRAVELING)
  {
    changeDroneStatus(DroneStatus::GOAL_SEEN);
    std::cout << "Status changed to GOAL_SEEN!" << std::endl;
    exists_previous_pwp_ = false;
  }

  // Don't plan if drone is not traveling
  if (drone_status_ == DroneStatus::GOAL_REACHED) 
  //  || (drone_status_ == DroneStatus::GOAL_SEEN)) //for now, keep planning until reach goal, this helps curb infeasibility problem
    // || (drone_status_ == DroneStatus::YAWING) 
  {
    // std::cout << "No replanning needed because" << std::endl;
    // print_status();
    // return false;  //for exp, stop planning when reach goal, comment out when needed
  }

  std::cout << bold << on_white << "**********************IN REPLAN CB*******************" << reset << std::endl;

  //////////////////////////////////////////////////////////////////////////
  ///////////////////////// Select mt::state A /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  mt::state A;
  mt::state next_immediate;
  int k_index_end, k_index;

  // If k_index_end=0, then A = plan_.back() = plan_[plan_.size() - 1]

  mtx_plan_.lock();

  mu::saturate(deltaT_, par_.lower_bound_runtime_snlopt / par_.dc, par_.upper_bound_runtime_snlopt / par_.dc);

  int future_index = (int)(plan_.size() - deltaT_);
  k_index_end = std::max(future_index, 0);

  if (plan_.size() < ceil(par_.T_span/par_.dc))
  {
    k_index_end = 0;
  }

  k_index = plan_.size() - 1 - k_index_end;
  A = plan_.get(k_index);
  next_immediate = plan_.front();

  mtx_plan_.unlock();

  if (future_index<0)
  {
    A.vel.setZero();
    A.accel.setZero();
  }  

  if ((next_immediate.pos - state_local.pos).norm() > 1.0) //tracking has some problem
  {
    A.pos = state_local.pos;  //for exp
  }
  // std::cout << red << "A pos is:" << A.pos << reset << std::endl;

  // std::cout << blue << "k_index:" << k_index << reset << std::endl;
  // std::cout << blue << "k_index_end:" << k_index_end << reset << std::endl;
  // std::cout << blue << "plan_.size():" << plan_.size() << reset << std::endl;

  double runtime_search;

  if (k_index_end != 0)
  {
    runtime_search = k_index * par_.dc - par_.runtime_opt;  // std::min(, par_.upper_bound_runtime_snlopt);
  }
  else
  {
    runtime_search = par_.upper_bound_runtime_snlopt;  // I'm stopped at the end of the trajectory --> take my
                                                       // time to replan
  }
  mu::saturate(runtime_search, par_.lower_bound_runtime_snlopt - par_.runtime_opt, 
                               par_.upper_bound_runtime_snlopt - par_.runtime_opt);
  frontEndPlanner_->setRunTime(runtime_search);

  ///////////////////////// Get point G ////////////////////////////////////
  double time_now = ros::Time::now().toSec();  // TODO this ros dependency shouldn't be here
  double t_start = k_index * par_.dc + time_now;
  double t_final = t_start + par_.T_span * par_.num_pol;

  std::vector<Eigen::Matrix<double, 4, 1>> coeff_z_init;
  getInitialZPwp(A, G_term.pos[2], coeff_z_init);
  frontEndPlanner_->setInitZCoeffs(coeff_z_init);

  time_init_opt_ = ros::Time::now().toSec();
  mt::ConvexHullsOfCurves_Std2d hull2d_noInflation;
  mtx_trajs_.lock();
  mt::SampledPointsofCurves SampledPointsForAll = SamplePointsOfCurves(t_start, t_final);
  mt::ConvexHullsOfCurves_Std2d hulls2d = convexHullsOfCurves2d( t_start, t_final, hull2d_noInflation);
  mtx_trajs_.unlock();

  mtx_bendpts_.lock();
  std::vector<std::vector<Eigen::Vector2d>> bendPtsForAgents = bendPtsForAgents_;  
  mtx_bendpts_.unlock();

  mtx_beta_.lock();
  eu::ent_state entangle_state_A;
  // PredictBetas(SampledPointsForAll, A, beta_at_A, beta2_at_A, beta3_at_A);
  // frontEndPlanner_->setUp(A, G_term.pos, hulls2d, SampledPointsForAll, beta_at_A, beta2_at_A, beta3_at_A);
  PredictAlphasBetas(SampledPointsForAll, A, entangle_state_A);
  frontEndPlanner_->setUp(A, G_term.pos, hulls2d, SampledPointsForAll, entangle_state_A, bendPtsForAgents);

  mtx_beta_.unlock();

  int search_result;
  std::vector<Eigen::Vector3d> currentSamplePath;
  bool got_solution = frontEndPlanner_->run(currentSamplePath, search_result);

  //////////////////////
  std::cout << on_cyan << bold << "Solved so far" << solutions_found_ << "/" << total_replannings_ << reset
            << std::endl;
  std::cout << on_cyan << bold <<  "Current dist to goal is "<<dist_to_goal << reset
            << std::endl;

  total_replannings_++;
  kino_compute_time_ = replanCB_t.ElapsedMs();

  if (search_result==0) //runtime reached
  {
    int states_last_replan = ceil(replanCB_t.ElapsedMs() / (par_.dc * 1000));  // Number of states that
                                                                               // would have been needed for
                                                                               // the last replan
    deltaT_ = std::max(par_.factor_alpha * states_last_replan, 1.0);
    deltaT_ = std::min(1.0 * deltaT_, 2.0 / par_.dc);
    // return false;
  }

  if (!got_solution )//|| search_result!=1) //add search_result!=1 if we do not want to execute local safe path
  {
    std::cout << on_cyan << bold <<  "returning with no solution "<< reset
            << std::endl;    
    return false;
  }
  if (search_result ==1) kino_success_ = true;
  if (search_result!=1 && par_.move_when_goal_occupied) //only local safe trajectory found, not goal-reaching
  {
    // check if it is because goal is occupied by some other agent, if so, we can still use the local safe traj
    bool goal_occupied = false;
    Eigen::Matrix<double, 2, 4> goal_hull;
    goal_hull << G_term.pos(0) + par_.drone_radius, G_term.pos(0) + par_.drone_radius,
                 G_term.pos(0) - par_.drone_radius, G_term.pos(0) - par_.drone_radius,
                 G_term.pos(1) + par_.drone_radius, G_term.pos(1) - par_.drone_radius,
                 G_term.pos(1) + par_.drone_radius, G_term.pos(1) - par_.drone_radius;

    for (int obs_index=0; obs_index<hulls2d.size(); obs_index++)
    {
      if (gjk::collision(hulls2d[obs_index].back(), goal_hull))
      {
        goal_occupied = true;
        break;
      }
    }

    if (!goal_occupied) return false;
  }
  solutions_found_++;

  // get pwp out
  MyROSTimer opt_time(true);
  mt::PieceWisePol pwp_now, pwp_init;
  std::vector<eu::ent_state> entStateVec;
  std::vector<mt::state> traj_out;

  frontEndPlanner_->getPwpOut_0tstart(pwp_init);
  frontEndPlanner_->getEntStateVector(entStateVec);

  if (!par_.use_direct)
  {
    backEndOptimizer_->setInitTrajectory(pwp_init);
    backEndOptimizer_->setHulls(hulls2d);
    backEndOptimizer_->setHullsNoInflation(hull2d_noInflation);
    backEndOptimizer_->setEntStateVector(entStateVec, bendPtsForAgents);

    opt_success_ = backEndOptimizer_->optimize(objective_value_);

    // if (!opt_success_ && dist_to_goal<1.0) //for exp, when it is already close to goal, 
    //                                        //only execute optimized traj
    // {
    //   return false;
    // }

    backEndOptimizer_->generatePwpOut(pwp_now, traj_out, t_start, par_.dc);
    opt_compute_time_ = opt_time.ElapsedMs();    
  }
  else //TODO: make these codes under the ddp_optimizer class
  {
    decomp_cvx_space::FlightCorridor corridor;
    int num_pol_init = pwp_init.coeff_x.size();
    for (int i=0; i<pwp_init.coeff_x.size(); i++)
    {
      corridor.appendTime(par_.T_span);
      decomp_cvx_space::Polytope empty_pltp;
      corridor.appendPolytope(empty_pltp);
    }

    double tspan4 = std::pow(par_.T_span, 4);
    double tspan3 = std::pow(par_.T_span, 3);
    double tspan2 = std::pow(par_.T_span, 2);

    Eigen::Matrix<double, 4, 4> Q_p_term;
    Eigen::Vector4d q_p_term;    
    Q_p_term << tspan3*tspan3, tspan3*tspan2, tspan4, tspan3,
               tspan3*tspan2, tspan4, tspan3, tspan2,
               tspan4, tspan3, tspan2, par_.T_span,
               tspan3, tspan2, par_.T_span, 1;   
    q_p_term << tspan3, tspan2, par_.T_span, 1; 

    MatrixXd pos = MatrixXd::Zero(2,3);
    MatrixXd vel = MatrixXd::Zero(2,3);
    MatrixXd acc = MatrixXd::Zero(2,3);
    MatrixXd jer = MatrixXd::Zero(2,3);

    pos.row(0) = A.pos.transpose();
    pos(1,0) = q_p_term.transpose() * pwp_init.coeff_x[num_pol_init-1];
    pos(1,1) = q_p_term.transpose() * pwp_init.coeff_y[num_pol_init-1];
    pos(1,2) = q_p_term.transpose() * pwp_init.coeff_z[num_pol_init-1];

    vel.row(0) = A.vel.transpose();
    acc.row(0) = A.accel.transpose();

    Eigen::DiagonalMatrix<double, 4> C_interval_convert;
    C_interval_convert.diagonal() << tspan3, tspan2, par_.T_span, 1.0;

    std::vector<Eigen::Matrix<double, 2, 4>> ctrlPtsInit_;
    for (int i=0; i<num_pol_init; i++)
    {
      Eigen::Matrix<double, 2, 4> P;
      Eigen::Matrix<double, 2, 4> Q;

      P.row(0)= pwp_init.coeff_x[i];
      P.row(1)= pwp_init.coeff_y[i];

      Q = P * C_interval_convert* A_rest_pos_basis_inverse_; //get position control points

      ctrlPtsInit_.push_back(Q);
    }
    for (int i=0; i<num_pol_init; i++)
    {
      for (int j=0; j<inflatedStaticObs_.size(); j++)
      {
        bool close_to_static_obs = false;
        double dist = (ctrlPtsInit_[i].col(0) - inflatedStaticObs_[j].col(0)).norm();
        //subtract edges of convex hull to estimate if they are close
        for (int k=0; k<3; k++)
        {
          dist -= (ctrlPtsInit_[i].col(k+1) - ctrlPtsInit_[i].col(k)).norm();
          if (dist < 0)
          {
            close_to_static_obs = true;
            break;
          }
        }
        for (int k=0; k<inflatedStaticObs_[j].cols()-1; k++)
        {
          dist -= (inflatedStaticObs_[j].col(k+1) - inflatedStaticObs_[j].col(k)).norm();
          if (dist < 0)
          {
            close_to_static_obs = true;
            break;
          }          
        }
        if (close_to_static_obs)
        {
          Eigen::Vector3d sp_line;
          bool solved = separator_solver_->solveModel(sp_line, inflatedStaticObs_[j], ctrlPtsInit_[i]);  
          if (solved)
          {
            corridor.polyhedrons[i].appendPlane(Eigen::Vector4d(sp_line(0), sp_line(1), 
                                                                0.0, sp_line(2) -1));
          }            
        }      
      }      
    }

    ddpTrajOptimizer * ddp_optimizer = new ddpTrajOptimizer();
    double w_snap = 1.0e0;
    double w_terminal = 1.0e2;
    double w_time = 1.0e2;
    int iter_max = 100;
    bool infeas = false;
    bool line_failed = true;
    int time_power = 2;
    bool minvo_flag = true;
    int rtn = ddp_optimizer->polyCurveGeneration(
                corridor, pos, vel, acc, jer, 
                par_.v_max(0), par_.a_max(0), par_.j_max*1.3, pwp_init, w_snap, w_terminal, 
                w_time, iter_max, infeas, false, false, line_failed, time_power, minvo_flag); 
    ddp_optimizer->generatePwpOut(pwp_now, traj_out, t_start, par_.dc);
    Eigen::VectorXd polytime = ddp_optimizer-> getPolyTime();
    // std::cout << bold << green << "KINODY initial time is" << 
    //                               num_pol_init*par_.T_span << std::endl;    
    std::cout << bold << green << "DIRECT optimized time is" << polytime.sum() << 
                              " vs "<<num_pol_init*par_.T_span <<std::endl;
  }

  mtx_trajs_.lock();
  mtx_beta_.lock();
  mtx_bendpts_.lock();  
  bool is_safe_after_replan = safetyCheckAfterReplan(pwp_now, SampledPointsForAll, A);
  mtx_trajs_.unlock();
  mtx_beta_.unlock();
  mtx_bendpts_.unlock();  
  // std::cout << bold << "Check Timer=" << check_t << std::endl;

  if (is_safe_after_replan == false)
  {
    ROS_ERROR_STREAM("safetyCheckAfterOpt is not satisfied, returning");
    return false;
  }

  M_ = G_term;

  //////////////////////////////////////////////////////////////////////////
  ///////////////////////// Append to plan /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  mtx_plan_.lock();

  int plan_size = plan_.size();

  if ((plan_size - 1 - k_index_end) < 0)
  {
    std::cout << bold << red << "Already published the point A" << reset << std::endl;
    // std::cout << "plan_size= " << plan_size << std::endl;
    // std::cout << "k_index_end= " << k_index_end << std::endl;
    mtx_plan_.unlock();
    return false;
  }

  // std::cout << "Appending" << std::endl;
  // std::cout << "before, plan_size=" << plan_.size() << std::endl;
  plan_.erase(plan_.end() - k_index_end - 1, plan_.end());  // this deletes also the initial condition...
  // std::cout << "middle, plan_size=" << plan_.size() << " sol.size()=" << (solver_->traj_solution_).size()
  // << std::endl;
  for (int i = 0; i < traj_out.size(); i++)  //... which is included in traj_solution_[0]
  {
    plan_.push_back(traj_out[i]);
    // std::cout << green << "traj_out pos is:" << traj_out[i].pos << reset << std::endl;

  }
  // std::cout << "after, plan_size=" << plan_.size() << std::endl;
  
  mtx_plan_.unlock();

  if (exists_previous_pwp_ == true)
  {
    pwp_out = mu::composePieceWisePol(time_now, par_.dc, pwp_prev_, pwp_now);
    pwp_prev_ = pwp_out;
  }
  else
  {  //
    pwp_out = pwp_now;
    pwp_prev_ = pwp_now;
    exists_previous_pwp_ = true;
  }

  X_safe_out = plan_.toStdVector();

  ///////////////       OTHER STUFF    //////////////////////

  // Check if we have planned until G_term
  mt::state F = plan_.back();  // Final point of the safe path (\equiv final point of the comitted path)
  double dist = (G_term_.pos - F.pos).norm();

  if (dist < par_.goal_radius)
  {
    changeDroneStatus(DroneStatus::GOAL_SEEN);
  }

  mtx_offsets.lock();

  int states_last_replan = ceil(replanCB_t.ElapsedMs() / (par_.dc * 1000));  // Number of states that
                                                                             // would have been needed for
                                                                             // the last replan
  deltaT_ = std::max(par_.factor_alpha * states_last_replan, 1.0);
  mtx_offsets.unlock();

  planner_initialized_ = true;

  return true;
}

bool Neptune::getInitialZPwp(mt::state initial_state, double z_final, std::vector<Eigen::Matrix<double, 4, 1>> &coeffs_z)
{
  ///////////////////////////
  // coeffs_z.clear();
  double p0 = initial_state.pos(2);
  double v0 = initial_state.vel(2);
  double a0 = initial_state.accel(2);

  // std::cout<<bold<<green<<"pos z initial is" <<p0<< std::endl;
  // here we saturate the value to ensure it is within the limits
  double v_max_z = par_.v_max[2];
  double a_max_z = par_.a_max[2];
  mu::saturate(v0, -v_max_z, v_max_z);
  mu::saturate(a0, -a_max_z, a_max_z);

  Eigen::VectorXd q(par_.num_pol+3);
  Eigen::VectorXd v(par_.num_pol+2);
  v(0) = 0;
  v(1) = 0;
  v(par_.num_pol) = 0;
  v(par_.num_pol+1) = 0;

  //////////////////////////////


  q(0) = p0;
  q(1) = p0 + par_.T_span * v0 / 3;
  q(2) = (3 * 3 * q(1) - 2 * par_.T_span * (-a0 * par_.T_span + v0) - 3 * (q(1) + (-2 * par_.T_span) * v0)) / 6;

  // q(par_.num_pol+2) = z_final;
  // q(par_.num_pol+1) = z_final;
  q(par_.num_pol) = z_final;

  double increment = (z_final - q(2))/(par_.num_pol-2);

  for (int i = 3; i<=par_.num_pol-1; i++)
  {
    q(i) = q(i-1) + increment;
  }


  for (int i = 3; i<= par_.num_pol; i++)
  { 
    v(i-1) = (q(i) - q(i-1))/par_.T_span;
    if (v(i-1)>v_max_z)
    {
      q(i) = q(i-1) + v_max_z * par_.T_span;
      v(i-1) = v_max_z;
    }
    else if (v(i-1)<-v_max_z)
    {
      q(i) = q(i-1) - v_max_z * par_.T_span;
      v(i-1) = -v_max_z;
    }
  }


  for (int i=2; i<= par_.num_pol-1; i++)
  {
    double a_i = (v(i) - v(i-1))/par_.T_span;
    if (a_i>a_max_z)
    {
      v(i) = v(i-1) + a_max_z * par_.T_span;
    }
    else if (a_i<-a_max_z)
    {
      v(i) = v(i-1) - a_max_z * par_.T_span;
    }
    q(i+1) = q(i) + v(i) * par_.T_span;    
  }
  q(par_.num_pol+1) = q(par_.num_pol);
  q(par_.num_pol+2) = q(par_.num_pol);  

  for (int i=0; i<par_.num_pol; i++)
  {
    Eigen::Matrix<double, 4, 1> coeff;
    coeff = (bspline_basis_[i] * q.segment(i,4)).reverse();
    coeffs_z.push_back(coeff);
  }

  // std::cout<<bold<<green<<"pos z initial poly coeff 0 is" <<coeffs_z[0](3)<< std::endl;

  return true;
}

void Neptune::getPlanningStats(bool& kino_success, double& kino_compute_time,
                             bool& opt_success, double& opt_compute_time, double& objective_value)
{
  kino_success = kino_success_;
  opt_success = opt_success_;
  kino_compute_time = kino_compute_time_;
  opt_compute_time = opt_compute_time_;
  objective_value = objective_value_;
}