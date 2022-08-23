/* ----------------------------------------------------------------------------
 * Nanyang Technological University
 * Authors: Cao Muqing, et al.
 * Acknowledgement: Jesus Tordesillas
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once
#ifndef NEPTUNE_HPP
#define NEPTUNE_HPP

#include <vector>
#include "cgal_utils.hpp"
#include "entangle_utils.hpp"

#include <mutex>

#include "mader_types.hpp"
#include "kinodynamic_search.hpp"
#include "solver_gurobi_poly.hpp"


// status_ : YAWING-->TRAVELING-->GOAL_SEEN-->GOAL_REACHED-->YAWING-->TRAVELING-->...

enum DroneStatus
{
  YAWING = 0,
  TRAVELING = 1,
  GOAL_SEEN = 2,
  GOAL_REACHED = 3
};

enum PlannerStatus
{
  FIRST_PLAN = 0,
  START_REPLANNING = 1,
  REPLANNED = 2
};

using namespace termcolor;

class Neptune
{
public:
  Neptune(mt::parameters par);

  bool replanKinodynamic(std::vector<mt::state>& X_safe_out, mt::PieceWisePol& pwp_out);
  bool replanFull(std::vector<mt::state>& X_safe_out, mt::PieceWisePol& pwp_out);

  void updateState(mt::state data);

  bool getNextGoal(mt::state& next_goal,  bool& last_point);
  void getState(mt::state& data);
  void getG(mt::state& G);
  void setTerminalGoal(mt::state& term_goal);
  void resetInitialization();

  bool IsTranslating();
  void updateTrajObstacles(mt::dynTraj traj);
  void setBetas(Eigen::VectorXi& beta, Eigen::VectorXi& beta2, Eigen::VectorXi& beta3, 
    std::vector<Eigen::Vector2d>& previousCheckingPos, std::vector<Eigen::Vector2d>& previousCheckingPosAgent);

  void setAlphasBetas(eu::ent_state& entangle_state, std::vector<Eigen::Vector2d>& previousCheckingPos,
                      std::vector<Eigen::Vector2d>& previousCheckingPosAgent);
  void setStaticObst(std::vector<mt::Polygon_Std>& convexHullOfStaticObs);

  void setStaticObstRep(std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, 
                         std::vector<Eigen::Vector2d> staticObsLongestDist);
  void setBendPtsForAgent(int agent_id, std::vector<Eigen::Vector2d>& bendPts);
  void getPlanningStats(bool& kino_success, double& kino_compute_time,
                         bool& opt_success, double& opt_compute_time, double& objective_value);

/********put to public for testing purpose only***********/
  mt::ConvexHullsOfCurves_Std2d convexHullsOfCurves2d(double t_start, double t_end, 
                                                      mt::ConvexHullsOfCurves_Std2d& result2);
  mt::ConvexHullsOfCurve_Std2d convexHullsOfCurve2d(mt::dynTrajCompiled& traj, double t_start, double t_end,
                                                    mt::ConvexHullsOfCurve_Std2d& convexHulls2);
  mt::Polygon_Std convexHullOfInterval2d(mt::dynTrajCompiled& traj, double t_start, double t_end,
                                         mt::Polygon_Std& polygon2);
  mt::Polygon_Std convexHullOfInterval2d(mt::dynTrajCompiled& traj, double t_start, double t_end);

  std::vector<Eigen::Vector2d> vertexesOfInterval2d(mt::PieceWisePol& pwp, double t_start, double t_end,
                                                    const Eigen::Vector3d& delta, std::vector<Eigen::Vector2d>& points2);
  std::vector<Eigen::Vector2d> vertexesOfInterval2d(mt::dynTrajCompiled& traj, double t_start, double t_end,
                                                    std::vector<Eigen::Vector2d>& points2);

  mt::SampledPointsofCurves SamplePointsOfCurves(double t_start, double t_end);
  mt::SampledPointsofIntervals SamplePointsOfIntervals(mt::dynTrajCompiled& traj, double t_start, double t_end);
  
  bool getInitialZPwp(mt::state initial_state, double z_final, std::vector<Eigen::Matrix<double, 4, 1>> &coeffs_z);


private:
  mt::state M_;
  mt::committedTrajectory plan_;

  double previous_yaw_ = 0.0;

  void dynTraj2dynTrajCompiled(const mt::dynTraj& traj, mt::dynTrajCompiled& traj_compiled);

  bool initializedStateAndTermGoal();

  bool safetyCheckAfterReplan(mt::PieceWisePol& pwp_optimized, mt::SampledPointsofCurves& SampledPointsForAll,
                              mt::state& current);

  bool trajsAndPwpAreInCollision2d(mt::dynTrajCompiled traj, mt::PieceWisePol pwp_generated, double t_start,
                                      double t_end);

  void PredictBetas(mt::SampledPointsofCurves& SampledPointsForAll, mt::state& current, 
                    Eigen::VectorXi& beta_at_A, Eigen::VectorXi& beta2_at_A, Eigen::VectorXi& beta3_at_A);

  void PredictAlphasBetas(mt::SampledPointsofCurves& SampledPointsForAll, mt::state& current,
                          eu::ent_state& entangle_state_A);
  void yaw(double diff, mt::state& next_goal);

  void getDesiredYaw(mt::state& next_goal);

  void changeDroneStatus(int new_status);


  mt::parameters par_;

  double t_;  // variable where the expressions of the trajs of the dyn obs are evaluated

  std::mutex mtx_trajs_;
  std::vector<mt::dynTrajCompiled> trajs_;

  bool state_initialized_ = false;
  bool got_first_goal_ = false;

  bool planner_initialized_ = false;

  int deltaT_ = 75;

  bool terminal_goal_initialized_ = false;

  int drone_status_ = DroneStatus::TRAVELING;  // status_ can be TRAVELING, GOAL_SEEN, GOAL_REACHED
  int planner_status_ = PlannerStatus::FIRST_PLAN;

  double dyaw_filtered_ = 0;

  std::mutex mtx_goals;

  std::mutex mtx_k;

  std::mutex mtx_planner_status_;
  std::mutex mtx_state;
  std::mutex mtx_offsets;
  std::mutex mtx_plan_;
  // std::mutex mtx_factors;

  std::mutex mtx_G;
  std::mutex mtx_G_term;
  std::mutex mtx_t_;
  std::mutex mtx_beta_;
  std::mutex mtx_beta2_;
  std::mutex mtx_bendpts_;

  mt::state state_;
  mt::state G_;       // This goal is always inside of the map
  mt::state G_term_;  // This goal is the clicked goal

  int solutions_found_ = 0;
  int total_replannings_ = 0;

  mt::PieceWisePol pwp_prev_;

  bool exists_previous_pwp_ = false;

  bool started_check_ = false;

  bool have_received_trajectories_while_checking_ = false;

  double time_init_opt_;

  // int id_;
  // int num_sample_per_interval_ = 1;
  int num_agents_=4;
  // double av_improvement_nlopt_ = 0.0;


  std::unique_ptr<KinodynamicSearch> frontEndPlanner_; 
  std::unique_ptr<PolySolverGurobi> backEndOptimizer_; 

  Eigen::Matrix<double, 4, 4> A_rest_pos_basis_;
  Eigen::Matrix<double, 4, 4> A_rest_pos_basis_inverse_;
  Eigen::Matrix<double, 4, 4> A_rest_pos_basis_t_inverse_;

  separator::Separator* separator_solver_;
  Eigen::VectorXi beta_, beta2_, beta3_;

  std::vector<Eigen::Vector2d> previousCheckingPos_;
  std::vector<Eigen::Vector2d> previousCheckingPosAgent_;
  std::vector<Eigen::Matrix<double, 4, 4>> bspline_basis_;

  std::vector<Eigen::Matrix<double, 2, 2>> staticObsRep_;

  std::vector<Eigen::Vector2i> alphas_;
  std::vector<Eigen::Vector2d> betas_; 
  std::vector<int> active_cases_;
  eu::ent_state entangle_state_;

  std::vector<std::vector<Eigen::Vector2d>> bendPtsForAgents_;
  bool kino_success_ = false;
  bool opt_success_ = false;
  double kino_compute_time_ = 100000.0;
  double opt_compute_time_ = 100000.0;
  double objective_value_ = 0.0;
  std::vector<mt::Polygon_Std> inflatedStaticObs_;

};

#endif