/* ----------------------------------------------------------------------------
 * Copyright 2021, Cao Muqing
 * Nanyang Technological University
 * All Rights Reserved
 * Authors: Cao Muqing, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#ifndef SOLVER_GUROBI_POLY_HPP
#define SOLVER_GUROBI_POLY_HPP
#include <Eigen/Dense>
#include "gurobi_c++.h"
#include <Eigen/StdVector>

#include <iomanip>  //set precision
#include "mader_types.hpp"
#include "utils.hpp"
#include "timer.hpp"
#include <decomp_geometry/polyhedron.h>  //For Polyhedron  and Hyperplane definition
#include "separator.hpp"
#include "entangle_utils.hpp"

#include "solver_params.hpp"

typedef MADER_timers::Timer MyTimer;

class PolySolverGurobi
{
public:
  PolySolverGurobi(int num_pol, int deg_pol, int id, double T_span,  
                   std::vector<Eigen::Vector2d> pb, double weight_term, double rad_term, 
                   bool use_linear_constraints);

  ~PolySolverGurobi();

  bool optimize(double& objective_value);

  // setters
  void setMaxRuntime(double runtime);
  void setMaxValues(double x_min, double x_max, double y_min, double y_max, double z_min,
                                  double z_max, double v_max, double a_max, double j_max);
  void setInitTrajectory(mt::PieceWisePol pwp_init);
  void setHulls(mt::ConvexHullsOfCurves_Std2d &hulls);
  void setHullsNoInflation(mt::ConvexHullsOfCurves_Std2d &hulls);
  void setBetasVector(std::vector<std::vector<Eigen::Vector3d>>& vecOfAgents);
  void setTetherLength(double tetherLength);
  void setStaticObstVert(std::vector<mt::Polygon_Std>& convexHullOfStaticObs);
  void setEntStateVector(std::vector<eu::ent_state>& entStateVec, 
                          std::vector<std::vector<Eigen::Vector2d>>& bendPtsForAgents);

  void generatePwpOut(mt::PieceWisePol& pwp_out, std::vector<mt::state>& traj_out, 
  					  double t_start, double dc);
protected:
private:
  void addObjective(bool add_terminal_va_cost = false);
  void addConstraints();
  void addEntangleConstraints();
  void addEntangleConstraintForIJCase(int agent_i, int pol_num, int ent_case, 
  								      std::vector<std::vector<GRBLinExpr>>& ctrl_pt_xyz);

  GRBEnv *env_ = new GRBEnv();
  GRBModel m_ = GRBModel(*env_);

  std::vector<std::vector<GRBVar>> q_var_;      // Each q_var_[i] has 3 elements (x,y,z)
  std::vector<std::vector<GRBLinExpr>> q_exp_;  // Each q_exp_[i] has 3 elements (x,y,z)
  std::vector<std::vector<std::vector<GRBLinExpr>>> poly_exp_;  // Each q_exp_[i] has 3 elements (x,y,z)
  std::vector<std::vector<std::vector<GRBLinExpr>>> collision_plane_exp_;  // Each q_exp_[i] has 3 elements (x,y,z)

  std::vector<Eigen::Vector3d> n_;  // Each n_[i] has 3 elements (nx,ny,nz)
  std::vector<double> d_;           // d_[i] has 1 element

  // void findCentroidHull(const mt::Polyhedron_Std &hull, Eigen::Vector3d &centroid);

  // void printIndexesConstraints();
  // void printIndexesVariables();

  mt::PieceWisePol solution_;

  // int basis_ = B_SPLINE;

  int deg_pol_ = 3;
  int num_pol_ = 5;
  int num_pol_init_ = 5;

  int p_ = 5;
  int i_min_;
  int i_max_;
  // int j_min_;
  // int j_max_;
  int k_min_;
  int k_max_;
  int M_;
  int N_;

  int num_of_normals_;

  int num_of_obst_;
  int num_of_segments_;

  // nlopt::algorithm solver_;

  std::vector<Hyperplane3D> planes_;

  double dc_;
  Eigen::RowVectorXd knots_;
  double t_init_;
  double t_final_;
  double deltaT_;
  // Eigen::Vector3d v_max_;
  // Eigen::Vector3d mv_max_;
  // Eigen::Vector3d a_max_;
  // Eigen::Vector3d ma_max_;

  double weight_ = 10000;
  double weight_modified_ = 10000;

  mt::state initial_state_;
  mt::state final_state_;

  Eigen::Vector3d q0_, q1_, q2_, qNm2_, qNm1_, qN_;

  mt::ConvexHullsOfCurves_Std2d hulls_;

  MyTimer opt_timer_;

  double max_runtime_ = 0.05;  //[seconds]

  // Guesses
  std::vector<Eigen::Vector3d> n_guess_;  // Guesses for the normals
  std::vector<Eigen::Vector3d> q_guess_;  // Guesses for the normals
  std::vector<double> d_guess_;           // Guesses for the normals

  double kappa_ = 0.2;  // kappa_*max_runtime_ is spent on the initial guess
  double mu_ = 0.5;     // mu_*max_runtime_ is spent on the optimization

  // double x_min_ = -std::numeric_limits<double>::max();
  // double x_max_ = std::numeric_limits<double>::max();

  // double y_min_ = -std::numeric_limits<double>::max();
  // double y_max_ = std::numeric_limits<double>::max();

  // double z_min_ = -std::numeric_limits<double>::max();
  // double z_max_ = std::numeric_limits<double>::max();

  std::vector<double> maxs_;
  std::vector<double> mins_;

  double v_max_ = std::numeric_limits<double>::max();
  double v_min_ = -std::numeric_limits<double>::max();

  double a_max_ = std::numeric_limits<double>::max();
  double a_min_ = -std::numeric_limits<double>::max();

  double j_max_ = std::numeric_limits<double>::max();
  double j_min_ = -std::numeric_limits<double>::max();

  int num_of_QCQPs_run_ = 0;

  int a_star_samp_x_ = 7;
  int a_star_samp_y_ = 7;
  int a_star_samp_z_ = 7;

  // transformation between the B-spline control points and other basis
  std::vector<Eigen::Matrix<double, 4, 4>> M_pos_bs2basis_;
  std::vector<Eigen::Matrix<double, 3, 3>> M_vel_bs2basis_;
  std::vector<Eigen::Matrix<double, 4, 4>> A_pos_bs_;

  double a_star_bias_ = 1.0;
  double a_star_fraction_voxel_size_ = 0.5;

  separator::Separator *separator_solver_;
  // OctopusSearch *octopusSolver_;
  // KinodynamicSearch *kinoSolver_;

  std::vector<Eigen::Matrix<double, 4, 4>> M_pos_bs2basis_inverse_;  // Mbs2basis_
  Eigen::Matrix<double, 4, 4> A_rest_pos_basis_;
  Eigen::Matrix<double, 4, 4> A_rest_pos_basis_inverse_;
  Eigen::Matrix<double, 3, 3> A_rest_vel_basis_;
  Eigen::Matrix<double, 3, 3> A_rest_vel_basis_inverse321_;  
  
  Eigen::Matrix<double, 4, 4> Q_v_term_weighted;
  Eigen::Matrix<double, 4, 4> Q_a_term_weighted;
  Eigen::Matrix<double, 4, 4> Q_p_term;
  Eigen::Vector4d q_p_term;
  Eigen::Vector4d q_v_term;
  Eigen::Vector4d q_a_term;

  // int num_pol_;

  std::vector<Eigen::Vector3d> result_;


  double Ra_ = 1e10;

  // SolverCvxgen cvxgen_solver_;

  double alpha_shrink_;

  Eigen::Matrix<double, 6, 1> initial_;
  double T_span_;
  Eigen::Vector2d basepoint_;
  double cablelength_;
  int id_;
  int num_of_agents_;
  Eigen::VectorXi beta_initial_;
  Eigen::VectorXi beta2_initial_;
  Eigen::VectorXi beta3_initial_;

  std::vector<Eigen::Vector2d> pb_;
  int num_sample_per_interval_ = 1;
  std::vector<Eigen::Matrix<double, 4, 1>> sampled_time_vector_;

  bool debug_setup_complete_ = false;
  bool debug_goal_reached_ = false;
  // NodeHashTable expanded_nodes_;
  mt::PieceWisePol pwp_out_;
  double initial_z_ = 0.0;
  double goal_z_ = 0.0;
  // std::vector<mt::state> traj_out_;
  double safe_factor_ = 1.0;
  // std::vector<KNodePtr> node_pool_;
  int node_num_max_ = 20000;
  int node_used_num_ = 0;
  std::vector<Eigen::Matrix<double, 4, 1>> coeffs_z_;  
  mt::PieceWisePol pwp_init_;

  Eigen::Vector3d final_pos_;
  double rad_term_;
  mt::ConvexHullsOfCurves_Std2d hullsNoInflation_;
  std::vector<std::vector<Eigen::Vector3d>> vecOfAgentsBetas_;
  double long_length_ = 100.0;

  bool use_linear_constraints_ = true;
  std::vector<Eigen::Matrix<double, 2, 4>> ctrlPtsInit_;
  int solutions_found_ = 0;
  int total_replannings_ = 0;

  int num_of_static_obst_ = 0;
  std::vector<mt::Polygon_Std> convexHullOfStaticObs_;  
  std::vector<eu::ent_state> entStateVec_;
  std::vector<std::vector<Eigen::Vector2d>> bendPtsForAgents_;
  GRBConstr* constrTerminalVA_;
};
#endif