/* ----------------------------------------------------------------------------
 * Nanyang Technological University
 * Authors: Cao Muqing, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once
#ifndef KINODYNAMIC_SEARCH_HPP
#define KINODYNAMIC_SEARCH_HPP

#include <vector>
#include <Eigen/Dense>
#include "mader_types.hpp"

#include "separator.hpp"
// #include "cgal_utils.hpp"
#include "entangle_utils.hpp"

#include <unordered_map>
#include <queue>
#include <tuple>
// #include <ruckig/ruckig.hpp>

//#include "solvers/cvxgen/solver_cvxgen.hpp"

typedef struct KNode KNode;  // needed to be able to have a pointer inside the struct

struct KNode
{
  KNode* previous = NULL;
  int state = 0;
  double g = 0;
  double h = 0;
  int index = 1;  // Start with q2_
  Eigen::Matrix<double, 6, 1> end; //p, v, a for two axis
  Eigen::Matrix<double, 4, 1> coeff_x; //a, b, c, d
  Eigen::Matrix<double, 4, 1> coeff_y; //a, b, c, d
  Eigen::Matrix<double, 2, 4> Q;
  Eigen::VectorXi beta;
  Eigen::VectorXi beta2;
  Eigen::VectorXi beta3; // TO DO: remove those betas

  eu::ent_state entangle_state;
};

typedef KNode* KNodePtr;

// Taken from https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
template <typename T>
struct matrix_hash_ : std::unary_function<T, size_t>
{
  std::size_t operator()(T const& matrix) const
  {
    // Note that it is obvious to the storage order of Eigen matrix (column- or
    // row-major). It will give you the same hash value for two different matrices if they
    // are the transpose of each other in different storage order.
    size_t seed = 0;
    for (size_t i = 0; i < matrix.size(); ++i)
    {
      auto elem = *(matrix.data() + i);
      seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

class NodeHashTable {
 private:
  /* data */
  std::unordered_map<Eigen::Vector3i, KNodePtr, matrix_hash_<Eigen::Vector3i>>
      data_3d_;
  std::unordered_map<Eigen::Vector2i, KNodePtr, matrix_hash_<Eigen::Vector2i>>
      data_2d_;
  std::unordered_map<Eigen::Matrix<int,6,1>, KNodePtr, matrix_hash_<Eigen::Matrix<int,6,1>>>
      data_6d_;

 public:
  NodeHashTable(/* args */) {}
  ~NodeHashTable() {}
  void insert(Eigen::Vector2i idx, KNodePtr node) {
    data_2d_.insert(std::make_pair(idx, node));
  }

  void insert(Eigen::Matrix<int,6,1> idx, KNodePtr node) {
    data_6d_.insert(std::make_pair(idx, node));
  }

  void insert(Eigen::Vector2i idx, int time_idx, KNodePtr node) {
    data_3d_.insert(std::make_pair(
        Eigen::Vector3i(idx(0), idx(1), time_idx), node));
  }

  KNodePtr find(Eigen::Vector2i idx) {
    auto iter = data_2d_.find(idx);
    return iter == data_2d_.end() ? NULL : iter->second;
  }

  KNodePtr find(Eigen::Matrix<int,6,1> idx) {
    auto iter = data_6d_.find(idx);
    return iter == data_6d_.end() ? NULL : iter->second;
  }  
  KNodePtr find(Eigen::Vector2i idx, int time_idx) {
    auto iter =
        data_3d_.find(Eigen::Vector3i(idx(0), idx(1), time_idx));
    return iter == data_3d_.end() ? NULL : iter->second;
  }

  void clear() {
    data_3d_.clear();
    data_2d_.clear();
    data_6d_.clear();
  }
};

class KinodynamicSearch
{
public:
  KinodynamicSearch(int num_pol, int deg_pol, int id,
                    double safe_factor, double T_span, int num_sample_per_interval,  
                    std::vector<Eigen::Vector2d> pb, bool use_not_reaching_soln, bool enable_entangle_check=true);  
  ~KinodynamicSearch();

  void setUp(double t_min, double t_max, const mt::ConvexHullsOfCurves_Std& hulls);
  void setUp(mt::state initial_state, Eigen::Vector3d& goal, const mt::ConvexHullsOfCurves_Std2d& hulls,
                              const mt::SampledPointsofCurves& SPoC, const Eigen::VectorXi& beta_initial);
  void setUp(mt::state initial_state, Eigen::Vector3d& goal, const mt::ConvexHullsOfCurves_Std2d& hulls,
                              const mt::SampledPointsofCurves& SPoC, const Eigen::VectorXi& beta_initial,
                              const Eigen::VectorXi& beta2_initial, const Eigen::VectorXi& beta3_initial);

  void setUp(mt::state initial_state, Eigen::Vector3d& goal, const mt::ConvexHullsOfCurves_Std2d& hulls, 
             const mt::SampledPointsofCurves& SPoC, eu::ent_state& entangle_state, 
             std::vector<std::vector<Eigen::Vector2d>>& bendPtsForAgents);

  void setMaxValuesAndSamples(double v_max, double a_max, double j_max, int num_samples);
  void setTetherLength(double tetherLength);

  void setGoal(Eigen::Vector3d& goal);

  void setRunTime(double max_runtime);
  void setGoalSize(double goal_size);

  void setBias(double bias);

  void setStaticObstRep(std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep,
                        std::vector<Eigen::Vector2d>& staticObsLongestDist);

  bool run(std::vector<Eigen::Vector3d>& result, int& status);

  bool runonce(std::vector<Eigen::Vector3d>& tmp_path, int &status);

  void recoverPath(KNodePtr node1_ptr);

  void setXYZMinMaxAndRa(double x_min, double x_max, double y_min, double y_max, 
                         double z_min, double z_max, double Ra, double voxel_size);
  void setInitZCoeffs(std::vector<Eigen::Matrix<double, 4, 1>>& coeffs_z);

  void setStaticObstVert(std::vector<mt::Polygon_Std>& convexHullOfStaticObs);

  int B_SPLINE = 1;  // B-Spline Basis
  int MINVO = 2;     // Minimum volume basis
  int BEZIER = 3;    // Bezier basis

  struct CompareCost
  {
    double bias;
    bool operator()(const KNodePtr left, const KNodePtr right)
    {
      double cost_left = left->g + bias * left->h;
      double cost_right = right->g + bias * right->h;
      if (fabs(cost_left - cost_right) < 1e-5)
      {
        return left->h > right->h;  // If two costs are ~the same, decide only upon heuristic
      }
      else
      {
        return cost_left > cost_right;
      }
    }
  };

  bool collidesWithObstaclesGivenVertexes(const Eigen::Matrix<double, 3, 4>& last4Cps, int index_lastCP);
  bool collidesWithObstacles2d(const Eigen::Matrix<double, 2, 4>& Cps, int index);
  bool collidesWithBases2d(const Eigen::Matrix<double, 2, 4>& Cps);
  bool collidesWithObstacles2dSolve(const Eigen::Matrix<double, 2, 4>& Cps, int index);

  // bool collidesWithObstacles(KNodePtr current);
  bool entanglesWithOtherAgents(KNodePtr current, double& arc_length);
  void clearProcess();

  double getCost();
  void generatePwpOut(mt::PieceWisePol& pwp_out, std::vector<mt::state>& traj_out, double t_start, double dc);
  void getPwpOut_0tstart(mt::PieceWisePol& pwp_out);
  void getBetasVector(std::vector<std::vector<Eigen::Vector3d>>& vecOfAgents);
  bool entangleCheckGivenPwp(mt::PieceWisePol& pwp, eu::ent_state& ent_state_begin);
  void updateSPocAndbendPtsForAgent(int idx, mt::SampledPointsofIntervals& SampledPtsForOne,
                                    std::vector<Eigen::Vector2d>& bendPtsForAgent);
  void getPath(std::vector<Eigen::Vector3d>& result, double& length2d);
  void getRuntime(double& runtime_this_round, double& time_spent_contact_pt, int& node_used_num);
  void getEntStateVector(std::vector<eu::ent_state>& entStateVec);

protected:
private:

  Eigen::Matrix<double, 3, 4> transformBSpline2otherBasis(const Eigen::Matrix<double, 3, 4>& Qbs, int interval);
  Eigen::Matrix<double, 3, 4> transformOtherBasis2BSpline(const Eigen::Matrix<double, 3, 4>& Qmv, int interval);


  // void plotExpandedNodesAndResult(std::vector<KNode>& expanded_nodes, KNode* result_ptr);
  // void expandAndAddToQueue(KNode& current, double constraint_xL, double constraint_xU, double constraint_yL,
                           // double constraint_yU, double constraint_zL, double constraint_zU);
  void expandAndAddToQueue(KNodePtr current);
  void expandAndAddToQueue(const Eigen::Matrix<double, 6, 1>& initstates);

  // void printPath(KNode& node1);
  // double h(KNode& node);
  double h(Eigen::Vector2d pos);
  // double g(KNode& node);
  // double weightEdge(KNode& node1, KNode& node2);

  void getCurrentSamplePath(const KNodePtr current, std::vector<Eigen::Vector3d>& output);
  void pvaToIndex(Eigen::Matrix<double, 6, 1>& input, Eigen::Matrix<int, 6, 1>& output);
  void recoverPwpOut(KNodePtr result_ptr);
  void recoverBetasVector(KNodePtr result_ptr);
  void recoverEntStateVector(KNodePtr result_ptr);

  bool entanglesWithStaticObs(KNodePtr current);
  void printEntangleState(eu::ent_state& entangle_state);
  unsigned int getIz(eu::ent_state& entangle_state);
  unsigned int power_int(unsigned int base, unsigned int exponent);
  double getH(KNodePtr node);
  double getG(KNodePtr node);

  // bias should be >=1.0
  double bias_;  // page 34 of https://www.cs.cmu.edu/~motionplanning/lecture/Asearch_v8.pdf

  int basis_ = B_SPLINE;

  Eigen::Vector3d goal_;
  Eigen::Vector2d goal2d_;

  mt::ConvexHullsOfCurves_Std2d hulls_;
  mt::SampledPointsofCurves SampledPtsForAll_;

  separator::Separator* separator_solver_;

  int num_samples_x_ = 3;
  int num_samples_y_ = 3;
  int num_samples_z_ = 3;
  int num_samples_ = 3;

  std::vector<int> indexes_samples_x_;
  std::vector<int> indexes_samples_y_;
  std::vector<int> indexes_samples_z_;

  std::vector<std::tuple<int, int>> all_combinations_;

  int p_;
  int N_;
  int M_;
  int num_of_obst_;
  int num_of_segments_;
  int num_of_normals_;

  Eigen::Vector3d orig_ = Eigen::Vector3d(0.0, 0.0, 0.0);  // origin of the search box (left bottom corner)

  Eigen::RowVectorXd knots_;

  // stores the closest node found
  double closest_dist_so_far_ = std::numeric_limits<double>::max();
  double smallest_dist_times_index_ = std::numeric_limits<double>::max();
  KNode* closest_result_so_far_ptr_ = NULL;

  // stores the closest node found that has full length (i.e. its index == (N_ - 2) )
  double complete_closest_dist_so_far_ = std::numeric_limits<double>::max();
  KNode* complete_closest_result_so_far_ptr_ = NULL;
  KNodePtr best_node_ptr_ = NULL;

  double goal_size_ = 0.5;    //[m]
  double max_runtime_ = 0.5;  //[s]

  double voxel_size_;

  double x_min_ = -std::numeric_limits<double>::max();
  double x_max_ = std::numeric_limits<double>::max();

  double y_min_ = -std::numeric_limits<double>::max();
  double y_max_ = std::numeric_limits<double>::max();

  double z_min_ = -std::numeric_limits<double>::max();
  double z_max_ = std::numeric_limits<double>::max();

  double v_max_ = std::numeric_limits<double>::max();
  double v_min_ = -std::numeric_limits<double>::max();

  double a_max_ = std::numeric_limits<double>::max();
  double a_min_ = -std::numeric_limits<double>::max();

  double j_max_ = std::numeric_limits<double>::max();
  double j_min_ = -std::numeric_limits<double>::max();

  // transformation between the B-spline control points and other basis (MINVO or Bezier)
  std::vector<Eigen::Matrix<double, 4, 4>> M_pos_bs2basis_;
  std::vector<Eigen::Matrix<double, 3, 3>> M_vel_bs2basis_;
  std::vector<Eigen::Matrix<double, 4, 4>> M_pos_bs2basis_inverse_;  // Mbs2basis_
  Eigen::Matrix<double, 4, 4> A_rest_pos_basis_;
  Eigen::Matrix<double, 4, 4> A_rest_pos_basis_inverse_;
  Eigen::Matrix<double, 3, 3> A_rest_vel_basis_;
  Eigen::Matrix<double, 3, 3> A_rest_vel_basis_inverse321_;  
  
  int num_pol_;

  std::vector<Eigen::Vector3d> result_;

  std::vector<KNodePtr> expanded_valid_nodes_;

  std::unordered_map<Eigen::Vector3i, bool, matrix_hash_<Eigen::Vector3i>> map_open_list_;

  std::priority_queue<KNodePtr, std::vector<KNodePtr>, CompareCost> openList_;  //= OpenSet, = Q

  double Ra_ = 1e10;

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
  NodeHashTable expanded_nodes_;
  mt::PieceWisePol pwp_out_;
  double initial_z_ = 0.0;
  double goal_z_ = 0.0;
  // std::vector<mt::state> traj_out_;
  double safe_factor_ = 1.0;
  std::vector<KNodePtr> node_pool_;
  int node_num_max_ = 20000;
  int node_used_num_ = 0;
  std::vector<Eigen::Matrix<double, 4, 1>> coeffs_z_;

  std::vector<std::vector<Eigen::Vector3d>> vecOfAgents_;

  int num_of_static_obst_ = 0;
  std::vector<mt::Polygon_Std> convexHullOfStaticObs_;
  std::vector<Eigen::Matrix<double, 2, 2>> staticObsRep_;

  std::vector<Eigen::Vector2i> alphas_init_;
  std::vector<Eigen::Vector2d> betas_init_; 
  std::vector<int> active_cases_init_;
  int ent_state_init_;

  eu::ent_state entangle_state_init_;
  std::vector<std::vector<Eigen::Vector2d>> bendPtsForAgents_;  
  int closest_ptr_switcher = 0, ignore_first_col_counter_ = 0;
  bool ignore_first_col_ = false;
  std::vector<Eigen::Vector2d> staticObsLongestDist_;
  // ruckig::Ruckig<2> timeOptimalGen {0.001}; 
  // ruckig::InputParameter<2> timeOptimalInput; // Number DoFs
  // ruckig::OutputParameter<2> timeOptimalOutput; // Number DoFs
  double runtime_this_round_ = 0;
  double time_spent_contact_pt_ = 0;
  bool use_not_reaching_solution_ = true;
  int ran_trigger = 0;
  bool goal_occupied_ = false;
  std::vector<eu::ent_state> entStateVec_;
  bool enable_entangle_check_ = true;  
};

#endif