/* ----------------------------------------------------------------------------
 * Nanyang Technological University
 * Authors: Cao Muqing, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */
#include <vector>
#include <random>
#include "ros/ros.h"

#include "timer.hpp"
#include "termcolor.hpp"
#include "kinodynamic_search.hpp"
#include "bspline_utils.hpp"
#include "utils.hpp"
// #include "cgal_utils.hpp"
# include "gjk.hpp"

using namespace termcolor;

typedef MADER_timers::Timer MyTimer;
typedef MADER_timers::ROSTimer MyROSTimer;

// template <typename T>
// int mu::sgn(T val)
// {
//   return (T(0) < val) - (val < T(0));
// }

KinodynamicSearch::KinodynamicSearch(int num_pol, int deg_pol, int id,double safe_factor, 
                                     double T_span, int num_sample_per_interval,  
                                     std::vector<Eigen::Vector2d> pb, bool use_not_reaching_soln,
                                     bool enable_entangle_check)
{
  p_ = deg_pol;
  M_ = num_pol + 2 * p_;
  N_ = M_ - p_ - 1;
  num_pol_ = num_pol;

  mt::basisConverter basis_converter; //convert from interval 0,1 to a user-defined one
  Eigen::Matrix<double, 4, 4> C_interval_convert;
  C_interval_convert << 
  1/std::pow(T_span,3), 0, 0, 0, 
  0, 1/std::pow(T_span,2), 0, 0, 
  0, 0, 1/T_span, 0,
  0, 0, 0, 1;

  Eigen::Matrix<double, 3, 3> Cv_interval_convert; //this one not checked yet
  Cv_interval_convert << 
  1/std::pow(T_span,2), 0, 0, 
  0, 1/T_span, 0, 
  0, 0, 1;

  std::string basis = "MINVO";
  T_span_ = T_span;
  if (basis == "MINVO")
  {
    // std::cout << green << bold << "A* is using MINVO" << reset << std::endl;
    // M_pos_bs2basis_ = basis_converter.getMinvoPosConverters(num_pol)*C_interval_convert;
    // M_vel_bs2basis_ = basis_converter.getBSplineVelConverters(num_pol)*Cv_interval_convert;  // getMinvoVelConverters TODO!!
    basis_ = MINVO;
    A_rest_pos_basis_ = basis_converter.getArestMinvo()*C_interval_convert;
    A_rest_vel_basis_ = basis_converter.getAvelrestMinvo()*Cv_interval_convert;


  }
  else if (basis == "BEZIER")
  {
    // std::cout << green << bold << "A* is using BEZIER" << reset << std::endl;
    M_pos_bs2basis_ = basis_converter.getBezierPosConverters(num_pol);
    M_vel_bs2basis_ = basis_converter.getBSplineVelConverters(num_pol);  // getBezierVelConverters TODO!!
    basis_ = BEZIER;
  }
  else if (basis == "B_SPLINE")
  {
    // std::cout << green << bold << "A* is using B_SPLINE" << reset << std::endl;

    M_pos_bs2basis_ = basis_converter.getBSplinePosConverters(num_pol);
    M_vel_bs2basis_ = basis_converter.getBSplineVelConverters(num_pol);

    basis_ = B_SPLINE;
  }
  else
  {
    std::cout << bold << red << "Basis not implemented yet" << std::endl;
    abort();
  }

  M_pos_bs2basis_inverse_.clear();
  for (auto matrix_i : M_pos_bs2basis_)
  {
    M_pos_bs2basis_inverse_.push_back(matrix_i.inverse());
  }

  A_rest_pos_basis_inverse_ = A_rest_pos_basis_.inverse();
  A_rest_vel_basis_inverse321_ = A_rest_vel_basis_.inverse();

  A_rest_vel_basis_inverse321_.row(0)=3*A_rest_vel_basis_inverse321_.row(0);
  A_rest_vel_basis_inverse321_.row(1)=2*A_rest_vel_basis_inverse321_.row(1);

  separator_solver_ = new separator::Separator();  // 0.0, 0.0, 0.0

  safe_factor_ = safe_factor;

  id_ = id;
  pb_ = pb;
  // std::cout << green << bold << "pb is" << pb<<reset << std::endl;
  // std::cout << green << bold << "id_ is" << id_<<reset << std::endl;

  basepoint_ = pb[id_-1];
  // std::cout << green << bold << "basepoint_ is" << basepoint_<<reset << std::endl;

  num_of_agents_ = pb.size();

  num_sample_per_interval_ = num_sample_per_interval;

  Eigen::Matrix<double, 4, 1> tt;
  tt.setZero();
  sampled_time_vector_.push_back(tt);

  for (int j=1; j<num_sample_per_interval; j++)
  {
    double _t = T_span*j/num_sample_per_interval;
    tt << _t*_t*_t, _t*_t, _t, 1;
    sampled_time_vector_.push_back(tt);
  }
  tt << T_span*T_span*T_span, T_span*T_span, T_span, 1;
  sampled_time_vector_.push_back(tt);
  // std::cout << green << bold << "sampled_time_vector_ is "<<sampled_time_vector_ << reset << std::endl;

  // Mbs2basis_inverse_ = Mbs2basis_.inverse();
  use_not_reaching_solution_ = use_not_reaching_soln;
  enable_entangle_check_ = enable_entangle_check;
  entangle_state_init_.active_cases.clear();
  for (int i=0; i<num_of_agents_; i++)
  {
    entangle_state_init_.active_cases.push_back((int)0);
  }
}

KinodynamicSearch::~KinodynamicSearch()
{
}

void KinodynamicSearch::setUp(mt::state initial_state, Eigen::Vector3d& goal, 
                              const mt::ConvexHullsOfCurves_Std2d& hulls,
                              const mt::SampledPointsofCurves& SPoC, const Eigen::VectorXi& beta_initial)
{
  num_of_obst_ = hulls.size();
  num_of_normals_ = num_of_segments_ * num_of_obst_;

  hulls_ = hulls;

  SampledPtsForAll_ = SPoC;

  goal2d_=goal.head(2);
  goal_z_ = goal(2);

  initial_ << initial_state.pos(0), initial_state.pos(1), 
              initial_state.vel(0), initial_state.vel(1), 
              initial_state.accel(0),initial_state.accel(1);
  beta_initial_ = beta_initial;

  initial_z_ = initial_state.pos(2);

}

void KinodynamicSearch::setUp(mt::state initial_state, Eigen::Vector3d& goal, 
                              const mt::ConvexHullsOfCurves_Std2d& hulls, const mt::SampledPointsofCurves& SPoC, 
                              const Eigen::VectorXi& beta_initial, const Eigen::VectorXi& beta2_initial, 
                              const Eigen::VectorXi& beta3_initial)
{
  num_of_obst_ = hulls.size();
  num_of_normals_ = num_of_segments_ * num_of_obst_;

  hulls_ = hulls;

  SampledPtsForAll_ = SPoC;

  goal2d_=goal.head(2);
  goal_z_ = goal(2);

  initial_ << initial_state.pos(0), initial_state.pos(1), 
              initial_state.vel(0), initial_state.vel(1), 
              initial_state.accel(0),initial_state.accel(1);
  beta_initial_ = beta_initial;
  beta2_initial_ = beta2_initial;
  beta3_initial_ = beta3_initial; 
}

void KinodynamicSearch::setUp(mt::state initial_state, Eigen::Vector3d& goal, 
                              const mt::ConvexHullsOfCurves_Std2d& hulls, 
                              const mt::SampledPointsofCurves& SPoC, eu::ent_state& entangle_state, 
                              std::vector<std::vector<Eigen::Vector2d>>& bendPtsForAgents)
{
  num_of_obst_ = hulls.size();
  num_of_normals_ = num_of_segments_ * num_of_obst_;

  hulls_ = hulls;

  SampledPtsForAll_ = SPoC;

  goal2d_=goal.head(2);
  goal_z_ = goal(2);

  initial_ << initial_state.pos(0), initial_state.pos(1), 
              initial_state.vel(0), initial_state.vel(1), 
              initial_state.accel(0),initial_state.accel(1);

  if (enable_entangle_check_) entangle_state_init_ = entangle_state;

  bendPtsForAgents_ = bendPtsForAgents;

    // check if goal is occupied by some other agent, if so, use nearest ptr to goal
  goal_occupied_ = false;
  Eigen::Matrix<double, 2, 4> goal_hull;
  double radius = 0.5;
  goal_hull << goal2d_(0) + radius, goal2d_(0) + radius,
               goal2d_(0) - radius, goal2d_(0) - radius,
               goal2d_(1) + radius, goal2d_(1) - radius,
               goal2d_(1) + radius, goal2d_(1) - radius;

  for (int obs_index=0; obs_index<hulls_.size(); obs_index++)
  {
    if (gjk::collision(hulls_[obs_index].back(), goal_hull))
    {
      goal_occupied_ = true;
      break;
    }
  }

  // ent_state_init_ = eu::NOT_ENTANGLED;
  // for (int i=0; i<num_of_agents_; i++)
  // {
  //   if (active_cases_init_[i]>1)
  //   {
  //     if (active_cases_init_[i]>2)
  //     {
  //       ent_state_init_ = eu::ENTANGLED;
  //       return;
  //     }
  //     ent_state_init_ = eu::RISK_ENTANGLING;
  //   }
  // }
  // std::cout<<bold<<blue<<"ent_state_init_ is "<<ent_state_init_<<std::endl;

  // timeOptimalInput.target_position = {goal2d_(0), goal2d_(1)};
  // timeOptimalInput.target_velocity = {0.0f, 0.0f};
  // timeOptimalInput.target_acceleration = {0.0f, 0.0f};  
}

void KinodynamicSearch::updateSPocAndbendPtsForAgent(int idx, 
      mt::SampledPointsofIntervals& SampledPtsForOne,
      std::vector<Eigen::Vector2d>& bendPtsForAgent)
{
  SampledPtsForAll_[idx] = SampledPtsForOne;
  bendPtsForAgents_[idx] = bendPtsForAgent;
}

void KinodynamicSearch::setTetherLength(double tetherLength)
{
  cablelength_ = tetherLength;
}


void KinodynamicSearch::setInitZCoeffs(std::vector<Eigen::Matrix<double, 4, 1>>& coeffs_z)
{
  coeffs_z_ = coeffs_z;
}



void KinodynamicSearch::setMaxValuesAndSamples(double v_max, double a_max, double j_max,
                                              int num_samples)
{
  all_combinations_.clear();
  indexes_samples_x_.clear();
  indexes_samples_y_.clear();
  // indexes_samples_z_.clear();

  v_max_ = v_max;
  v_min_ = -v_max;

  a_max_ = a_max;
  a_min_ = -a_max;

  j_max_ = j_max;
  j_min_ = -j_max;

  // ensure they are odd numbers (so that vx=0 is included in the samples)
  num_samples_ = (num_samples % 2 == 0) ? ceil(num_samples) : num_samples;
  // num_samples_y_ = (num_samples_y % 2 == 0) ? ceil(num_samples_y) : num_samples_y;
  // num_samples_z_ = (num_samples_z % 2 == 0) ? ceil(num_samples_z) : num_samples_z;

  for (int i = 0; i < num_samples; i++)
  {
    indexes_samples_x_.push_back(i);
  }

  for (int i = 0; i < num_samples; i++)
  {
    indexes_samples_y_.push_back(i);
  }

  // for (int i = 0; i < num_samples; i++)
  // {
  //   indexes_samples_z_.push_back(i);
  // }

  for (int jx : indexes_samples_x_)
  {
    for (int jy : indexes_samples_y_)
    {
      // for (int jz : indexes_samples_z_)
      // {
        // std::cout << "Pushing combination " << jx << ", " << jy << ", " << jz << std::endl;
        all_combinations_.push_back(std::tuple<int, int>(jx, jy));
      // }
    }
  }

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  shuffle(all_combinations_.begin(), all_combinations_.end(), std::default_random_engine(seed));

  // double min_voxel_size;
  // double max_voxel_size;

}

void KinodynamicSearch::setXYZMinMaxAndRa(double x_min, double x_max, double y_min, double y_max, double z_min,
                                      double z_max, double Ra, double voxel_size)
{
  x_min_ = x_min;
  x_max_ = x_max;

  y_min_ = y_min;
  y_max_ = y_max;

  z_min_ = z_min;
  z_max_ = z_max;
  Ra_ = Ra;

  voxel_size_ = voxel_size;

}

void KinodynamicSearch::setBias(double bias)
{
  bias_ = bias;
}

void KinodynamicSearch::setGoal(Eigen::Vector3d& goal)
{
  goal_ = goal;
}

void KinodynamicSearch::setRunTime(double max_runtime)
{
  max_runtime_ = max_runtime;
}
void KinodynamicSearch::setGoalSize(double goal_size)
{
  goal_size_ = goal_size;
}

void KinodynamicSearch::setStaticObstVert(std::vector<mt::Polygon_Std>& convexHullOfStaticObs)
{
  num_of_static_obst_ = convexHullOfStaticObs.size();
  convexHullOfStaticObs_ = convexHullOfStaticObs;

  node_num_max_ = ceil((x_max_-x_min_)*(y_max_-y_min_)/(voxel_size_*voxel_size_) 
                        * 15); //15 is just emperical
  
  node_pool_.resize(node_num_max_);
  for (int i = 0; i < node_num_max_; i++) //pre-allocate the memory for the nodes
  {
    node_pool_[i] = new KNode;
    // node_pool_[i]->entangle_state.alphas.reserve((num_of_agents_+num_of_static_obst_)*3);
    // node_pool_[i]->entangle_state.betas.reserve((num_of_agents_+num_of_static_obst_)*3);
    // node_pool_[i]->entangle_state.bendPointsIdx.reserve((num_of_agents_+num_of_static_obst_)*3);
    // node_pool_[i]->entangle_state.active_cases.reserve(num_of_agents_+num_of_static_obst_);
  }
  std::cout << green << bold << "Finished initializing kinodynamic_search " << std::endl;  
}

void KinodynamicSearch::setStaticObstRep(std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep,
                                          std::vector<Eigen::Vector2d>& staticObsLongestDist)
{
  staticObsRep_ = staticObsRep;
  staticObsLongestDist_ = staticObsLongestDist;
}

// returns the minimum voxel_size needed (if bigger, all the neighbous from q2_ are inside the voxel of q2_)
// min_voxel_size: if voxel_size < min_voxel_size, all the neighbous from q2_ will be expanded correctly
// max_voxel_size: if voxel_size > max_voxel_size, there will be no nodes expanded from q2_ (they are on the same cell
// as q2_)





double KinodynamicSearch::getH(KNodePtr node)
{
  // MyTimer timer_expand(true);
  // timeOptimalInput.current_position = {node->end(0), node->end(1)};
  // timeOptimalInput.current_velocity = {node->end(2), node->end(3)};
  // timeOptimalInput.current_acceleration = {node->end(4), node->end(5)}; 

  // if (timeOptimalGen.update(timeOptimalInput, timeOptimalOutput) == ruckig::Result::Working)
  // {
  //   double dist = (node->end.head<2>()-goal2d_).norm();
  //   if (dist<7.0)
  //   {
  //     return std::min(dist/(v_max_*sqrt(2.0)), timeOptimalOutput.trajectory.get_duration());
  //   }
  //   // std::cout<<bold<<red<<"time taken for ruckig is "<<timer_expand.ElapsedUs()<<"Us \n";
  //   return timeOptimalOutput.trajectory.get_duration();  
  // }
  // else
  // {
  //   // std::cout<<bold<<red<<"RUcikg failed! \n";
  //   return 100000.0; //a large number
  // }
  return (node->end.head<2>() - goal2d_).norm();
}

double KinodynamicSearch::getG(KNodePtr node)
{
  // return node->previous->g + T_span_;
  return node->previous->g + (node->previous->end.head<2>() - node->end.head<2>()).norm();

}

double KinodynamicSearch::h(Eigen::Vector2d pos)
{
  return (pos - goal2d_).norm();
}

// double KinodynamicSearch::g(KNode& node)
// {
//   double cost = 0;
//   KNodePtr tmp = &node;
//   while (tmp->previous != NULL)
//   {
//     cost = cost + weightEdge(*tmp->previous, *tmp);
//     tmp = tmp->previous;
//   }

//   return cost;
// }

// double KinodynamicSearch::weightEdge(KNode& node1, KNode& node2)  // edge cost when adding node 2
// {
//   return (node2.qi - node1.qi).norm();
// }

// void KinodynamicSearch::printPath(KNode& node1)
// {
//   KNode tmp = node1;
//   while (tmp.previous != NULL)
//   {
//     // std::cout << tmp.index << ", ";  // qi.transpose().x()
//     std::cout << tmp.previous->qi.transpose() << std::endl;
//     tmp = *tmp.previous;
//   }
//   std::cout << std::endl;
// }

void KinodynamicSearch::recoverPath(KNodePtr result_ptr)
{
  // std::cout << "Recovering path" << std::endl;
  result_.clear();

  if (result_ptr == NULL)
  {
    return;
  }

  KNodePtr tmp = result_ptr;

  while (tmp != NULL)
  {
    std::vector<Eigen::Vector3d> tmpPath;
    getCurrentSamplePath(tmp, tmpPath);
    for (auto i : tmpPath) result_.push_back(i);

    tmp = tmp->previous;
  }

  std::reverse(std::begin(result_), std::end(result_)); 
}

void KinodynamicSearch::getPath(std::vector<Eigen::Vector3d>& result, double& length2d)
{
  result.clear();
  length2d = 0.0;

  if (best_node_ptr_ == NULL)
  {
    return;
  }

  KNodePtr tmp = best_node_ptr_;

  while (tmp != NULL)
  {
    std::vector<Eigen::Vector3d> tmpPath;
    getCurrentSamplePath(tmp, tmpPath);
    for (auto i : tmpPath) result.push_back(i);

    tmp = tmp->previous;
  }

  std::reverse(std::begin(result), std::end(result));  
  for (int i=0; i<result.size()-1; i++)
  {
    length2d += (result[i+1].head(2) - result[i].head(2)).norm();
  }

}

void KinodynamicSearch::recoverPwpOut(KNodePtr result_ptr)
{
  // std::cout << "Recovering path" << std::endl;
  pwp_out_.clear();

  if (result_ptr == NULL)
  {
    return;
  }

  KNodePtr tmp = result_ptr;
  // Eigen::Matrix<double, 4, 1> _coeff_z;
  // _coeff_z << 0.0, 0.0, 0.0, goal_z_; //z is constant

  while (tmp != NULL)
  {
    if (tmp->index <= num_pol_)
    {
      pwp_out_.times.push_back(tmp->index*T_span_);
      pwp_out_.coeff_x.push_back(tmp->coeff_x);  
      pwp_out_.coeff_y.push_back(tmp->coeff_y);  
      pwp_out_.coeff_z.push_back(coeffs_z_[tmp->index-1]);    
    }

    tmp = tmp->previous;
  }
  pwp_out_.times.push_back(0.0);
  std::reverse(std::begin(pwp_out_.times), std::end(pwp_out_.times));  // result_ is [q0 q1 q2 q3 ...]
  std::reverse(std::begin(pwp_out_.coeff_x), std::end(pwp_out_.coeff_x));  // result_ is [q0 q1 q2 q3 ...]
  std::reverse(std::begin(pwp_out_.coeff_y), std::end(pwp_out_.coeff_y));  // result_ is [q0 q1 q2 q3 ...]
  std::reverse(std::begin(pwp_out_.coeff_z), std::end(pwp_out_.coeff_z));  // result_ is [q0 q1 q2 q3 ...]

}

void KinodynamicSearch::recoverBetasVector(KNodePtr result_ptr)
{
  // std::cout << "Recovering path" << std::endl;
  vecOfAgents_.clear();

  if (result_ptr == NULL)
  {
    return;
  }

  for (int i=0; i<num_of_agents_; i++)
  {
    KNodePtr tmp = result_ptr;
    std::vector<Eigen::Vector3d> vecOfEachAgent;
    while (tmp != NULL)
    {
      if (tmp->index <= num_pol_)
      {
        vecOfEachAgent.push_back(Vector3d(tmp->beta(i), tmp->beta2(i), tmp->beta3(i)));
      }
      tmp = tmp->previous;
    }
    std::reverse(std::begin(vecOfEachAgent), std::end(vecOfEachAgent));  // result_ is [q0 q1 q2 q3 ...]
    vecOfAgents_.push_back(vecOfEachAgent);
  }
}

void KinodynamicSearch::recoverEntStateVector(KNodePtr result_ptr)
{
  // std::cout << "Recovering path" << std::endl;
  entStateVec_.clear();

  if (result_ptr == NULL)
  {
    return;
  }

  KNodePtr tmp = result_ptr;
  while (tmp != NULL)
  {
    if (tmp->index <= num_pol_)
    {
      entStateVec_.push_back(tmp->entangle_state);
    }
    tmp = tmp->previous;
  }
  entStateVec_.push_back(entangle_state_init_);
  std::reverse(std::begin(entStateVec_), std::end(entStateVec_));  // result_ is [q0 q1 q2 q3 ...]
}

void KinodynamicSearch::getBetasVector(std::vector<std::vector<Eigen::Vector3d>>& vecOfAgents)
{
  vecOfAgents = vecOfAgents_;
}

void KinodynamicSearch::getEntStateVector(std::vector<eu::ent_state>& entStateVec)
{
  entStateVec = entStateVec_;
}

void KinodynamicSearch::getPwpOut_0tstart(mt::PieceWisePol& pwp_out)
{
  pwp_out = pwp_out_;
}


void KinodynamicSearch::generatePwpOut(mt::PieceWisePol& pwp_out, 
                                      std::vector<mt::state>& traj_out, double t_start, double dc)
{
  pwp_out = pwp_out_;
  int pwp_t_size = pwp_out.times.size();
  std::vector<Eigen::Matrix<double, 3, 4>> poly;

  for (int i=0; i<pwp_t_size; i++) //front end searcher is not aware of the time
  {
    pwp_out.times[i] += t_start;
    Eigen::Matrix<double, 3, 4> polyp;

    if (i == pwp_t_size-1) continue;
    polyp.row(0) = pwp_out.coeff_x[i];
    polyp.row(1) = pwp_out.coeff_y[i];
    polyp.row(2) = pwp_out.coeff_z[i];
    poly.push_back(polyp);

  }

  traj_out.clear();

  double _t = 0;
  int i = 0;
  double delta_t;
  while (i<pwp_t_size-1)
  { 
    delta_t = _t - i*T_span_;
    if (delta_t<0 ||delta_t>T_span_)
    {
      std::cout<<bold<<yellow<<"delta_t is not correct, something is wrong here!"<<std::endl;
    }
    Eigen::Vector4d tp(delta_t*delta_t*delta_t, delta_t*delta_t, delta_t, 1);
    Eigen::Vector3d tv(3*delta_t*delta_t, 2*delta_t, 1);
    Eigen::Vector2d ta(6*delta_t, 2);

    mt::state state_i;
    state_i.setPos(poly[i] * tp);  // First column
    state_i.setVel(poly[i].block<3, 3>(0, 0) * tv);
    state_i.setAccel(poly[i].block<3, 2>(0, 0) * ta);
    state_i.setJerk(poly[i].col(0)*6);
    traj_out.push_back(state_i);

    _t += dc;
    if (_t>(i+1)*T_span_) i++;
  }

}

Eigen::Matrix<double, 3, 4> KinodynamicSearch::transformBSpline2otherBasis(const Eigen::Matrix<double, 3, 4>& Qbs,
                                                                       int interval)
{
  return Qbs * M_pos_bs2basis_[interval];
}

Eigen::Matrix<double, 3, 4> KinodynamicSearch::transformOtherBasis2BSpline(const Eigen::Matrix<double, 3, 4>& Qmv,
                                                                       int interval)
{
  return Qmv * M_pos_bs2basis_inverse_[interval];
}


void KinodynamicSearch::getCurrentSamplePath(const KNodePtr current, std::vector<Eigen::Vector3d>& output)
{
  Eigen::Matrix<double, 3, 4> P;
  P.setZero();
  P.row(0)= current->coeff_x;
  P.row(1)= current->coeff_y;
  int num_sample_per_int = 10;
  // std::cout << red << "P is "<<P<<reset<<std::endl;

  std::vector<Eigen::Vector3d> output_tmp;

  Eigen::Matrix<double, 4, 1> tt;

  for (int j=num_sample_per_int; j>=0; j--) //start from the end point instead
  {
    double _t = (double)j/num_sample_per_int*T_span_;
    tt << _t*_t*_t, _t*_t, _t, 1;
    Eigen::Vector3d tmp = P * tt;
    // std::cout << red << "pushing back point "<<tmp<<reset<<std::endl;
    output_tmp.push_back(tmp);
  }
  output = output_tmp;
}

bool KinodynamicSearch::entanglesWithOtherAgents(KNodePtr current, double& arc_length)
{
  // Eigen::VectorXi beta_old = current->beta;
  // for (int i=1; i<=num_of_agents_; i++) //first id is 1
  // {
  //   // std::cout<<bold<<green<<"going to check entangle for ID "<<i<<std::endl;
  //   if (i==id_) continue;
  //   if (SampledPtsForAll_[i-1].empty()) continue;
  //   // std::cout<<bold<<blue<<"ID "<<i<<" has non-zero sampled points"<<std::endl;
  //   // std::cout<<bold<<green<<"current->ent_state is "<<current->ent_state<<std::endl;
  //   int active_cases_i = current->active_cases[i-1];

  //    Eigen::Matrix<double, 2, 4> P;
  //    P.row(0)= current->coeff_x;
  //    P.row(1)= current->coeff_y;

  //   if (current->index>num_pol_)
  //   {
  //     Eigen::Vector2d lastPt=SampledPtsForAll_[i-1][num_pol_-1].rightCols(1); //find the last position
  //     Eigen::Vector2d pk(current->coeff_x(3), current->coeff_y(3));
  //     Eigen::Vector2d pkplus1;
  //     // std::cout<<bold<<blue<<"lastpt is "<<lastPt<<std::endl;

  //     for (int j=1; j<num_sample_per_interval_; j++)
  //     {
  //       pkplus1 = P * sampled_time_vector_[j];

  //       eu::entangleHSigUpdate3(i, current->alphas, current->betas, current->active_cases[i-1], pk, pkplus1, 
  //                               lastPt, lastPt, basepoint_, pb_[i-1]);
  //       // if (eu::checkEntangleUpdate3(current->beta(i-1), current->beta2(i-1), current->beta3(i-1), pk, pkplus1, 
  //       //                              lastPt, lastPt, basepoint_, pb_[i-1], cablelength_, coeffs_z_[num_pol_-1](3)))
  //       if (active_cases_i>=2 && current->active_cases[i-1]>active_cases_i)                                                  
  //       {
  //         std::cout<<bold<<red<<"for agent "<<i<<"_state is "<<current->active_cases[i-1]<<std::endl;
  //         return true;
  //       }        
  //       pk = pkplus1;   
  //     }

  //     pkplus1 << current->end(0), current->end(1); //check last point

  //     eu::entangleHSigUpdate3(i, current->alphas, current->betas, current->active_cases[i-1], pk, pkplus1, 
  //                             lastPt, lastPt, basepoint_, pb_[i-1]);

  //     // if (eu::checkEntangleUpdate3(current->beta(i-1), current->beta2(i-1), current->beta3(i-1), pk, pkplus1, 
  //     //                              lastPt, lastPt, basepoint_, pb_[i-1], cablelength_, coeffs_z_[num_pol_-1](3)))  
  //      if (active_cases_i>=2 && current->active_cases[i-1]>active_cases_i)                                                  
  //     {
  //       std::cout<<bold<<red<<"for agent "<<i<<"_state is "<<current->active_cases[i-1]<<std::endl;
  //       return true;
  //     }
  //   }
  //   else
  //   {
  //     Eigen::Vector2d pik=SampledPtsForAll_[i-1][current->index-1].leftCols(1); //find the first position
  //     Eigen::Vector2d pikplus1;
  //     Eigen::Vector2d pk(current->coeff_x(3), current->coeff_y(3)); //initial point
  //     Eigen::Vector2d pkplus1;

  //     for (int j=1; j<num_sample_per_interval_; j++)
  //     {
  //       pikplus1 = SampledPtsForAll_[i-1][current->index-1].col(j);
  //       pkplus1 = P * sampled_time_vector_[j];
        
  //       eu::entangleHSigUpdate3(i, current->alphas, current->betas, current->active_cases[i-1], pk, pkplus1, 
  //                               pik, pikplus1, basepoint_, pb_[i-1]);
  //       // if (eu::checkEntangleUpdate3(current->beta(i-1), current->beta2(i-1), current->beta3(i-1), pk, pkplus1, 
  //       //                              pik, pikplus1, basepoint_, pb_[i-1], cablelength_, coeffs_z_[current->index-1](3))) 
  //       if (active_cases_i>=2 && current->active_cases[i-1]>active_cases_i)                                                  
  //       {
  //         std::cout<<bold<<red<<"for agent "<<i<<"_state is "<<current->active_cases[i-1]<<std::endl;
  //         return true;
  //       }     
  //       pk = pkplus1;   
  //       pik = pikplus1;
  //     }

  //     pkplus1 << current->end(0), current->end(1); //check last point
  //     pikplus1 = SampledPtsForAll_[i-1][current->index-1].col(num_sample_per_interval_);
  
  //     eu::entangleHSigUpdate3(i, current->alphas, current->betas, current->active_cases[i-1], pk, pkplus1, 
  //                             pik, pikplus1, basepoint_, pb_[i-1]);      
  //     // if (eu::checkEntangleUpdate3(current->beta(i-1), current->beta2(i-1), current->beta3(i-1), pk, pkplus1, 
  //     //                              pik, pikplus1, basepoint_, pb_[i-1], cablelength_, coeffs_z_[current->index-1](3)))
  //     if (active_cases_i>=2 && current->active_cases[i-1]>active_cases_i)                                                  
  //     {
  //       std::cout<<bold<<red<<"for agent "<<i<<"_state is "<<current->active_cases[i-1]<<std::endl;
  //       return true;
  //     }

  //   }    

  //   if (current->active_cases[i-1]<active_cases_i) //check further reduction 
  //   {
  //     eu::shrinkAlphaBetas(current->alphas, current->betas, current->active_cases, num_of_agents_);
  //   }
  // }

  Eigen::Matrix<double, 2, 4> P;
  P.row(0)= current->coeff_x;
  P.row(1)= current->coeff_y;  
  Eigen::Vector2d pk(current->coeff_x(3), current->coeff_y(3));
  Eigen::Vector2d pkplus1;   
  Eigen::Vector2d pik;
  Eigen::Vector2d pikplus1;

  std::vector<int> active_cases_old = current->entangle_state.active_cases;

  for (int j=1; j<=num_sample_per_interval_; j++)
  {   
    std::vector<Eigen::Vector2i> alphasToAdd;
    if (j<num_sample_per_interval_) pkplus1 = P * sampled_time_vector_[j];
    else pkplus1 << current->end(0), current->end(1); //check last point
    arc_length += (pkplus1-pk).norm();
    // std::cout<<bold<<red<<"entangleHSigToAddAgentInd for agent "<<std::endl;

    for (int i=0; i<num_of_agents_; i++) //first id is 1
    {
      if (i==id_-1) continue;
      if (SampledPtsForAll_[i].empty()) continue; 
      if (current->index>num_pol_)
      {
        pik = SampledPtsForAll_[i][num_pol_-1].rightCols(1); //find the last position
        pikplus1 = pik;
      }
      else
      {
        pik = SampledPtsForAll_[i][current->index-1].col(j-1); //find the first position
        pikplus1 = SampledPtsForAll_[i][current->index-1].col(j);
      }

      eu::entangleHSigToAddAgentInd(alphasToAdd, pk, pkplus1, pik, pikplus1, basepoint_, 
                                    bendPtsForAgents_[i], i+1);   
    }

    eu::entangleHSigToAddStatic(alphasToAdd, pk, pkplus1, staticObsRep_, num_of_agents_);

    if (current->entangle_state.alphas.size() + alphasToAdd.size() 
        >(num_of_agents_+num_of_static_obst_))
    {
      return true;
    }

    eu::addAlphaBetaToList(alphasToAdd, current->entangle_state, pk, pb_, basepoint_, 
                           staticObsRep_, num_of_agents_, bendPtsForAgents_);     

    for (int i=0; i<num_of_agents_; i++)  // check entanglement
    { 
      if (active_cases_old[i] < 2 && current->entangle_state.active_cases[i] >= 2)
      {
        // std::cout<<bold<<red<<"for agent "<<i<<"active_cases_old is "
        //          <<active_cases_old[i]<<std::endl;
        // std::cout<<bold<<red<<"for agent "<<i<<"entangle_state.active_cases is "
        //          <<current->entangle_state.active_cases[i]<<std::endl;                 
        return true;

      }
      else if (active_cases_old[i] >=2 && 
               current->entangle_state.active_cases[i] > active_cases_old[i])
      {
        // std::cout<<bold<<red<<"for agent "<<i<<"active_cases_old is "
        //          <<active_cases_old[i]<<std::endl;
        // std::cout<<bold<<red<<"for agent "<<i<<"entangle_state.active_cases is "
        //          <<current->entangle_state.active_cases[i]<<std::endl;           
        return true;
      }
    }
    
    MyROSTimer timer_function(true);
    eu::updateBendPts(current->entangle_state, pkplus1, pb_, 
                      basepoint_, staticObsRep_, num_of_agents_);
    time_spent_contact_pt_ += (double)timer_function.ElapsedMs();

    active_cases_old = current->entangle_state.active_cases;
    pk = pkplus1;   
  }

  MyROSTimer timer_function(true);
  if (eu::getTetherLength(current->entangle_state, pb_, basepoint_, pkplus1, 
      staticObsRep_, staticObsLongestDist_, num_of_agents_) > cablelength_)
  {
    // std::cout<<bold<<red<<"failing cable length check!"<<std::endl;     
    time_spent_contact_pt_ += (double)timer_function.ElapsedMs();     
    return true;
  }
  time_spent_contact_pt_ += (double)timer_function.ElapsedMs();

  return false;
}

bool KinodynamicSearch::entangleCheckGivenPwp(mt::PieceWisePol& pwp, eu::ent_state& ent_state_begin)
{
  for (int pwp_idx = 0; pwp_idx<pwp.coeff_x.size(); pwp_idx++)
  {
    Eigen::Matrix<double, 2, 4> P;
    P.row(0)= pwp.coeff_x[pwp_idx];
    P.row(1)= pwp.coeff_y[pwp_idx];  
    Eigen::Vector2d pk(pwp.coeff_x[pwp_idx](3), pwp.coeff_y[pwp_idx](3));
    Eigen::Vector2d pkplus1;   
    Eigen::Vector2d pik;
    Eigen::Vector2d pikplus1;

    std::vector<int> active_cases_old = ent_state_begin.active_cases;

    for (int j=1; j<=num_sample_per_interval_; j++)
    {   
      std::vector<Eigen::Vector2i> alphasToAdd;
      pkplus1 = P * sampled_time_vector_[j];
      // std::cout<<bold<<red<<"entangleHSigToAddAgentInd for agent "<<std::endl;

      for (int i=0; i<num_of_agents_; i++) //first id is 1
      {
        if (i==id_-1) continue;
        if (SampledPtsForAll_[i].empty()) continue; 
        if (pwp_idx>num_pol_-1)
        {
          pik = SampledPtsForAll_[i][num_pol_-1].rightCols(1); //find the last position
          pikplus1 = pik;
        }
        else
        {
          pik = SampledPtsForAll_[i][pwp_idx].col(j-1); //find the first position
          pikplus1 = SampledPtsForAll_[i][pwp_idx].col(j);
        }

        eu::entangleHSigToAddAgentInd(alphasToAdd, pk, pkplus1, pik, pikplus1, basepoint_, 
                                      bendPtsForAgents_[i], i+1);   
      }

      eu::entangleHSigToAddStatic(alphasToAdd, pk, pkplus1, staticObsRep_, num_of_agents_);

      if (ent_state_begin.alphas.size() + alphasToAdd.size() 
          >(num_of_agents_+num_of_static_obst_)*3)
      {
        return true;
      }

      eu::addAlphaBetaToList(alphasToAdd, ent_state_begin, pk, pb_, basepoint_, 
                             staticObsRep_, num_of_agents_, bendPtsForAgents_);     

      for (int i=0; i<num_of_agents_; i++)  // check entanglement
      { 
        if (active_cases_old[i] < 2 && ent_state_begin.active_cases[i] >= 2)
        {
          // std::cout<<bold<<red<<"for agent "<<i<<"active_cases_old is "
          //          <<active_cases_old[i]<<std::endl;
          // std::cout<<bold<<red<<"for agent "<<i<<"entangle_state.active_cases is "
          //          <<current->entangle_state.active_cases[i]<<std::endl;                 
          return true;

        }
        else if (active_cases_old[i] >=2 && 
                 ent_state_begin.active_cases[i] > active_cases_old[i])
        {
          // std::cout<<bold<<red<<"for agent "<<i<<"active_cases_old is "
          //          <<active_cases_old[i]<<std::endl;
          // std::cout<<bold<<red<<"for agent "<<i<<"entangle_state.active_cases is "
          //          <<current->entangle_state.active_cases[i]<<std::endl;           
          return true;
        }
      }

      eu::updateBendPts(ent_state_begin, pkplus1, pb_, 
                        basepoint_, staticObsRep_, num_of_agents_);
      active_cases_old = ent_state_begin.active_cases;
      pk = pkplus1;   
    }

    // if (eu::getTetherLength(ent_state_begin, pb_, basepoint_, pkplus1, 
    //     staticObsRep_, num_of_agents_) > cablelength_)
    // {
    //   // std::cout<<bold<<red<<"failing cable length check!"<<std::endl;       
    //   return true;
    // }

    return false;    
  }

}

bool KinodynamicSearch::entanglesWithStaticObs(KNodePtr current)
{
  //  Eigen::Matrix<double, 2, 4> P;
  //  P.row(0)= current->coeff_x;
  //  P.row(1)= current->coeff_y;

  //  int index_z;
  // if (current->index>num_pol_) index_z = num_pol_-1;
  // else index_z = current->index-1;

  // Eigen::Vector2d pk(current->coeff_x(3), current->coeff_y(3));
  // Eigen::Vector2d pkplus1;

  // for (int j=1; j<num_sample_per_interval_; j++)
  // {
  //   // std::cout<<bold<<blue<<"checking for obstacle "<<i<<" at interval "<<j<<std::endl;

  //   pkplus1 = P * sampled_time_vector_[j];
    
  //   // std::cout<<"current->beta(num_of_agents_+i) is "<<current->beta(num_of_agents_+i)<<std::endl;
  //   // std::cout<<"current->beta2(num_of_agents_+i) is "<<current->beta2(num_of_agents_+i)<<std::endl;
  //   // std::cout<<"current->beta3(num_of_agents_+i) is "<<current->beta3(num_of_agents_+i)<<std::endl;
  //   // std::cout<<"staticObsRep_[i].col(1) is "<<staticObsRep_[i].col(1)<<std::endl;
  //   // std::cout<<"staticObsRep_[i].col(0) is "<<staticObsRep_[i].col(0)<<std::endl;
  //   // std::cout<<"coeffs_z_[index_z](3) is "<<coeffs_z_[index_z](3)<<std::endl;
  //   // std::cout<<" is "<<<<std::endl;
  //   // std::cout<<" is "<<<<std::endl;
  //   // std::cout<<" is "<<<<std::endl;
  //   // if (eu::checkEntangleStaticUpdate(current->beta(num_of_agents_+i), current->beta2(num_of_agents_+i), 
  //   //                              current->beta3(num_of_agents_+i), pk, pkplus1, staticObsRep_[i].col(1), 
  //   //                              basepoint_, staticObsRep_[i].col(0), cablelength_, coeffs_z_[index_z](3))
  //   //                               == eu::ENTANGLED)  

  //   eu::entangleHSigUpdate3(num_of_agents_+i+1, current->alphas, current->betas, 
  //                           current->active_cases[num_of_agents_+i], pk, pkplus1, staticObsRep_[i].col(1), 
  //                           staticObsRep_[i].col(1), basepoint_, staticObsRep_[i].col(0));      
  //   if (current->active_cases[num_of_agents_+i]>2)                                      
  //   {
  //     // std::cout<<bold<<blue<<"not satisfying entangle constraint 2!"<<std::endl;        
  //     return true;
  //   }
  //   pk = pkplus1;           
  // }

  // pkplus1 << current->end(0), current->end(1); //check last point

  // eu::entangleHSigUpdate3(num_of_agents_+i+1, current->alphas, current->betas, 
  //                         current->active_cases[num_of_agents_+i], pk, pkplus1, staticObsRep_[i].col(1), 
  //                         staticObsRep_[i].col(1), basepoint_, staticObsRep_[i].col(0));
  // if (current->active_cases[num_of_agents_+i]>2)                                      
  // {
  //   // std::cout<<bold<<blue<<"not satisfying entangle constraint 2!"<<std::endl;        
  //   return true;
  // }        

  return false;
}

void KinodynamicSearch::expandAndAddToQueue(KNodePtr current)
{
  MyTimer timer_expand(true);


  int jx, jy;
  double tau = T_span_;

  double delta_x = ((j_max_ - j_min_) / (num_samples_ - 1));
  Eigen::Matrix<double, 6, 1> initstates = current->end;
  // std::cout << green << "Initial to be expanded state is" << initstates<<reset << std::endl;  

  Eigen::Vector2d ji;
  for (auto comb : all_combinations_)
  {
    if (node_used_num_ == node_num_max_-1)
    {
      std::cout <<bold<<red<< "run out of memory." << std::endl;
      return;
    } 

    KNodePtr neighbor = node_pool_[node_used_num_];

    neighbor->index = current->index + 1;
    neighbor->previous = current;

    jx = std::get<0>(comb);
    jy = std::get<1>(comb);

    ji << j_min_ + jx * delta_x, j_min_ + jy * delta_x;
    // std::cout << blue << "sampled inputs are " <<ji<<reset << std::endl;  


    neighbor->end.head<2>() = initstates.head<2>() + initstates.segment<2>(2)*tau +
                       initstates.tail<2>()*tau*tau/2 + ji * tau*tau*tau/6;

    neighbor->end.segment<2>(2) = initstates.segment<2>(2) + initstates.tail<2>()*tau +
                             ji * tau*tau/2;

    neighbor->end.tail<2>() = initstates.tail<2>()+ ji * tau;
    if ((neighbor->end-initstates).norm()<0.00001) continue;

    if (neighbor->end(5)>a_max_ || neighbor->end(5)<a_min_ ||
      neighbor->end(4)>a_max_ || neighbor->end(4)<a_min_ ) //check acceleration bound
    {
        // std::cout << red << "not satisfying a max condition" <<reset << std::endl;  
        // std::cout << red << "Values are"<< neighbor->end.tail(2) <<reset << std::endl;        
      continue;
    }
    // std::cout << green << "Expanded end state is" << neighbor->end<<reset << std::endl;  

    neighbor->coeff_x << ji(0)/6, initstates(4)/2, initstates(2), initstates(0);
    neighbor->coeff_y << ji(1)/6, initstates(5)/2, initstates(3), initstates(1);

    // std::cout << green << "neighbor->coeff_x is" << neighbor->coeff_x<<reset << std::endl;  
    // std::cout << green << "neighbor->coeff_y is" << neighbor->coeff_y<<reset << std::endl;  

    Eigen::Matrix<double, 2, 4> P;
    Eigen::Matrix<double, 2, 4> Q;

    P.row(0)= neighbor->coeff_x;
    P.row(1)= neighbor->coeff_y;

    Q = P * A_rest_pos_basis_inverse_; //get position control points
    Eigen::Matrix<double, 2, 3> Qv;
    Eigen::Matrix<int, 6, 1> outIndex;

    for (int i = 0; i<4; i++) { //check position control points
      if (Q(0,i) <x_min_ || Q(0,i) >x_max_ || Q(1,i) <y_min_ || Q(1,i) >y_max_||
        (Q.col(i)-basepoint_).norm()>cablelength_)
      {
        // std::cout << red << "not satisfying x y max condition" <<reset << std::endl;  
        // std::cout << red << "Values are"<< Q <<reset << std::endl;  
        goto cont;
      }

    }

    Qv= P.block<2, 3>(0, 0) * A_rest_vel_basis_inverse321_;

    for (int i = 0; i<3; i++) { //check velocity control points
      if (Qv(0,i) <v_min_ || Qv(0,i) >v_max_ || Qv(1,i) <v_min_ || Qv(1,i) >v_max_)
      {
        // std::cout << red << "not satisfying v max condition" <<reset << std::endl;  
        // std::cout << red << "Values are"<< Qv <<reset << std::endl;  
        goto cont;
      }

    }

    //check for future violation of constraints
    if (neighbor->end(4)>0 && neighbor->end(2)-0.5*neighbor->end(4)*neighbor->end(4)/j_min_>v_max_)
    {
      // std::cout << red << "not satisfying future vel max condition for x!" <<reset << std::endl;      
      continue;
    }
    else if (neighbor->end(4)<0 && neighbor->end(2)-0.5*neighbor->end(4)*neighbor->end(4)/j_max_<v_min_)
    {
      // std::cout << red << "not satisfying future vel min condition for x!" <<reset << std::endl;          
      continue;
    }
    if (neighbor->end(5)>0 && neighbor->end(3)-0.5*neighbor->end(5)*neighbor->end(5)/j_min_>v_max_)
    {
      // std::cout << red << "not satisfying future vel max condition for y!" <<reset << std::endl;          
      continue;
    }
    else if (neighbor->end(5)<0 && neighbor->end(3)-0.5*neighbor->end(5)*neighbor->end(5)/j_max_<v_min_)
    {
      // std::cout << red << "not satisfying future vel min condition for y!" <<reset << std::endl;            
      continue;
    }    

    neighbor->entangle_state = current->entangle_state;

    double arc_length;
    arc_length = 0;    
    if (enable_entangle_check_)
    {
      if (entanglesWithOtherAgents(neighbor, arc_length)) continue;
    }
    else arc_length = (neighbor->end.head<2>()-initstates.head<2>()).norm();

    // if (entanglesWithStaticObs(neighbor)) continue;


    unsigned int ix, iy, iz;
    iz = getIz(neighbor->entangle_state);
    ix = round((neighbor->end(0) - orig_.x()) / voxel_size_);
    iy = round((neighbor->end(1) - orig_.y()) / voxel_size_);

    // Eigen::Vector2i _idx;
    // _idx = Eigen::Vector2i(ix, iy);

    // pvaToIndex(neighbor->end, outIndex);
    neighbor->Q = Q;
    neighbor->g = current-> g + arc_length;
    neighbor->h = getH(neighbor) + 0.3 * (double)neighbor->entangle_state.alphas.size() +
                  1.0 * (double)neighbor->entangle_state.bendPointsIdx.size();

    KNodePtr nodeptr;
    nodeptr = expanded_nodes_.find(Eigen::Vector2i(ix, iy), iz);
    // nodeptr = expanded_nodes_.find(Eigen::Vector2i(ix, iy));
    // nodeptr = expanded_nodes_.find(outIndex);
    if (nodeptr != NULL) 
    {
      if (nodeptr->state == 1 && nodeptr->index == neighbor->index)
      {
        //update the node if the new one is better
        if (neighbor->g + bias_*neighbor->h < nodeptr->g + bias_*nodeptr->h
            && ran_trigger%2==0)
        {
          nodeptr->previous = current;
          nodeptr->g = neighbor->g;
          nodeptr->h = neighbor->h;
          nodeptr->end = neighbor->end;
          nodeptr->coeff_x = neighbor->coeff_x;
          nodeptr->coeff_y = neighbor->coeff_y;
          nodeptr->Q = Q;

          //nodeptr->entangle_state;          
        }
        ran_trigger ++;
        continue;
      }
      else
      {
        continue;
      }
    }
    // std::cout << green << neighbor->qi.transpose() << " cost=" << neighbor->g + bias_ * neighbor->h << reset <<
    // std::endl;
    neighbor-> state = 1;
    openList_.push(neighbor);
    // std::cout << red << "pushing into open list!" <<reset << std::endl;      
    expanded_nodes_.insert(Eigen::Vector2i(ix, iy), iz, neighbor);
    // expanded_nodes_.insert(Eigen::Vector2i(ix, iy), neighbor);
    node_used_num_ += 1;   
    // expanded_nodes_.insert(outIndex, neighbor);
    // std::cout<<" "<<node_used_num_<<" ";

    cont:;
  }

}

void KinodynamicSearch::pvaToIndex(Eigen::Matrix<double, 6, 1>& input, Eigen::Matrix<int, 6, 1>& output)
{
  output(0) = round((input(0) - orig_.x()) / voxel_size_);
  output(1) = round((input(1) - orig_.y()) / voxel_size_);
  output(2) = round((input(2) - orig_.x()) / voxel_size_);
  output(3) = round((input(3) - orig_.y()) / voxel_size_);
  output(4) = round((input(4) - orig_.x()) / voxel_size_);
  output(5) = round((input(5) - orig_.y()) / voxel_size_);  
}

void KinodynamicSearch::expandAndAddToQueue(const Eigen::Matrix<double, 6, 1>& initstates)
{
  MyTimer timer_expand(true);

  double tau = T_span_;
  int jx, jy;

  double delta_x = ((j_max_ - j_min_) / (num_samples_ - 1));
  // double delta_y = ((constraint_yU - constraint_yL) / (num_samples_y_ - 1));
    std::cout << green << "Initial to be expanded state is" << initstates<<reset << std::endl;  

  // MyTimer timer_forLoop(true);
  Eigen::Vector2d ji;
  for (auto comb : all_combinations_)
  {
    KNodePtr neighbor = node_pool_[node_used_num_];

    neighbor->index = 1;
    neighbor->previous = NULL;

    jx = std::get<0>(comb);
    jy = std::get<1>(comb);

    ji << j_min_ + jx * delta_x, j_min_ + jy * delta_x;
    // std::cout << blue << "sampled inputs are " <<ji<<reset << std::endl;  

    neighbor->end.head<2>() = initstates.head<2>() + initstates.segment<2>(2)*tau +
                       initstates.tail<2>()*tau*tau/2 + ji * tau*tau*tau/6;

    neighbor->end.segment<2>(2) = initstates.segment<2>(2) + initstates.tail<2>()*tau +
                             ji * tau*tau/2;

    neighbor->end.tail<2>() = initstates.tail<2>()+ ji * tau;
    if ((neighbor->end-initstates).norm()<0.00001) continue;

    if (neighbor->end(5)>a_max_ || neighbor->end(5)<a_min_ ||
      neighbor->end(4)>a_max_ || neighbor->end(4)<a_min_ ) //check acceleration bound
    {
        // std::cout << red << "not satisfying a max condition" <<reset << std::endl;  
        // std::cout << red << "Values are"<< neighbor->end.tail(2) <<reset << std::endl;           
      continue;
    }
    // std::cout << green << "Expanded end state is" << neighbor->end<<reset << std::endl;  

    neighbor->coeff_x << ji(0)/6, initstates(4)/2, initstates(2), initstates(0);
    neighbor->coeff_y << ji(1)/6, initstates(5)/2, initstates(3), initstates(1);

    Eigen::Matrix<double, 2, 4> P;
    Eigen::Matrix<double, 2, 4> Q;

    P.row(0)= neighbor->coeff_x;
    P.row(1)= neighbor->coeff_y;

    Q = P * A_rest_pos_basis_inverse_; //get position control points
    Eigen::Matrix<double, 2, 3> Qv;
    // Eigen::Matrix<int, 6, 1> outIndex;

    for (int i = 0; i<4; i++) { //check position control points
      if (Q(0,i) <x_min_ || Q(0,i) >x_max_ || Q(1,i) <y_min_ || Q(1,i) >y_max_||
        (Q.col(i)-basepoint_).norm()>cablelength_)
      {
        // std::cout << red << "not satisfying x y max condition" <<reset << std::endl;  
        // std::cout << red << "Values are"<< Q <<reset << std::endl;          
        goto cont;
      }
    }

    // Qv= P.block(0, 0, 2, 3) * A_rest_vel_basis_inverse321_;
    // for (int i = 0; i<3; i++) { //check velocity control points
    //   if (Qv(0,i) <v_min_ || Qv(0,i) >v_max_ || Qv(1,i) <v_min_ || Qv(1,i) >v_max_)
    //   {
    //     std::cout << red << "not satisfying v max condition" <<reset << std::endl;  
    //     std::cout << red << "Values are"<< Qv <<reset << std::endl;          
    //     goto cont;
    //   }

    // }

    // for initial sampling we do not check velocity points because if the initial velocity is close to max,
    // control points methods may report exceeding the limit
    // instead, we check for future violation of velocity constraints
    if (neighbor->end(4)>0 && neighbor->end(2)-0.5*neighbor->end(4)*neighbor->end(4)/j_min_>v_max_)
    {
      // std::cout << red << "not satisfying future vel max condition for x!" <<reset << std::endl;      
      continue;
    }
    else if (neighbor->end(4)<0 && neighbor->end(2)-0.5*neighbor->end(4)*neighbor->end(4)/j_max_<v_min_)
    {
      // std::cout << red << "not satisfying future vel min condition for x!" <<reset << std::endl;          
      continue;
    }
    if (neighbor->end(5)>0 && neighbor->end(3)-0.5*neighbor->end(5)*neighbor->end(5)/j_min_>v_max_)
    {
      // std::cout << red << "not satisfying future vel max condition for y!" <<reset << std::endl;          
      continue;
    }
    else if (neighbor->end(5)<0 && neighbor->end(3)-0.5*neighbor->end(5)*neighbor->end(5)/j_max_<v_min_)
    {
      continue;
    }

    neighbor->entangle_state = entangle_state_init_;
    // neighbor->ent_state = ent_state_init_;
    double arc_length;
    arc_length = 0;
    if (enable_entangle_check_) 
    {
      if (entanglesWithOtherAgents(neighbor, arc_length))
      {
        // std::cout << red << "not satisfying entanglesWithOtherAgents condition!" <<reset << std::endl;                  
        continue;        
      } 
    }
    else arc_length = (neighbor->end.head<2>()-initstates.head<2>()).norm();

    // pvaToIndex(neighbor->end, outIndex);

    neighbor->Q = Q;
    neighbor->g = arc_length;//(neighbor->end.head<2>()-initstates.head<2>()).norm();
    neighbor->h = getH(neighbor) + 0.3 * (double)neighbor->entangle_state.alphas.size() +
                  1.0 * (double)neighbor->entangle_state.bendPointsIdx.size();

    // std::cout << green << neighbor->qi.transpose() << " cost=" << neighbor->g + bias_ * neighbor->h << reset <<
    // std::endl;
    neighbor->state = 1;          
    openList_.push(neighbor);
    // std::cout << red << "pushing into open list!" <<reset << std::endl;  
    unsigned int ix, iy, iz;
    iz = getIz(neighbor->entangle_state);
    ix = round((neighbor->end(0) - orig_.x()) / voxel_size_);
    iy = round((neighbor->end(1) - orig_.y()) / voxel_size_);

    // for (int i=0; i<num_of_agents_; i++)
    // {
    //   iz += (i+1)*neighbor->entangle_state.active_cases[i]*3;
    // }

    expanded_nodes_.insert(Eigen::Vector2i(ix, iy), iz, neighbor);
    // expanded_nodes_.insert(Eigen::Vector2i(ix, iy), neighbor);
    node_used_num_ += 1;
    // expanded_nodes_.insert(outIndex, neighbor);

    cont:;   
  }

}

bool KinodynamicSearch::collidesWithObstaclesGivenVertexes(const Eigen::Matrix<double, 3, 4>& last4Cps, int index_lastCP)
{
  // MyTimer timer_function(true);
  // std::cout << "In collidesWithObstaclesGivenVertexes, index_lastCP= " << index_lastCP << std::endl;

  int interval = index_lastCP - 3;

  Eigen::Matrix<double, 3, 4> last4Cps_new_basis = transformBSpline2otherBasis(last4Cps, interval);

  bool satisfies_LP = true;
  Eigen::Vector3d n_i;
  double d_i;

  for (int obst_index = 0; obst_index < num_of_obst_; obst_index++)
  {
    mt::ConvexHullsOfCurves_Std dummyhull;
    satisfies_LP = separator_solver_->solveModel(n_i, d_i, dummyhull[obst_index][interval], last4Cps_new_basis);
    if (satisfies_LP == false)
    {
      goto exit;
    }
  }

exit:

  return (!satisfies_LP);
}

void KinodynamicSearch::clearProcess()
{
  // expanded_valid_nodes_.clear();
  debug_goal_reached_ = 0;
  debug_setup_complete_ = 0;
  expanded_nodes_.clear();

  // stores the closest node found
  closest_dist_so_far_ = std::numeric_limits<double>::max();
  smallest_dist_times_index_ = std::numeric_limits<double>::max();
  closest_result_so_far_ptr_ = NULL;
  best_node_ptr_ = NULL;

  // stores the closest node found that has full length (i.e. its index == (N_ - 2) )
  complete_closest_dist_so_far_ = std::numeric_limits<double>::max();
  complete_closest_result_so_far_ptr_ = NULL;

  result_.clear();
  coeffs_z_.clear();

  CompareCost comparer;
  comparer.bias = bias_;
  std::priority_queue<KNodePtr, std::vector<KNodePtr>, CompareCost> openList_tmp(comparer);
  openList_.swap(openList_tmp);  

  // for (int i = 0; i < node_used_num_; i++)
  // {
  //   KNodePtr node = node_pool_[i];
  //   node->previous = NULL;
  // }
  node_used_num_ = 0;  

  closest_ptr_switcher++;
  if (closest_ptr_switcher>=100) closest_ptr_switcher = 0;

  ignore_first_col_counter_++;
  if (ignore_first_col_counter_ == 100)
  {
    ignore_first_col_counter_ = 0;
    ignore_first_col_ = false;
  }
  else ignore_first_col_ = false;

  runtime_this_round_ = 0.0;
  time_spent_contact_pt_ = 0.0;
  ran_trigger = 0;

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  shuffle(all_combinations_.begin(), all_combinations_.end(), std::default_random_engine(seed));
}
//this is a conservative but (hopefully) fast collision check using control points
bool KinodynamicSearch::collidesWithObstacles2d(const Eigen::Matrix<double, 2, 4>& Cps, int index)
{
  // MyTimer timer_function(true);
  // std::cout << "In collidesWithObstaclesGivenVertexes, index_lastCP= " << index_lastCP << std::endl;

  //consider obstacles to be static after the short planning horizon

  if (index>num_pol_) index = num_pol_; 

  for (int obst_index = 0; obst_index < num_of_obst_; obst_index++)
  {
    // std::cout<< bold<<red<<"checking collisions"<<std::endl;
    // std::cout<< bold<<red<<"Size of hulls_ is "<<hulls_.size()<<std::endl;
    // std::cout<< bold<<red<<"Size of hulls_ for an agent is "<<hulls_[obst_index].size()<<std::endl;

    double d1=(Cps.col(0)-hulls_[obst_index][index-1].col(0)).norm();
    for (int i=0; i<3; i++)
    {
      d1=d1-(Cps.col(i)-Cps.col(i+1)).norm()*safe_factor_;
      if (d1<0) 
      {
        // std::cout<<bold<<blue<<"Cps"<<std::endl;
        // std::cout<<bold<<blue<<Cps<<std::endl;     
        // std::cout<<bold<<blue<<"hulls_[obst_index][index-1]"<<std::endl;
        // std::cout<<bold<<blue<<hulls_[obst_index][index-1]<<std::endl;        
        goto exit;
      }
    }

    for (int i=0; i<hulls_[obst_index][index-1].cols()-1; i++)
    {
      d1=d1-(hulls_[obst_index][index-1].col(i)-hulls_[obst_index][index-1].col(i+1)).norm()*safe_factor_;
      if (d1<0) 
      {
        // std::cout<<bold<<blue<<"hulls_[obst_index][index-1]"<<std::endl;
        // std::cout<<bold<<blue<<hulls_[obst_index][index-1]<<std::endl;
        goto exit;
      }  
    }
  }
  return false;

exit:
  // std::cout<<bold<<blue<<"not satisfying collision constraint!"<<std::endl;
  return true;
}

//this is a not conservative but slower collision check using control points
bool KinodynamicSearch::collidesWithObstacles2dSolve(const Eigen::Matrix<double, 2, 4>& Cps, int index)
{
  // MyTimer timer_function(true);
  // std::cout << "In collidesWithObstaclesGivenVertexes, index_lastCP= " << index_lastCP << std::endl;

  //consider obstacles to be static after the short planning horizon
  if (index>num_pol_) index = num_pol_; 

  for (int obst_index = 0; obst_index < num_of_obst_; obst_index++)
  {
    // early skip if it satisfy a simple but conservative collision check
    // double d1=(Cps.col(0)-hulls_[obst_index][index-1].col(0)).norm();
    // for (int i=0; i<3; i++)
    // {
    //   d1=d1-(Cps.col(i)-Cps.col(i+1)).norm();
    // }

    // for (int i=0; i<hulls_[obst_index][index-1].cols()-1; i++)
    // {
    //   d1=d1-(hulls_[obst_index][index-1].col(i)-hulls_[obst_index][index-1].col(i+1)).norm();
    // }

    // if (d1>0) continue;

    // MyTimer timer_function(true);

    // bool satisfies_LP = separator_solver_->solveModel(hulls_[obst_index][index-1], Cps);
    // std::cout<<green<<"time ElapsedUs for solvemodel is "<< timer_function.ElapsedUs()<<std::endl;
    // timer_function.Reset();
    // satisfies_LP = cu::intersect(hulls_[obst_index][index-1], Cps);
    // std::cout<<red<<"time ElapsedUs for intersect is "<< timer_function.ElapsedUs()<<std::endl;
    bool gjk_collision = gjk::collision(hulls_[obst_index][index-1], Cps);
    // std::cout<<green<<"time ElapsedUs for gjk collision check is "<< timer_function.ElapsedUs()<<std::endl;

    // if (gjk_collision==satisfies_LP)
    // {
    //   std::cout<<red<<"gjk collision check returns different result as separator! "<<std::endl;
    //   std::cout<<red<<"gjk collision check: "<<gjk_collision<<std::endl;

    // }
    // std::cout << "intersection returns " << satisfies_LP << std::endl;

    if (gjk_collision == true)
    {
      // std::cout<<bold<<red<<"collides with agent "<< obst_index+1<<" !"<<std::endl;
      goto exit;
    }    

  }

  for (int obst_index = 0; obst_index < num_of_static_obst_; obst_index++)
  {

    bool gjk_collision = gjk::collision(convexHullOfStaticObs_[obst_index], Cps);

    if (gjk_collision == true)
    {
      // std::cout<<bold<<red<<"collides with static obst "<< obst_index+1<<" !"<<std::endl;
      goto exit;
    }    
  }  
  return false;

exit:
  // std::cout<<bold<<blue<<"not satisfying collision constraint!"<<std::endl;
  return true;
}


bool KinodynamicSearch::collidesWithBases2d(const Eigen::Matrix<double, 2, 4>& Cps)
{
  if (!enable_entangle_check_) return false;

  double radius = 0.7;
  double safe_dist = T_span_*v_max_*2;
  for (int _index = 0; _index < num_of_agents_; _index++)
  {
    // std::cout<< bold<<red<<"checking collisions"<<std::endl;
    // std::cout<< bold<<red<<"Size of hulls_ is "<<hulls_.size()<<std::endl;
    // std::cout<< bold<<red<<"Size of hulls_ for an agent is "<<hulls_[obst_index].size()<<std::endl;
    if (_index == id_-1) continue; // no collision check with its own base

    double d1=(Cps.col(0)- pb_[_index]).norm();
    if (d1>safe_dist) continue;
    // for (int i=0; i<3; i++)
    // {
    //   d1=d1-(Cps.col(i)-Cps.col(i+1)).norm()*safe_factor_;
    //   // if (d1<0)
    //   // {
    //   //   // std::cout<<bold<<blue<<"collides with base of agent "<<_index+1<<std::endl;        
    //   //   // std::cout<<bold<<blue<<"Cps"<<std::endl;
    //   //   // std::cout<<bold<<blue<<Cps<<std::endl;           
    //   //   goto exit;
    //   // } 
    // }
    // if (d1>0) continue;
    Eigen::Matrix<double, 2, 4> cps_base;
    cps_base << pb_[_index](0) + radius, pb_[_index](0) + radius, pb_[_index](0) - radius, pb_[_index](0) - radius, 
                pb_[_index](1) + radius, pb_[_index](1) - radius, pb_[_index](1) - radius, pb_[_index](1) + radius;
 
    // bool satisfies_LP = separator_solver_->solveModel(cps_base, Cps);
    bool gjk_collision = gjk::collision(cps_base, Cps);

    if (gjk_collision == true)
    {
      goto exit;
    }                    
  }
  return false;

exit:
  // std::cout<<bold<<blue<<"not satisfying base collision constraint!"<<std::endl;
  return true;
}

bool KinodynamicSearch::run(std::vector<Eigen::Vector3d>& result, int& status)
{
  std::cout << "[A*] Running..." << std::endl;

  MyTimer timer_astar(true);

  expandAndAddToQueue(initial_);

  int RUNTIME_REACHED = 0;
  int GOAL_REACHED = 1;
  int EMPTY_OPENLIST = 2;
  KNodePtr current_ptr;

  while (openList_.size() > 0)
  {
    // std::cout << "[A*] current open list size is "<<openList_.size() << std::endl;
    // Check if runtime is over
    if (timer_astar.ElapsedMs() > (max_runtime_ * 1000))
    {
      std::cout << "[A*] Max Runtime was reached" << std::endl;
      std::cout << "[A*] Run time is " << timer_astar.ElapsedMs()<<std::endl;      
      status = RUNTIME_REACHED;
      goto exitloop;
    }

    current_ptr = openList_.top();  // copy the first element onto current_ptr
    openList_.pop();                 // remove it from the list
    current_ptr-> state = -1;

    double dist = (current_ptr->end.head<2>() - goal2d_).norm();
    double dist_from_initial = (current_ptr->end.head<2>() - initial_.head(2)).norm();

    /////////////////////
    // MyTimer timer_collision_check(true);
    // // bool collides = collidesWithObstacles(*current_ptr);
    bool constraintViolated = false;
    if (current_ptr->index > 1 || current_ptr->index == 1 && !ignore_first_col_) 
    //only check collision if it is not the initial segment
    //this is to prevent two agents end up in a dead zone
    {
      constraintViolated = collidesWithObstacles2dSolve(current_ptr->Q, current_ptr->index);
    }

    // // std::cout << "collision check took " << timer_collision_check << std::endl;
    if (constraintViolated) continue;

    constraintViolated = collidesWithBases2d(current_ptr->Q);
    if (constraintViolated) continue;

    // constraintViolated = entanglesWithOtherAgents(current_ptr);
    // if (constraintViolated) continue;

    // expanded_valid_nodes_.push_back(current_ptr);

    // if ((dist < goal_size_) && dist < std::max(0.0, complete_closest_dist_so_far_ - 1e-4))  // the 1e-4 is to avoid numerical issues of paths
    //                                                                  // essentially with the same dist to the goal. In
    //                                                                  // those cases, this gives priority to the
    //                                                                  // trajectories found first
    // {
    //   complete_closest_dist_so_far_ = dist;
    //   complete_closest_result_so_far_ptr_ = current_ptr;
    // }

    // if (dist < closest_dist_so_far_ && current_ptr->beta3.head(num_of_agents_).isZero())

    bool valid_endpoint = true;
    for (int i=0; i<num_of_agents_; i++)
    {
      if (current_ptr->entangle_state.active_cases[i]>1)
      {
        valid_endpoint = false;
        break;
      }
    }
    if (current_ptr->index == 1 && ignore_first_col_ && 
        collidesWithObstacles2dSolve(current_ptr->Q, current_ptr->index))
    {
      valid_endpoint = false;
    }

    double dist_to_compare;
    if (!goal_occupied_) dist_to_compare = dist_from_initial;
    else dist_to_compare = dist*dist;

    double dist_times_index = dist_to_compare * (double)current_ptr->index;
    if (dist_times_index < smallest_dist_times_index_ && valid_endpoint) //while safe, try to stay near the start
    {
      smallest_dist_times_index_ = dist_times_index;
      closest_dist_so_far_ = dist_to_compare;
      closest_result_so_far_ptr_ = current_ptr;
    }

    // check if we are already in the goal
    if ((dist < goal_size_))
    {
      // if (!current_ptr->beta3.head(num_of_agents_).isZero()) continue; //not a valid path if it is risk entangling
      if (valid_endpoint) //not a valid path if it is risk entangling
      {
        //publish the goal status for debugging
        // printEntangleState(current_ptr->entangle_state);

        runtime_this_round_ = (double)timer_astar.ElapsedUs()/1000.0;
        std::cout << "[A*] Goal was reached!" << std::endl;
        std::cout << bold << green <<"[A*] solving time is" << runtime_this_round_ << "ms !!"<<std::endl;
        std::cout << bold << green <<"[A*] num of expansion is" << node_used_num_ << " !!"<<std::endl;
        status = GOAL_REACHED;
        goto exitloop;
      }
    }

    expandAndAddToQueue(current_ptr);
  }

  std::cout << "[A*] openList_ is empty" << std::endl;
  std::cout << "[A*] Run time is " << timer_astar.ElapsedMs()<<std::endl;      
  status = EMPTY_OPENLIST;
  goto exitloop;

exitloop:

  // std::cout << "status= " << status << std::endl;
  // std::cout << "expanded_nodes_.size()= " << expanded_nodes_.size() << std::endl;
  // std::cout << "complete_closest_dist_so_far_= " << complete_closest_dist_so_far_ << std::endl;

  KNodePtr best_node_ptr = NULL;

  bool have_a_solution = (closest_result_so_far_ptr_ != NULL);

  if (status == GOAL_REACHED)
  {
    std::cout << "[A*] choosing current_ptr as solution" << std::endl;
    best_node_ptr = current_ptr;
  }
  else if (status == RUNTIME_REACHED || status == EMPTY_OPENLIST)
  {
    if (have_a_solution && use_not_reaching_solution_)
    {
      // if (closest_dist_so_far_<voxel_size_) 
      // {
      //   while (closest_result_so_far_ptr_->index >1) //only execute the first step of the safe trajectory
      //   {
      //     closest_result_so_far_ptr_ = closest_result_so_far_ptr_->previous;
      //   }

      //   closest_result_so_far_ptr_->coeff_x << 0.0, 0.0, 0.0, initial_(0); 
      //   closest_result_so_far_ptr_->coeff_y << 0.0, 0.0, 0.0, initial_(1); 
      //   best_node_ptr = closest_result_so_far_ptr_;
      // }
      // else
      // {
        best_node_ptr = closest_result_so_far_ptr_;        
      // }
      std::cout << "[A*] choosing closest safe path as solution" << std::endl;
      std::cout << bold << blue << "closest result index is " << closest_result_so_far_ptr_->index << reset
                << std::endl;
      std::cout << bold << blue << "closest_dist_so_far_= " << closest_dist_so_far_ << reset
                << std::endl;
      std::cout << bold << green <<"[A*] num of expansion is" << node_used_num_ << " !!"<<std::endl;

    }
    else
    {
      std::cout << "[A*] not finding a feasible solution!" << std::endl;
      // best_node_ptr = closest_result_so_far_ptr_;
      if (status == RUNTIME_REACHED)
      {
        std::cout <<bold<<red<< "[A*] RUN TIME REACHED!" << std::endl;
      }
      std::cout << "[A*] num of nodes expanded is " <<node_used_num_<<std::endl;
      return false;
    }
  }
  else
  {
    std::cout << red << "This state should never occur" << reset << std::endl;
    abort();
    return false;
  }

  // Fill (until we arrive to N_-2), with the same qi
  // Note that, by doing this, it's not guaranteed feasibility wrt a dynamic obstacle
  // and hence the need of the function checkFeasAndFillND()

  bool path_found_is_not_complete = (best_node_ptr->index < num_pol_);

  if (path_found_is_not_complete)
  {

  }

  // recoverPath(best_node_ptr);
  recoverPwpOut(best_node_ptr);
  recoverEntStateVector(best_node_ptr);

  result = result_;
  best_node_ptr_ = best_node_ptr;

  return true;
}

bool KinodynamicSearch::runonce(std::vector<Eigen::Vector3d>& tmp_path, int &status)
{
  if (debug_goal_reached_)
  {
    status=1;
    return true;
  }

  // std::cout << red << "Running a star" << reset << std::endl;  

  if (!debug_setup_complete_)
  {
    debug_setup_complete_=true;

    std::cout << "[A*] Running..." << std::endl;

    expandAndAddToQueue(initial_);
  }

  KNodePtr current_ptr;
  std::cout << "[A*] current open list size is "<<openList_.size() << std::endl;

  if (openList_.size() == 0)
  {
    status = -1;
  }
  else
  {
    current_ptr = new KNode;

    current_ptr = openList_.top();  // copy the first element onto (*current_ptr)
    openList_.pop();                 // remove it from the list

    double dist = (current_ptr->end.head<2>() - goal2d_).norm();
    double dist_from_initial = (current_ptr->end.head<2>() - initial_.head(2)).norm();

    /////////////////////
    // MyTimer timer_collision_check(true);
    bool constraintViolated = collidesWithObstacles2dSolve(current_ptr->Q, current_ptr->index);
    // // std::cout << "collision check took " << timer_collision_check << std::endl;
    if (constraintViolated) goto endcurrent;

    constraintViolated = collidesWithBases2d(current_ptr->Q);
    if (constraintViolated) goto endcurrent;
    
    // constraintViolated = entanglesWithOtherAgents(current_ptr);
    // if (constraintViolated) goto endcurrent;

    // expanded_valid_nodes_.push_back(current_ptr);

    //publish the goal status for debugging
    printEntangleState(current_ptr->entangle_state);

    bool valid_endpoint;
    valid_endpoint = true;

    for (int i=0; i<num_of_agents_; i++)
    {
      if (current_ptr->entangle_state.active_cases[i]>1)
      {
        valid_endpoint = false;
        break;
      }
    }

    double dist_to_compare;
    if (!goal_occupied_) dist_to_compare = dist_from_initial;
    else dist_to_compare = dist;

    if (dist_to_compare < closest_dist_so_far_ && valid_endpoint)
    {
      closest_dist_so_far_ = dist_to_compare;
      closest_result_so_far_ptr_ = current_ptr;
    }

    // check if we are already in the goal
    if ((dist < goal_size_))
    {
      if (valid_endpoint) //not a valid path if it is risk entangling 
      {
        //publish the goal status for debugging
        printEntangleState(current_ptr->entangle_state);

        std::cout << "[A*] Goal was reached!" << std::endl;
        status = 1;
        debug_goal_reached_=true;
        return true; //exitloop not used for now
        goto exitloop;
      }
    }

    expandAndAddToQueue(current_ptr);

    endcurrent:
    status = 2; //ongoing
    getCurrentSamplePath(current_ptr, tmp_path);    
  }


  return true;

exitloop:


  KNodePtr best_node_ptr = NULL;

  bool have_a_solution = (closest_result_so_far_ptr_ != NULL);

  if (status == 1)
  {
    std::cout << "[A*] choosing current_ptr as solution" << std::endl;
    best_node_ptr = current_ptr;
  }
  else if (status == -1)
  {

    if (have_a_solution)
    {
      std::cout << "[A*] choosing closest safe path as solution" << std::endl;
      best_node_ptr = closest_result_so_far_ptr_;
      std::cout << bold << blue << "closest_dist_so_far_= " << closest_dist_so_far_ << reset
                << std::endl;
    }
    else
    {
      std::cout << "[A*] not even finding a safe path" << std::endl;
      return false;
    }
  }
  else
  {
    std::cout << red << "This state should never occur lol" << reset << std::endl;
    abort();
    return false;
  }

  bool path_found_is_not_complete = (best_node_ptr->index < num_pol_);

  if (path_found_is_not_complete)
  {

  }

  recoverPath(best_node_ptr);  // saved in result_

  return true;
}

void KinodynamicSearch::getRuntime(double& runtime_this_round, double& time_spent_contact_pt,
                                  int& node_used_num)
{
  runtime_this_round = runtime_this_round_;
  time_spent_contact_pt = time_spent_contact_pt_;
  node_used_num = node_used_num_;
}

void KinodynamicSearch::printEntangleState(eu::ent_state& entangle_state)
{
    //publish the goal status for debugging
  for (int i=0; i<entangle_state.alphas.size(); i++)
  {
    std::cout<<bold<<blue<<"CURRENT STATE alphas "<<i<<"th entry is" 
              <<entangle_state.alphas[i]<<std::endl;
  }
  for (int i=0; i<entangle_state.active_cases.size(); i++)
  {
    std::cout<<bold<<blue<<"CURRENT STATE active cases for "<<i<<"th agent is" 
             <<entangle_state.active_cases[i]<<std::endl;
  }
  // for (int i=0; i<entangle_state.bendPointsIdx.size(); i++)
  // {
  //   std::cout<<bold<<blue<<"CURRENT STATE "<<i<<"th bendPointsIdx is " <<
  //             entangle_state.alphas[entangle_state.bendPointsIdx[i]](0)
  //             <<std::endl;
  // }
}

unsigned int KinodynamicSearch::getIz(eu::ent_state& entangle_state)
{ 
    int iz = 0;
    for (int i=0; i<entangle_state.alphas.size(); i++)
    {
      iz += (i+1)*power_int(entangle_state.alphas[i](0), entangle_state.alphas[i](1));
    }
    return iz;
}

unsigned int KinodynamicSearch::power_int(unsigned int base, unsigned int exponent)
{
  if (exponent == 0) { return 1;    }
  if (base < 2)      { return base; }

  unsigned int result = 1;

  for (unsigned int term = base; ; term = term * term)
  { 
    if (exponent % 2 != 0) { result *= term; }
    exponent /= 2;
    if (exponent == 0)     { break; }
  }

  return result;
}