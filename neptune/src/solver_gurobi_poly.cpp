/* ----------------------------------------------------------------------------
 * Copyright 2021, Cao Muqing
 * Nanyang Technological University
 * All Rights Reserved
 * Authors: Cao Muqing, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "solver_gurobi_poly.hpp"
#include "termcolor.hpp"
#include "bspline_utils.hpp"
#include "ros/ros.h"
#include "solver_gurobi_utils.hpp"

#include <decomp_util/ellipsoid_decomp.h>  //For Polyhedron definition
#include <unsupported/Eigen/Splines>
#include <iostream>
#include <list>
#include <random>
#include <iostream>
#include <vector>
using namespace termcolor;


PolySolverGurobi::PolySolverGurobi(int num_pol, int deg_pol, int id, double T_span, 
                                   std::vector<Eigen::Vector2d> pb, double weight_term, 
                                   double rad_term, bool use_linear_constraints)
{
  // std::cout << green << bold << "entering point 1" << reset << std::endl;
  p_ = deg_pol;
  num_pol_ = num_pol;
  rad_term_ = rad_term;
  use_linear_constraints_ = use_linear_constraints;

  double tspan4 = std::pow(T_span, 4);
  double tspan3 = std::pow(T_span, 3);
  double tspan2 = std::pow(T_span, 2);

  mt::basisConverter basis_converter; //convert from interval 0,1 to a user-defined one
  Eigen::Matrix<double, 4, 4> C_interval_convert;
  C_interval_convert << 
  1/tspan3, 0, 0, 0, 
  0, 1/tspan2, 0, 0, 
  0, 0, 1/T_span, 0,
  0, 0, 0, 1;

  Eigen::Matrix<double, 3, 3> Cv_interval_convert; //this one not checked yet
  Cv_interval_convert << 
  1/tspan2, 0, 0, 
  0, 1/T_span, 0, 
  0, 0, 1;

  std::string basis = "MINVO";
  T_span_ = T_span;
  if (basis == "MINVO")
  {
    // std::cout << green << bold << "A* is using MINVO" << reset << std::endl;
    // M_pos_bs2basis_ = basis_converter.getMinvoPosConverters(num_pol)*C_interval_convert;
    // M_vel_bs2basis_ = basis_converter.getBSplineVelConverters(num_pol)*Cv_interval_convert;  // getMinvoVelConverters TODO!!
    // basis_ = MINVO;
    A_rest_pos_basis_ = basis_converter.getArestMinvo()*C_interval_convert;
    A_rest_vel_basis_ = basis_converter.getAvelrestMinvo()*Cv_interval_convert;

  }
  else if (basis == "BEZIER")
  {
    // std::cout << green << bold << "A* is using BEZIER" << reset << std::endl;
    M_pos_bs2basis_ = basis_converter.getBezierPosConverters(num_pol);
    M_vel_bs2basis_ = basis_converter.getBSplineVelConverters(num_pol);  // getBezierVelConverters TODO!!
    // basis_ = BEZIER;
  }
  else if (basis == "B_SPLINE")
  {
    // std::cout << green << bold << "A* is using B_SPLINE" << reset << std::endl;

    M_pos_bs2basis_ = basis_converter.getBSplinePosConverters(num_pol);
    M_vel_bs2basis_ = basis_converter.getBSplineVelConverters(num_pol);

    // basis_ = B_SPLINE;
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
  // std::cout << green << bold << "entering point 4" << reset << std::endl;



  id_ = id;	
  pb_ = pb;


  basepoint_ = pb[id_-1];
  num_of_agents_ = pb.size();


  Q_v_term_weighted << 9*tspan4, 6*tspan3, 3*tspan2, 0,
                       6*tspan3, 4*tspan2, 2*T_span, 0,
                       3*tspan2, 2*T_span, 1, 0, 
                       0, 0, 0, 0;

  Q_a_term_weighted << 36*tspan2, 12*T_span, 0, 0, 
                       12*T_span, 4, 0, 0, 
                       0, 0, 0, 0,
                       0, 0, 0, 0;
  Q_v_term_weighted *= weight_term;
  Q_a_term_weighted *= weight_term;

  Q_p_term << tspan3*tspan3, tspan3*tspan2, tspan4, tspan3,
             tspan3*tspan2, tspan4, tspan3, tspan2,
             tspan4, tspan3, tspan2, T_span,
             tspan3, tspan2, T_span, 1;
  q_p_term << tspan3, tspan2, T_span, 1;
  // q_p_term *= 2;
  q_v_term << 3*tspan2, 2*T_span, 1, 0;
  q_a_term << 6*T_span, 2, 0, 0;

  weight_ = weight_term;
  separator_solver_ = new separator::Separator();  // 0.0, 0.0, 0.0

}

PolySolverGurobi::~PolySolverGurobi()
{
}

void PolySolverGurobi::setMaxValues(double x_min, double x_max, double y_min, double y_max, double z_min,
                                      double z_max, double v_max, double a_max, double j_max)
{

  maxs_.clear();
  mins_.clear();

  maxs_.push_back(x_max);
  maxs_.push_back(y_max);
  maxs_.push_back(z_max);

  mins_.push_back(x_min);
  mins_.push_back(y_min);
  mins_.push_back(z_min);

  // x_min_ = x_min;
  // x_max_ = x_max;

  // y_min_ = y_min;
  // y_max_ = y_max;

  // z_min_ = z_min;
  // z_max_ = z_max;

  v_max_ = v_max;
  v_min_ = -v_max;

  a_max_ = a_max;
  a_min_ = -a_max;

  j_max_ = j_max;
  j_min_ = -j_max;

  long_length_ = sqrt((x_max-x_min)*(x_max-x_min)+(y_max-y_min)*(y_max-y_min));

}

void PolySolverGurobi::setTetherLength(double tetherLength)
{
  cablelength_ = tetherLength;
}

void PolySolverGurobi::setMaxRuntime(double max_runtime)
{
  max_runtime_ = max_runtime;
}

void PolySolverGurobi::setInitTrajectory(mt::PieceWisePol pwp_init)
{
  pwp_init_.clear();
  pwp_out_.clear();

	pwp_init_ = pwp_init;
  pwp_out_.times = pwp_init.times;
  num_pol_init_ = pwp_init_.coeff_x.size();

  std::vector<std::string> coords = { "x", "y", "z" };

  poly_exp_.clear();

  for (int k=0; k<3; k++) //x, y, x
  {  
    std::vector<std::vector<GRBLinExpr>> poly_mat;
    for (int i=0; i<num_pol_init_; i++)
    {
      std::vector<GRBLinExpr> poly_vec;
    	for (int j=0; j<=p_; j++)
    	{
        GRBVar tmp;
        if (j!=p_)
        {  
          tmp = m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "q_" + 
                          std::to_string(i) +  std::to_string(j) + coords[k]);
        }
        else
        {
          tmp = m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "q_" + 
                          std::to_string(i) + std::to_string(j) + coords[k]);
        }
        poly_vec.push_back(GRBLinExpr(tmp));
    	}
      poly_mat.push_back(poly_vec);
    }
    poly_exp_.push_back(poly_mat);
  }

  final_pos_(0) = q_p_term.transpose() * pwp_init_.coeff_x[num_pol_init_-1];
  final_pos_(1) = q_p_term.transpose() * pwp_init_.coeff_y[num_pol_init_-1];
  final_pos_(2) = q_p_term.transpose() * pwp_init_.coeff_z[num_pol_init_-1];

  ctrlPtsInit_.clear();

  for (int i=0; i<num_pol_init_; i++)
  {
    Eigen::Matrix<double, 2, 4> P;
    Eigen::Matrix<double, 2, 4> Q;

    P.row(0)= pwp_init_.coeff_x[i];
    P.row(1)= pwp_init_.coeff_y[i];

    Q = P * A_rest_pos_basis_inverse_; //get position control points

    ctrlPtsInit_.push_back(Q);
  }
}

void PolySolverGurobi::setHulls(mt::ConvexHullsOfCurves_Std2d &hulls)
{
  hulls_.clear();
  hulls_ = hulls;

  num_of_obst_ = hulls_.size();


  // Create the variables

  collision_plane_exp_.clear();

  for (int k=0; k<num_of_obst_; k++) //
  {  
    std::vector<std::vector<GRBLinExpr>> poly_mat;
    for (int i=0; i<num_pol_init_; i++)
    {
      std::vector<GRBLinExpr> poly_vec;
      GRBVar tmp1, tmp2, tmp3;
      tmp1 = m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "a_" + 
                          std::to_string(i) +  std::to_string(k));
      tmp2 = m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "b_" + 
                          std::to_string(i) +  std::to_string(k));
      tmp3 = m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "c_" + 
                          std::to_string(i) +  std::to_string(k));

      poly_vec.push_back(GRBLinExpr(tmp1));
      poly_vec.push_back(GRBLinExpr(tmp2));
      poly_vec.push_back(GRBLinExpr(tmp3));
      poly_mat.push_back(poly_vec);
    }
    collision_plane_exp_.push_back(poly_mat);
  }
  std::cout << green << "Gurobi finish setting convex hulls" << std::endl;

}

void PolySolverGurobi::setHullsNoInflation(mt::ConvexHullsOfCurves_Std2d &hulls)
{
  hullsNoInflation_.clear();
  hullsNoInflation_ = hulls;

}

void PolySolverGurobi::setBetasVector(std::vector<std::vector<Eigen::Vector3d>>& vecOfAgents)
{
  vecOfAgentsBetas_.clear();
  vecOfAgentsBetas_ = vecOfAgents;

  // for (int i=0; i<num_of_agents_; i++)
  // {
  //   if (i==id_-1) continue;
  //   if (vecOfAgentsBetas_[i].empty()) continue;
  //   for (int j=0; j<num_pol_init_; j++)
  //   {
  //     std::cout << green << "betas for agent "<<i<<" at "<<j<<"th seg is "
  //               <<vecOfAgentsBetas_[i][j]<< std::endl;
  //   }
  // }
}

void PolySolverGurobi::setEntStateVector(std::vector<eu::ent_state>& entStateVec,
                                        std::vector<std::vector<Eigen::Vector2d>>& bendPtsForAgents)
{
  entStateVec_.clear();
  entStateVec_ = entStateVec;
  bendPtsForAgents_.clear();
  bendPtsForAgents_ = bendPtsForAgents;
}

void PolySolverGurobi::setStaticObstVert(std::vector<mt::Polygon_Std>& convexHullOfStaticObs)
{
  num_of_static_obst_ = convexHullOfStaticObs.size();
  convexHullOfStaticObs_ = convexHullOfStaticObs;
}

void PolySolverGurobi::addObjective(bool add_terminal_va_cost)
{
  // std::cout << green << "Gurobi adding objective" << std::endl;

  GRBQuadExpr cost = 0.0;

  double q_jerk = 36 * T_span_; //(3*2*1)^2*t

  for (int i = 0; i < num_pol_init_; i++) //jerk cost
  {

    cost += q_jerk * (poly_exp_[0][i][0] * poly_exp_[0][i][0] + 
                      poly_exp_[1][i][0] * poly_exp_[1][i][0] +
                      poly_exp_[2][i][0] * poly_exp_[2][i][0]);

  }

  // terminal vel and acc as cost
  if (add_terminal_va_cost)
  {
    for (int i=0; i<3; i++)
    {
      for (int k1=0; k1<3; k1++)// for terminal vel, only first 3 terms
      {
        for (int k2=0; k2<3; k2++) 
        {
          cost += Q_v_term_weighted(k1, k2) * poly_exp_[i][num_pol_init_-1][k1] * poly_exp_[i][num_pol_init_-1][k2];
        }
      }

      for (int k1=0; k1<2; k1++)
      {
        for (int k2=0; k2<2; k2++) // for terminal acc, only first 2 terms
        {
          cost += Q_a_term_weighted(k1, k2) * poly_exp_[i][num_pol_init_-1][k1] * poly_exp_[i][num_pol_init_-1][k2];
        }
      }
    }
  }

  //terminal pos as a cost
  GRBQuadExpr term_p_cost = 0.0;
  for (int i=0; i<3; i++)
  {
    for (int k1=0; k1<4; k1++)// for terminal pos
    {
      for (int k2=0; k2<4; k2++) 
      {
        term_p_cost += Q_p_term(k1, k2) * poly_exp_[i][num_pol_init_-1][k1] * poly_exp_[i][num_pol_init_-1][k2];
      }

      term_p_cost -= 2* q_p_term(k1) * poly_exp_[i][num_pol_init_-1][k1] * final_pos_(i);
    }
    term_p_cost += final_pos_(i)*final_pos_(i);
  }  
  term_p_cost = term_p_cost * weight_;

  cost += term_p_cost;
  m_.setObjective(cost, GRB_MINIMIZE);
  // std::cout << green << "Gurobi finish adding objective" << std::endl;

}

void PolySolverGurobi::addConstraints()
{
  // std::cout << green << "Gurobi adding constraints" << std::endl;

  //initial p v a constraint
  for (int k1=1; k1<4; k1++)// for init p, v, a
  {
    m_.addConstr(poly_exp_[0][0][k1] == pwp_init_.coeff_x[0](k1));
    m_.addConstr(poly_exp_[1][0][k1] == pwp_init_.coeff_y[0](k1));
    m_.addConstr(poly_exp_[2][0][k1] == pwp_init_.coeff_z[0](k1));

  }
  // std::cout << green << "Gurobi finish adding init constraints" << std::endl;

  // continuity constraint
  for (int i=0; i<num_pol_init_-1; i++)
  {
    for (int j=0; j<3; j++) // x, y, z
    {
      GRBLinExpr tmp;
      for (int k=0; k<4; k++)
      {
        tmp += q_p_term(k) * poly_exp_[j][i][k];
      }
      m_.addConstr(tmp - poly_exp_[j][i+1][3]==0);

      GRBLinExpr tmp_v;
      for (int k=0; k<3; k++)
      {
        tmp_v += q_v_term(k) * poly_exp_[j][i][k];
      }
      m_.addConstr(tmp_v - poly_exp_[j][i+1][2]==0);

      GRBLinExpr tmp_a;
      for (int k=0; k<2; k++)
      {
        tmp_a += q_a_term(k) * poly_exp_[j][i][k];
      }
      m_.addConstr(tmp_a - 2*poly_exp_[j][i+1][1]==0);      
    }
  }

  // std::cout << green << "Gurobi finish adding continuity constraints" << std::endl;

  double tspan_times_6 = T_span_ * 6;

  // we put all constraints that need the control points together, so we do not need to 
  // save global variables for control points linear expression
  for (int i=0; i<num_pol_init_; i++)
  {
    std::vector<std::vector<GRBLinExpr>> ctrl_pt_xyz;

    for (int j=0; j<3; j++) // x, y, z
    {
      std::vector<GRBLinExpr> ctrl_pt;

      for (int k=0; k<4; k++)
      {
        GRBLinExpr tmp;
        tmp = poly_exp_[j][i][0] * A_rest_pos_basis_inverse_(0, k) +
              poly_exp_[j][i][1] * A_rest_pos_basis_inverse_(1, k) +
              poly_exp_[j][i][2] * A_rest_pos_basis_inverse_(2, k) +
              poly_exp_[j][i][3] * A_rest_pos_basis_inverse_(3, k);
        m_.addConstr(tmp <=maxs_[j]);      
        m_.addConstr(tmp >=mins_[j]);  

        ctrl_pt.push_back(tmp); 
      }

      ctrl_pt_xyz.push_back(ctrl_pt);

      for (int k=0; k<3; k++)
      {
        GRBLinExpr tmp;
        tmp = poly_exp_[j][i][0] * A_rest_vel_basis_inverse321_(0, k) +
              poly_exp_[j][i][1] * A_rest_vel_basis_inverse321_(1, k) +
              poly_exp_[j][i][2] * A_rest_vel_basis_inverse321_(2, k);

        m_.addConstr(tmp <= v_max_);      
        m_.addConstr(tmp >= v_min_);      
      }

      GRBLinExpr tmp;
      tmp = poly_exp_[j][i][0] * tspan_times_6 + poly_exp_[j][i][1] * 2;
      m_.addConstr(tmp <= a_max_);      
      m_.addConstr(tmp >= a_min_);      
    }
    // std::cout << green << "Gurobi finish adding pos vel acc constraints" << std::endl;

    // adding inter-agent collision constraints
    if (use_linear_constraints_)
    {
      for (int j=0; j<num_of_obst_; j++)  
      {
        Eigen::Vector3d sp_line;
        bool solved = separator_solver_->solveModel(sp_line, hulls_[j][i], ctrlPtsInit_[i]); 
        // std::cout << green << "adding for agent " <<j<< std::endl;

        if (solved)
        {
          for (int k=0; k<4; k++) //add linear constraint
          {
            m_.addConstr(ctrl_pt_xyz[0][k] * sp_line(0) + ctrl_pt_xyz[1][k] * sp_line(1) +
                          sp_line(2) -1 <= 0);
          }          
        }
        else
        {
          std::cout << red << "separation solver not getting solution, conflict with the gjk algo!" << std::endl;
        }
      }    
    }
    else
    {
      //for collision constraint using separation line
      for (int j=0; j<num_of_obst_; j++)  
      {
        // std::cout << green << "adding for agent " <<j<< std::endl;

        for (int k=0; k<4; k++) //for our control points, add quadratic constraint
        {
          m_.addQConstr(ctrl_pt_xyz[0][k] * collision_plane_exp_[j][i][0] +
                        ctrl_pt_xyz[1][k] * collision_plane_exp_[j][i][1] +
                        collision_plane_exp_[j][i][2] <=-1e-7);
        }

        for (int k=0; k<hulls_[j][i].cols(); k++) //for agents points, add linear constraint
        {
          m_.addConstr(hulls_[j][i](0,k) * collision_plane_exp_[j][i][0] +
                       hulls_[j][i](1,k) * collision_plane_exp_[j][i][1] +
                       collision_plane_exp_[j][i][2] >= 1e-7);        
        }
      }       
    }
    // std::cout << green << "Gurobi finish adding inter-agent collision constraints" << std::endl;

    // adding base collision constraint
    double base_radius = 0.7;
    for (int j=0; j<num_of_agents_; j++)
    {
      bool close_to_base = false;
      for (int k=0; k<4; k++) //add linear constraint
      {
        if ((ctrlPtsInit_[i].col(k)-pb_[j]).norm()<base_radius*3)
        {
          close_to_base = true;
          break;
        }
      }
      if (close_to_base)
      {
        Eigen::Matrix<double, 2, 4> base_hull;
        base_hull << pb_[j](0) + base_radius, pb_[j](0) + base_radius,
                     pb_[j](0) - base_radius, pb_[j](0) - base_radius,
                     pb_[j](1) + base_radius, pb_[j](1) - base_radius,
                     pb_[j](1) + base_radius, pb_[j](1) - base_radius;

        Eigen::Vector3d sp_line;
        bool solved = separator_solver_->solveModel(sp_line, base_hull, ctrlPtsInit_[i]);   
        if (solved)
        {
          for (int k=0; k<4; k++) //add linear constraint
          {
            m_.addConstr(ctrl_pt_xyz[0][k] * sp_line(0) + ctrl_pt_xyz[1][k] * sp_line(1) +
                          sp_line(2) -1 <= 0);
          } 
        }  
      }              
    }

    // adding static obstacle collision constraint
    for (int j=0; j<num_of_static_obst_; j++)
    {
      bool close_to_static_obs = false;
      double dist = (ctrlPtsInit_[i].col(0) - convexHullOfStaticObs_[j].col(0)).norm();
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
      for (int k=0; k<convexHullOfStaticObs_[j].cols()-1; k++)
      {
        dist -= (convexHullOfStaticObs_[j].col(k+1) - convexHullOfStaticObs_[j].col(k)).norm();
        if (dist < 0)
        {
          close_to_static_obs = true;
          break;
        }          
      }
      if (close_to_static_obs)
      {
        if (use_linear_constraints_)
        {
          Eigen::Vector3d sp_line;
          bool solved = separator_solver_->solveModel(sp_line, convexHullOfStaticObs_[j], ctrlPtsInit_[i]);  
          if (solved)
          {
            for (int k=0; k<4; k++) //add linear constraint
            {
              m_.addConstr(ctrl_pt_xyz[0][k] * sp_line(0) + ctrl_pt_xyz[1][k] * sp_line(1) +
                            sp_line(2) -1 <= 0);
            } 
          }            
        }
        else
        {
          GRBLinExpr var_a, var_b, var_c;
          var_a = GRBLinExpr(m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "sta_a_" + 
                        std::to_string(i) +  std::to_string(j) ));
          var_b = GRBLinExpr(m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "sta_b_" + 
                        std::to_string(i) +  std::to_string(j) ));
          var_c = GRBLinExpr(m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "sta_c_" + 
                        std::to_string(i) +  std::to_string(j) ));            
          for (int k=0; k<4; k++) //for our control points, add quadratic constraint
          {
            m_.addQConstr(ctrl_pt_xyz[0][k] * var_a + ctrl_pt_xyz[1][k] * var_b + var_c <=-1e-7);
          }

          for (int k=0; k<convexHullOfStaticObs_[j].cols(); k++) //for obs points, add linear constraint
          {
            m_.addConstr(convexHullOfStaticObs_[j](0,k) * var_a + 
                         convexHullOfStaticObs_[j](1,k) * var_b + var_c >=1e-7);        
          }
        }                
      }
    } 

    // std::cout << green << "Gurobi finish adding static obsts/bases constraints" << std::endl;

    //adding non-entangling constraint using separation line
    for (int j=0; j<num_of_agents_; j++)
    {
      if (j==id_-1) continue;

      if (entStateVec_[i].active_cases[j] == 1) 
      {
        int case_id = 0;
        for (int jj=0; jj<entStateVec_[i].alphas.size(); jj++)
        {
          if (entStateVec_[i].alphas[jj](0) == j+1) case_id = entStateVec_[i].alphas[jj](1);
        }
        if (case_id == 0) continue; //not our concern atm
        for (int k=1; k<bendPtsForAgents_[j].size()+1; k++)
        {
          if (k==case_id) continue;
          addEntangleConstraintForIJCase(j, i, k, ctrl_pt_xyz);
        }
      }
      // std::cout << green << "adding entanglement constraint for agent "<< j<< std::endl;



    }

    //adding tether length constraint, not quite useful atm
    // double pbxsq = basepoint_(0) * basepoint_(0);
    // double pbysq = basepoint_(1) * basepoint_(1);
    // double tethersq = cablelength_ * cablelength_;
    // for (int k=0; k<4; k++) //for our control points, add quadratic constraint
    // {
    //   GRBQuadExpr qd1 = ctrl_pt_xyz[0][k] * ctrl_pt_xyz[0][k] - 2 * basepoint_(0) * ctrl_pt_xyz[0][k]
    //                     + pbxsq;
    //   GRBQuadExpr qd2 = ctrl_pt_xyz[1][k] * ctrl_pt_xyz[1][k] - 2 * basepoint_(1) * ctrl_pt_xyz[1][k]
    //                     + pbysq;                        
    //   m_.addQConstr(qd1 + qd2 + ctrl_pt_xyz[2][k] * ctrl_pt_xyz[2][k] <= tethersq);
    // }        
  }
  // std::cout << green << "Gurobi finish adding all control point related constraints" << std::endl;

  // terminal v a constraint
  constrTerminalVA_ = new GRBConstr[6];
  for (int i=0; i<3; i++)
  {
    GRBLinExpr term_v, term_a;

    for (int k=0; k<3; k++) 
    {
      term_v += q_v_term(k) * poly_exp_[i][num_pol_init_-1][k];
    }
    constrTerminalVA_[i*2] = m_.addConstr(term_v ==0);


    for (int k=0; k<2; k++) // for terminal acc, only first 2 terms
    {
      term_a += q_a_term(k) * poly_exp_[i][num_pol_init_-1][k];
    }
    constrTerminalVA_[i*2+1] = m_.addConstr(term_a ==0);

  }

  // terminal position constraint
  GRBQuadExpr terminal_p;
  for (int i=0; i<3; i++)
  {
    for (int k1=0; k1<4; k1++)// for terminal pos
    {
      for (int k2=0; k2<4; k2++) 
      {
        terminal_p += Q_p_term(k1, k2) * poly_exp_[i][num_pol_init_-1][k1] * poly_exp_[i][num_pol_init_-1][k2];
      }

      terminal_p -= 2* q_p_term(k1) * poly_exp_[i][num_pol_init_-1][k1] * final_pos_(i);
    }
    terminal_p += final_pos_(i)*final_pos_(i);
  }  
  // std::cout << green << "Gurobi trying adding terminal constraints" << std::endl;

  Eigen::Vector3d init_pos(pwp_init_.coeff_x[0](3), pwp_init_.coeff_y[0](3), pwp_init_.coeff_z[0](3));
  //if initial and final are very close, then set termial constraint to be much smaller
  if ((init_pos - final_pos_).norm() < 1.0) 
  {
    m_.addQConstr(terminal_p - 0.10*0.10 <=0); //for exp the value is 0.15
  }
  // else
  // {
  //   m_.addQConstr(terminal_p - rad_term_*rad_term_ <=0);
  // }

  std::cout << green << "Gurobi finish adding constraints" << std::endl;

}

//case 1: no tresspassing between base and agent
//case 2: no tresspassing beyond agent
//case 3: no tresspassing beyond base
void PolySolverGurobi::addEntangleConstraintForIJCase(int agent_i, int pol_num, int ent_case,
                                                      std::vector<std::vector<GRBLinExpr>>& ctrl_pt_xyz)
{
  Eigen::Vector2d pointA, pointB;
  if (ent_case==1) 
  {
    pointA = (1-long_length_) * bendPtsForAgents_[agent_i].back() + 
             long_length_ * hullsNoInflation_[agent_i][pol_num].col(0);  
    pointB = hullsNoInflation_[agent_i][pol_num].col(0);      
  }
  else if (ent_case>1 && ent_case <= bendPtsForAgents_[agent_i].size())
  {
    //add a far point beyond agent i
    pointA = bendPtsForAgents_[agent_i][ent_case-2];
    pointB = bendPtsForAgents_[agent_i][ent_case-1];
  }
  else if (ent_case == bendPtsForAgents_[agent_i].size()+1)  
  {
    pointA = bendPtsForAgents_[agent_i].back();
    pointB = hullsNoInflation_[agent_i][pol_num].col(0);  
  }
  else return;

  double hulldist = 0.0;
  for (int k=0; k<3; k++) //add linear constraint
  {
    hulldist += (ctrlPtsInit_[pol_num].col(k+1) - ctrlPtsInit_[pol_num].col(k)).norm();
  }   
  if ((pointA - ctrlPtsInit_[pol_num].col(0)).norm() - hulldist > 0 &&
      (pointB - ctrlPtsInit_[pol_num].col(0)).norm() - hulldist > 0 )
    return;

  if (use_linear_constraints_)
  {
    Eigen::Vector3d sep_line;
    bool solved;
    solved = separator_solver_->solveModel(sep_line, pointA, pointB, ctrlPtsInit_[pol_num]);
    if (solved)
    {
      for (int k=0; k<4; k++) //add linear constraint
      {
        m_.addConstr(ctrl_pt_xyz[0][k] * sep_line(0) + ctrl_pt_xyz[1][k] * sep_line(1) +
                      sep_line(2) -1 <= 0);
      } 
    }
    else
    {
      std::cout << red << "not getting a initial separation line for entanglement constraint!" << std::endl;           
    }
  }
  else
  {
    GRBLinExpr var_a, var_b, var_c;
    var_a = GRBLinExpr(m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "ent_a_" + 
                  std::to_string(agent_i) +  std::to_string(pol_num) + std::to_string(ent_case)));
    var_b = GRBLinExpr(m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "ent_b_" + 
                  std::to_string(agent_i) +  std::to_string(pol_num) + std::to_string(ent_case)));
    var_c = GRBLinExpr(m_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "ent_c_" + 
                  std::to_string(agent_i) +  std::to_string(pol_num) + std::to_string(ent_case)));      
    for (int k=0; k<4; k++) //for our control points, add quadratic constraint
    {
      m_.addQConstr(ctrl_pt_xyz[0][k] * var_a + ctrl_pt_xyz[1][k] * var_b + var_c <= 0);
    }
    //add point A
    m_.addConstr(pointA(0) * var_a + pointA(1) * var_b + var_c >= 0);    
    //add point B                                     
    m_.addConstr(pointB(0) * var_a + pointB(1) * var_b + var_c >= 0);           
  }

}

void PolySolverGurobi::addEntangleConstraints()
{
  std::cout << green << "Gurobi adding entangle constraints" << std::endl;

  for (int i=0; i<num_of_agents_; i++)
  {
    if (i==id_-1) continue;
    if (vecOfAgentsBetas_[i].empty()) continue;
    for (int j=0; j<num_pol_init_; j++)
    {
      std::cout << green << "betas for agent "<<i<<" at "<<j<<"th seg is "
                <<vecOfAgentsBetas_[i][j]<< std::endl;
    }
  }

  std::cout << green << "Gurobi finish adding Entangle constraints" << std::endl;
}

bool PolySolverGurobi::optimize(double& objective_value)
{
  // reset some stuff
  MyTimer timer_function_setup(true);

  resetCompleteModel(m_);

  m_.set("OutputFlag", std::to_string(0));  // verbose (1) or not (0)
  m_.set("TimeLimit", std::to_string(max_runtime_));
  if (!use_linear_constraints_) m_.set("NonConvex", std::to_string(2));

  addObjective();
  addConstraints();
  m_.update();  // needed due to the lazy evaluation
  // std::cout<<bold<<white<<"setup time for gurobi is "<< timer_function_setup.ElapsedMs()<<" ms"<<std::endl;

  std::cout << "Starting optimization, allowing time = " << max_runtime_ * 1000 << " ms" << std::endl;
  MyTimer timer_function(true);

  m_.optimize();
  int optimstatus = m_.get(GRB_IntAttr_Status);
  printGurobiStatus(optimstatus);
  std::cout<<bold<<white<<"time optimization for gurobi is "<< timer_function.ElapsedMs()<<" ms"<<std::endl;
  total_replannings_++;

  int number_of_stored_solutions = m_.get(GRB_IntAttr_SolCount);
  bool found_soln = true;

  if ((optimstatus != GRB_OPTIMAL && optimstatus != GRB_TIME_LIMIT &&
       optimstatus != GRB_USER_OBJ_LIMIT &&                                    ///////////////
       optimstatus != GRB_ITERATION_LIMIT && optimstatus != GRB_NODE_LIMIT &&  ///////////////
       optimstatus != GRB_SOLUTION_LIMIT) ||
      number_of_stored_solutions <= 0) //optimization fails
  {
    for (int i=0; i<6; i++)
    {
      m_.remove(constrTerminalVA_[i]); //REMOVE terminal V A constraint
    }
    // m_.set("NumObj", std::to_string(0));
    m_.update();  
    addObjective(true);
    m_.update();    
    m_.optimize();
    optimstatus = m_.get(GRB_IntAttr_Status);    
    number_of_stored_solutions = m_.get(GRB_IntAttr_SolCount);
    printGurobiStatus(optimstatus);    
    if ((optimstatus != GRB_OPTIMAL && optimstatus != GRB_TIME_LIMIT &&
       optimstatus != GRB_USER_OBJ_LIMIT &&                                    ///////////////
       optimstatus != GRB_ITERATION_LIMIT && optimstatus != GRB_NODE_LIMIT &&  ///////////////
       optimstatus != GRB_SOLUTION_LIMIT) ||
      number_of_stored_solutions <= 0) //optimization fails again
    {
      std::cout << red << "Gurobi failed to find a solution, using initial guess" << reset
                << std::endl;
      pwp_out_ = pwp_init_;
      return false;
    }
  }
  std::cout << green << "Gurobi found a solution" << reset << std::endl;
  solutions_found_++;

  // copy the solution
  for (int i=0; i<num_pol_init_; i++)
  {
    pwp_out_.coeff_x.push_back(Eigen::Vector4d(poly_exp_[0][i][0].getValue(), poly_exp_[0][i][1].getValue(),
                                               poly_exp_[0][i][2].getValue(), poly_exp_[0][i][3].getValue()));

    pwp_out_.coeff_y.push_back(Eigen::Vector4d(poly_exp_[1][i][0].getValue(), poly_exp_[1][i][1].getValue(),
                                               poly_exp_[1][i][2].getValue(), poly_exp_[1][i][3].getValue()));

    pwp_out_.coeff_z.push_back(Eigen::Vector4d(poly_exp_[2][i][0].getValue(), poly_exp_[2][i][1].getValue(),
                                               poly_exp_[2][i][2].getValue(), poly_exp_[2][i][3].getValue()));      
  }

  //give more ascending speed if goal is too near
  Eigen::Vector2d init_pos(pwp_init_.coeff_x[0](3), pwp_init_.coeff_y[0](3));
  if ((init_pos - final_pos_.head(2)).norm() < 1.0) pwp_out_.coeff_z = pwp_init_.coeff_z;

  objective_value = m_.getObjective().getValue();

  std::cout << on_cyan << bold << "GUROBI Solved so far" << solutions_found_ << "/" << total_replannings_ << reset
            << std::endl;
  return true;
}

void PolySolverGurobi::generatePwpOut(mt::PieceWisePol& pwp_out, 
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