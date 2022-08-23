/* ----------------------------------------------------------------------------
 * Copyright 2021, Cao Muqing 
 * Nanyang Technological University
 * All Rights Reserved
 * Authors: Cao Muqing, et al.
 * -------------------------------------------------------------------------- */

#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <vector>

namespace eu  // cgal utils
{
enum
{
	NOT_ENTANGLED=0,
	RISK_ENTANGLING=1,
	ENTANGLED=2
};

struct ent_state
{
  std::vector<Eigen::Vector2i> alphas; //(agent_id, case_no)
  std::vector<double> betas; //stores the w for each ent_case
  std::vector<int> bendPointsIdx; //list of points where the tether bend
  std::vector<int> active_cases;
};

struct ent_state_orig
{
  std::vector<Eigen::Vector2i> hsig; //(agent_id, case_no)
  // std::vector<int> bendPointsIdx; //list of points where the tether bend
  std::vector<Eigen::Vector2d> contPts;
};

double vectorWedge2(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c);
double vectorWedge2(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c, 
                        Eigen::Vector2d& ab, Eigen::Vector2d& ac);

bool lineSegIntersect(Eigen::Vector2d& a1, Eigen::Vector2d& a2, Eigen::Vector2d& b1, Eigen::Vector2d& b2);

int checkMaxPointAnchor(double ll, double cableLength, double maxheight);

// bool checkEntangleUpdate(int& state, Eigen::Vector2d& pk, Eigen::Vector2d& pkplus1, Eigen::Vector2d& pik, Eigen::Vector2d& pikplus1);
// bool checkEntangleUpdate(int& state, Eigen::Vector2d& pk, Eigen::Vector2d& pkplus1, Eigen::Vector2d& pik, Eigen::Vector2d& pikplus1, 
// 						Eigen::Vector2d& pb, Eigen::Vector2d& pbi, double cableLength, double maxheight=2.0);

bool checkEntangleUpdate(int& state, int& state2, Eigen::Vector2d& pk, Eigen::Vector2d& pkplus1, Eigen::Vector2d& pik,
                               Eigen::Vector2d& pikplus1, Eigen::Vector2d& pb, Eigen::Vector2d& pbi, double cableLength,
                               double maxheight);
bool checkEntangleUpdate2(int& state, int& state2, Eigen::Vector2d& pk, Eigen::Vector2d& pkplus1, Eigen::Vector2d& pik,
                               Eigen::Vector2d& pikplus1, Eigen::Vector2d& pb, Eigen::Vector2d& pbi, double cableLength,
                               double maxheight);
int checkEntangleUpdate3(int& state, int& state2, int& state3, const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1, 
                             const Eigen::Vector2d& pik, const Eigen::Vector2d& pikplus1, const Eigen::Vector2d& pb, 
                             const Eigen::Vector2d& pbi, double cableLength, double maxheight);

int entangleStatesUpdate3(int& state, int& state2, int& state3, const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1, 
                             const Eigen::Vector2d& pik, const Eigen::Vector2d& pikplus1, const Eigen::Vector2d& pb, 
                             const Eigen::Vector2d& pbi, double& d1, double& e1);

int checkEntangleStaticUpdate(int& state, int& state2, int& state3, const Eigen::Vector2d& pk, 
                              const Eigen::Vector2d& pkplus1, const Eigen::Vector2d& pik, 
                              const Eigen::Vector2d& pb, const Eigen::Vector2d& pbi, double cableLength, double maxheight);

void entangleHSigUpdate3(int id, std::vector<Eigen::Vector2i>& alphas, std::vector<Eigen::Vector2d>& betas, 
                        int& active_case_num, const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1, 
                        const Eigen::Vector2d& pik, const Eigen::Vector2d& pikplus1, const Eigen::Vector2d& pb, 
                        const Eigen::Vector2d& pbi);

void addOrRemoveLastElement(int agent_id, int case_id, std::vector<Eigen::Vector2i>& alphas, 
                            std::vector<Eigen::Vector2d>& betas, int& active_case_num, double d1, double e1);

void shrinkAlphaBetas(std::vector<Eigen::Vector2i>& alphas, std::vector<Eigen::Vector2d>& betas, 
                               std::vector<int>& active_cases, int num_of_agents);

void entangleHSigToAddStatic(std::vector<Eigen::Vector2i>& alphasToAdd,  
                            const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1,  
                            std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, int num_of_agents);

void addAlphaBetaToList(std::vector<Eigen::Vector2i>& alphasToAdd, ent_state& entangle_state, 
                        const Eigen::Vector2d& pk, std::vector<Eigen::Vector2d>& pb, Eigen::Vector2d& pb_self,
                        std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, int num_of_agents,
                        std::vector<std::vector<Eigen::Vector2d>>& bendpts);

void entangleHSigToAddAgentInd(std::vector<Eigen::Vector2i>& alphasToAdd,  
                            const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1,  
                            const Eigen::Vector2d& pik, Eigen::Vector2d& pikplus1,
                            const Eigen::Vector2d& pb, std::vector<Eigen::Vector2d>& bendpts, int agent_id);

void entangleHSigToAddAgentInd(std::vector<Eigen::Vector2i>& alphasToAdd, const Eigen::Vector2d& pk, 
                                   const Eigen::Vector2d& pkplus1, const Eigen::Vector2d& pik, 
                                   Eigen::Vector2d& pikplus1, const Eigen::Vector2d& pb, 
                                   std::vector<Eigen::Vector2d>& bendpts, std::vector<Eigen::Vector2d>& bendpts_prev, 
                                   int agent_id);

bool breakcondition(Eigen::Vector2i& alphaToAdd, Eigen::Vector2i& alphaInList, int num_of_agents,
                        int idx_to_check, int idx_last_bend,
                        std::vector<std::vector<Eigen::Vector2d>>& bendpts);

void getBendPt2d(Eigen::Vector2d& bp, eu::ent_state& entangle_state, std::vector<Eigen::Vector2d>& pb, 
                     Eigen::Vector2d& pb_self, std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, 
                     int num_of_agents);

void getBendPt2dwIdx(Eigen::Vector2d& bp, eu::ent_state& entangle_state, std::vector<Eigen::Vector2d>& pb, 
                     std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, 
                     std::vector<Eigen::Vector2d>& staticObsLongestDist,
                     int num_of_agents, int idxFromBeg, double& compensation);

double calculateBetaForCase(const Eigen::Vector2i& alphaToAddOrUpdate, const Eigen::Vector2d& pk, 
                            std::vector<Eigen::Vector2d>& pb, Eigen::Vector2d& bp, 
                            std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, int num_of_agents);

void updateBendPts(ent_state& entangle_state, const Eigen::Vector2d& pkplus1, 
                   std::vector<Eigen::Vector2d>& pb, Eigen::Vector2d& pb_self,
                   std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, int num_of_agents);

double getTetherLength(eu::ent_state& entangle_state, std::vector<Eigen::Vector2d>& pb, 
                         Eigen::Vector2d pb_self, Eigen::Vector2d& pkplus1, 
                         std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep, 
                         std::vector<Eigen::Vector2d>& staticObsLongestDist,
                         int num_of_agents);

void updateHSig_orig(ent_state_orig& ent_state,  
                    const Eigen::Vector2d& pk, const Eigen::Vector2d& pkplus1,  
                    std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep);

void updateContPts_orig(std::vector<Eigen::Vector2d>& contPts,  
                        std::vector<Eigen::Matrix<double, 2, -1>>& convexHullOfStaticObs);

double tetherLength_orig(std::vector<Eigen::Vector2d>& contPts);

}  