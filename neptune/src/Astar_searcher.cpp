/* ----------------------------------------------------------------------------
 * A implementation of the algorithm presented in
 * Path Planning for a Tethered Mobile Robot, ICRA 2014
 * Copyright 2020, Cao Muqing
 * Nanyang Technology University
 * All Rights Reserved
 * Author: Cao Muqing
 * See LICENSE file for the license information
 * 
 * -------------------------------------------------------------------------- */

#include "Astar_searcher.h"
# include "gjk.hpp"

using namespace std;
using namespace Eigen;

void AstarPathFinder::init(double _resolution, double x_min, double x_max, double y_min, 
                            double y_max, double z_min, double z_max, double cableLength)
{   
    gl_xl = x_min;
    gl_yl = y_min;
    gl_zl = z_min;

    gl_xu = x_max;
    gl_yu = y_max;
    gl_zu = z_max;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    double max_x_id = (int)((x_max-x_min) * inv_resolution);
    double max_y_id = (int)((y_max-y_min) * inv_resolution);
    double max_z_id = (int)((z_max-z_min) * inv_resolution);

    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
    
    GridNodeMap = new GridNodePtr[GLXYZ_SIZE];
    for(int i = 0; i < GLXYZ_SIZE; i++){
        Vector2i tmpIdx(i,i);
        Vector2d pos(0.0, 0.0);
        GridNodeMap[i] = new GridNode(tmpIdx, pos);
    }
    cableLength_ = cableLength;
}

void AstarPathFinder::resetGrid(GridNodePtr ptr)
{
    ptr->id = 0;
    ptr->cameFrom = NULL;
    ptr->gScore = inf;
    ptr->fScore = inf;
}

void AstarPathFinder::resetUsedGrids()
{   
    for(int i=0; i < GLXYZ_SIZE; i++)
        resetGrid(GridNodeMap[i]);
}

void AstarPathFinder::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      

    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

vector<Vector2d> AstarPathFinder::getVisitedNodes()
{   
    vector<Vector2d> visited_nodes;

    for (int i=0; i<node_used_num_; i++)
    {
        visited_nodes.push_back(GridNodeMap[i]->coord);
    }

    // ROS_WARN("visited_nodes size : %d", (int)visited_nodes.size());
    return visited_nodes;
}

Vector3d AstarPathFinder::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector2d AstarPathFinder::gridIndex2coord(const Vector2i & index) 
{
    Vector2d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;

    return pt;
}

Vector2i AstarPathFinder::coord2gridIndex(const Vector2d & pt) 
{
    Vector2i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector2d AstarPathFinder::coordRounding(const Eigen::Vector2d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

inline bool AstarPathFinder::isFree(const Eigen::Vector3i & index) const
{
    return isFree(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isFree(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

inline bool AstarPathFinder::isOccupied(const Eigen::Vector2i index)
{
    if (index(0) < 0 || index(0) >= GLX_SIZE || index(1) < 0 || index(1) >= GLY_SIZE)
        return true;
    Vector2d coord_center = gridIndex2coord(index);
    Eigen::Matrix<double, 2, 4> Cps;
    double res = 0.05;
    Cps << coord_center(0) + res, coord_center(0) + res,
           coord_center(0) - res, coord_center(0) - res,
           coord_center(1) + res, coord_center(1) - res,
           coord_center(1) + res, coord_center(1) - res;

      for (int obst_index = 0; obst_index < num_of_static_obst_; obst_index++)
      {

        bool gjk_collision = gjk::collision(inflatedHullOfStaticObs_[obst_index], Cps);

        if (gjk_collision == true)
        {
            return true;
        }    
      }      
    return false;
}

double AstarPathFinder::tetherLength(std::vector<Eigen::Vector2d>& contPts)
{
    double length = 0.0;
    for (int i=0; i<contPts.size()-1; i++){
        length += (contPts[i+1]-contPts[i]).norm();
    }
    return length;
}

unsigned int AstarPathFinder::getIz(eu::ent_state_orig& entangle_state)
{ 
    int iz;
    for (int i=0; i<entangle_state.hsig.size(); i++)
    {
      iz += (i+1)*power_int(2, entangle_state.hsig[i](0)) + entangle_state.hsig[i](1);
    }
    return iz;
}

unsigned int AstarPathFinder::power_int(unsigned int base, unsigned int exponent)
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

inline void AstarPathFinder::AstarGetSucc(GridNodePtr currentPtr)
{   

    int m = 0;
    for (int i=-1; i<=1; i++)
        for (int j=-1; j<=1; j++){
            if (i==0 && j==0) continue; //skip the current point
            Vector2i offsetIndex(i, j);
            Vector2i neighbourIndex = currentPtr->index + offsetIndex;
            double offsetCost = offsetIndex.cast<float>().lpNorm<2>();
            if (isOccupied(neighbourIndex)) continue;

            //update entangle state
            ros::Time time1 = ros::Time::now();
            eu::ent_state_orig ent_state_next = currentPtr->entangle_state;
            Vector2d pkplus1 = gridIndex2coord(neighbourIndex);
            eu::updateHSig_orig(ent_state_next, currentPtr->coord, pkplus1, staticObsRep_);            
            ent_state_next.contPts.push_back(pkplus1);
            ros::Time time2 = ros::Time::now();
            time_homotopy_update += (time2 - time1);

            //check cable length first
            ros::Time time11 = ros::Time::now();
            double _lengthUsed = currentPtr->lengthUsed + offsetCost;
            if (_lengthUsed > cableLength_)
            {
                eu::updateContPts_orig(ent_state_next.contPts, noInflatedHullOfStaticObs_);
                _lengthUsed = eu::tetherLength_orig(ent_state_next.contPts);
                if (_lengthUsed > cableLength_) continue;
            }                
            ros::Time time22 = ros::Time::now();
            time_length_check += (time22 - time11);

            int homotopy_index = getIz(ent_state_next);
            // Vector3i locateIndex(neighbourIndex(0), neighbourIndex(1), homotopy_index);
            GridNodePtr nodeptr;
            nodeptr = expanded_nodes_.find(neighbourIndex, homotopy_index);        
            if (nodeptr == NULL) { //not expanded at all

                if (node_used_num_ >= GLXYZ_SIZE) return;
                nodeptr = GridNodeMap[node_used_num_++];
                nodeptr-> id = 1;
                nodeptr-> index = neighbourIndex;                
                nodeptr-> gScore = offsetCost + currentPtr->gScore; //get accumulated g cost
                nodeptr-> fScore = nodeptr->gScore + bias_ * getHeu(nodeptr, endPtr);
                nodeptr-> cameFrom = currentPtr;
                nodeptr-> coord = pkplus1;
                nodeptr-> entangle_state = ent_state_next;
                nodeptr-> lengthUsed = _lengthUsed;
                // nodeptr->nodeMapIt = openSet.insert(make_pair(nodeptr->fScore, nodeptr)); //push the neighbour node into the queue
                openList_.push(nodeptr);
                expanded_nodes_.insert(neighbourIndex, homotopy_index, nodeptr);

                continue;                

            } else if (nodeptr->id==1) { //in open list already
                double newgScore = offsetCost + currentPtr->gScore; //get accumulated g cost
                if (nodeptr->gScore > newgScore){       //compare new g cost with previous g cost, if lower then update the openset
                    nodeptr->gScore = newgScore;
                    nodeptr->fScore = newgScore + bias_ * getHeu(nodeptr, endPtr);
                    nodeptr->cameFrom = currentPtr;
                    nodeptr-> entangle_state = ent_state_next;
                    nodeptr-> lengthUsed = _lengthUsed;

                    // openSet.erase(nodeptr->nodeMapIt);
                    // nodeptr->nodeMapIt = openSet.insert(make_pair(nodeptr->fScore, nodeptr)); //push the neighbour node into the queue
                }                

            } else { // in closed list already, dun care
                // double newgScore = offsetCost + currentPtr->gScore; //get accumulated g cost
                // if (nodeptr->gScore > newgScore){       //compare new g cost with previous g cost, if lower then update the openset
                //     nodeptr->gScore = newgScore;
                //     nodeptr->fScore = newgScore + getHeu(nodeptr, endPtr);
                //     nodeptr->cameFrom = currentPtr;

                //     nodeptr-> id = 1;
                //     nodeptr->nodeMapIt = openSet.insert(make_pair(nodeptr->fScore, nodeptr)); //push the neighbour node into the queue
                // }  
            }
    
            m++;

        }
}

double AstarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2, int Heu_type, bool tie_break)
{

    //int Heu_type = 3; // 0-dijkstra, 1-manhattan, 2-euclidean, 3-diagonal
    Vector2d diff = node1->index.cast<double>() -node2->index.cast<double>();
    // std::cout<<"node1 is "<<node1->index<<std::endl;
    // std::cout<<"node2 is "<<node2->index<<std::endl;

    // std::cout<<"diff is "<<diff<<std::endl;
    double heu_cost = 0;
    if (!Heu_type){
    } else if (Heu_type == 1){     //manhattan distance
        // Vector3d diff = node1->coord - node2->coord;
        // double distance = diff.cwiseAbs().sum();
        heu_cost = diff.lpNorm<1>();
    } else if (Heu_type == 2){
        heu_cost = diff.lpNorm<2>();
    } else if (Heu_type == 3){
        diff = diff.cwiseAbs();
        heu_cost = diff.sum()-1.26794919*diff.minCoeff();
    } else {
        ROS_WARN("Unrecognized heuristic type, using 0 value for H function instead");
    }
    if (tie_break) heu_cost*=1.01;
    // ROS_WARN("heuristic cost is %f", heu_cost);
    return heu_cost;
}

void AstarPathFinder::AstarGraphSearch(Vector2d start_pt, Vector2d end_pt, int& search_result)
{   
    CompareCost comparer;
    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, CompareCost> openList_tmp(comparer);
    openList_.swap(openList_tmp);      

    resetUsedGrids();
    time_homotopy_update = ros::Duration(0);
    time_length_check = ros::Duration(0);    
    time_astar = 0.0; 
    ros::Time time_1 = ros::Time::now();    

    //index of start_point and end_point
    Vector2i start_idx = coord2gridIndex(start_pt);
    Vector2i end_idx   = coord2gridIndex(end_pt);
    goalIdx = end_idx;

    //position of start_point and end_point
    // start_pt = gridIndex2coord(start_idx);
    end_pt   = gridIndex2coord(end_idx);

    //Initialize the pointers of struct GridNode which represent start node and goal node
    node_used_num_ = 0;
    GridNodePtr startPtr = GridNodeMap[node_used_num_++];
    startPtr ->index = start_idx;
    endPtr   = new GridNode(end_idx,   end_pt);

    //openSet is the open_list implemented through multimap in STL library
    openSet.clear();
    expanded_nodes_.clear();

    // currentPtr represents the node with lowest f(n) in the open_list
    GridNodePtr currentPtr  = NULL;
    GridNodePtr neighborPtr = NULL;

    //put start node in open set
    startPtr -> gScore = 0;
    startPtr -> fScore = bias_ * getHeu(startPtr, endPtr);   
    //STEP 1: finish the AstarPathFinder::getHeu , which is the heuristic function
    startPtr -> id = 1; 
    startPtr -> coord = start_pt;
    startPtr -> entangle_state = entangle_state_init_;
    startPtr -> lengthUsed = eu::tetherLength_orig(entangle_state_init_.contPts);
    startPtr -> cameFrom = NULL; 

    openList_.push(startPtr);
    // openSet.insert( make_pair(startPtr -> fScore, startPtr) );
    int homotopy_index = 0;
    expanded_nodes_.insert(start_idx, homotopy_index, startPtr);    

    vector<GridNodePtr> neighborPtrSets;
    vector<double> edgeCostSets;

    // this is the main loop
    // while ( !openSet.empty() ){
    while ( !openList_.empty() ){
        if ((ros::Time::now() - time_1).toSec() > 30)
        {
            search_result = 0;
            return;
        }

        // multimap<double, GridNodePtr> :: iterator it = openSet.begin(); 
        // currentPtr = (*it).second;
        currentPtr = openList_.top();
        openList_.pop();
        // ROS_WARN("Currently expand a node with coord %f, %f", currentPtr->coord(0), currentPtr->coord(1));        
        // ROS_WARN("Currently expand a node with Fscore %f", currentPtr->fScore);
        // openSet.erase(it);
        currentPtr->id = -1; //set it to be in closed set

        // if the current node is the goal 
        if( currentPtr->index == goalIdx ){
            search_result = 1;
            ros::Time time_2 = ros::Time::now();
            terminatePtr = currentPtr;
            terminal_ent_state_ = terminatePtr->entangle_state;
            time_astar = (time_2 - time_1).toSec() * 1000.0;
            time_length_check_ms = time_length_check.toSec() * 1000.0;
            ROS_WARN("[A*]{sucess}  Time in A*  is %f ms, path cost if %f m", time_astar, currentPtr->gScore * resolution );
            ROS_WARN("[A*]{sucess}  Time in homotopy update is %f ms", time_homotopy_update.toSec() * 1000.0);
            ROS_WARN("[A*]{sucess}  Time in Length check is %f ms", time_length_check_ms);
            ROS_WARN("[A*]{sucess}  Number of expansion is %d", node_used_num_ );

            return;
        }
        //get the succetion
        AstarGetSucc(currentPtr);     
      
    }
    //if search fails
    search_result = -1;
    ros::Time time_2 = ros::Time::now();
    if((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("[A*]{Fail} Time consume in Astar path finding is %f", (time_2 - time_1).toSec() );
}


vector<Vector3d> AstarPathFinder::getPath(double& length) 
{   
    vector<Vector3d> path;
    vector<GridNodePtr> gridPath;
    length = 0;

    GridNodePtr gptr=terminatePtr;
    while(gptr->cameFrom!=NULL){
        gridPath.push_back(gptr);
        gptr=gptr->cameFrom;

    }
    gridPath.push_back(gptr);

    GridNodePtr lastnode = gridPath.back();
    for (auto ptr: gridPath)
    {
        path.push_back(Vector3d(ptr->coord(0), ptr->coord(1),0));
    }
        
    reverse(path.begin(),path.end());
    for (int i = 0; i<path.size()-1; i++)
    {
        length += (path[i+1] - path[i]).norm();
    }
    return path;
}

void AstarPathFinder::setStaticObstVert(std::vector<mt::Polygon_Std>& inflatedHullOfStaticObs, 
                       std::vector<mt::Polygon_Std>& noInflatedHullOfStaticObs)
{
  num_of_static_obst_ = inflatedHullOfStaticObs.size();
  inflatedHullOfStaticObs_ = inflatedHullOfStaticObs;
  noInflatedHullOfStaticObs_ = noInflatedHullOfStaticObs;
}

void AstarPathFinder::setHsig(eu::ent_state_orig& entangle_state_orig_)
{
    entangle_state_init_ = entangle_state_orig_;
}

void AstarPathFinder::getHsigTerminal(eu::ent_state_orig& entangle_state_orig_)
{
    entangle_state_orig_ = terminal_ent_state_;
}

void AstarPathFinder::setStaticObstRep(std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep)
{
  staticObsRep_ = staticObsRep;
}

void AstarPathFinder::getRuntime(double& runtime_this_round, double& time_spent_contact_pt, 
                                int& node_used_num)
{
  runtime_this_round = time_astar;
  time_spent_contact_pt = time_length_check_ms;
  node_used_num = node_used_num_;
}