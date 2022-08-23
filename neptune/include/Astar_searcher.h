#ifndef _ASTART_SEARCHER_H
#define _ASTART_SEARCHER_H

#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include "backward.hpp"
#include <unordered_map>
#include <queue>
#include "mader_types.hpp"
#include "entangle_utils.hpp"

#define inf 1<<20
struct GridNode;
typedef GridNode* GridNodePtr;

struct GridNode
{     
    int id;        // 1--> open set, -1 --> closed set
    Eigen::Vector2d coord; 
    Eigen::Vector2i dir;   // direction of expanding
    Eigen::Vector2i index;
	
    double gScore, fScore;
    GridNodePtr cameFrom;
    std::multimap<double, GridNodePtr>::iterator nodeMapIt;
	eu::ent_state_orig entangle_state;
	double lengthUsed;

    GridNode(Eigen::Vector2i _index, Eigen::Vector2d _coord){  
		id = 0;
		index = _index;
		coord = _coord;
		dir   = Eigen::Vector2i::Zero();

		gScore = inf;
		fScore = inf;
		cameFrom = NULL;
    }

    GridNode(){};
    ~GridNode(){};
};

// Taken from https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
template <typename T>
struct matrix_hash1_ : std::unary_function<T, size_t>
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

class GridNodeHashTable {
 private:
  /* data */
  std::unordered_map<Eigen::Vector3i, GridNodePtr, matrix_hash1_<Eigen::Vector3i>>
      data_3d_;
  std::unordered_map<Eigen::Vector2i, GridNodePtr, matrix_hash1_<Eigen::Vector2i>>
      data_2d_;

 public:
  GridNodeHashTable(/* args */) {}
  ~GridNodeHashTable() {}
  void insert(Eigen::Vector2i idx, GridNodePtr node) {
    data_2d_.insert(std::make_pair(idx, node));
  }

  void insert(Eigen::Vector2i idx, int time_idx, GridNodePtr node) {
    data_3d_.insert(std::make_pair(
        Eigen::Vector3i(idx(0), idx(1), time_idx), node));
  }

  GridNodePtr find(Eigen::Vector2i idx) {
    auto iter = data_2d_.find(idx);
    return iter == data_2d_.end() ? NULL : iter->second;
  }
 
  GridNodePtr find(Eigen::Vector2i idx, int time_idx) {
    auto iter =
        data_3d_.find(Eigen::Vector3i(idx(0), idx(1), time_idx));
    return iter == data_3d_.end() ? NULL : iter->second;
  }

  void clear() {
    data_3d_.clear();
    data_2d_.clear();
  }
};

class AstarPathFinder
{	
	private:
		 struct CompareCost
		 {
		   double bias;
		   bool operator()(const GridNodePtr left, const GridNodePtr right)
		   {
		     if (fabs(left->fScore - right->fScore) < 1e-5)
		     {
		       return left->fScore-left->gScore > right->fScore-right->gScore;  // If two costs are ~the same, decide only upon heuristic
		     }
		     else
		     {
		       return left->fScore > right->fScore;
		     }
		   }
		 };
	    int num_of_static_obst_ = 0;
	    std::vector<mt::Polygon_Std> inflatedHullOfStaticObs_;
  		std::vector<mt::Polygon_Std> noInflatedHullOfStaticObs_;

	    eu::ent_state_orig entangle_state_init_;
		std::vector<Eigen::Matrix<double, 2, 2>> staticObsRep_;
		ros::Duration time_homotopy_update;
		ros::Duration time_length_check;
		double time_astar;		
		double time_length_check_ms;
  		std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, CompareCost> openList_;  //= OpenSet, = Q
  		double bias_ = 1.1;
		eu::ent_state_orig terminal_ent_state_;

	protected:
		uint8_t * data;
		GridNodePtr * GridNodeMap;
		Eigen::Vector2i goalIdx;
		int GLX_SIZE, GLY_SIZE, GLZ_SIZE;
		int GLXYZ_SIZE, GLYZ_SIZE;	
		double cableLength_ = 1000;

		double resolution, inv_resolution;
		double gl_xl, gl_yl, gl_zl;
		double gl_xu, gl_yu, gl_zu;
    	GridNodePtr endPtr; 

		GridNodePtr terminatePtr;
		std::multimap<double, GridNodePtr, std::less<double>> openSet;
  		GridNodeHashTable expanded_nodes_;
  		int node_used_num_ = 0;

		double getHeu(GridNodePtr node1, GridNodePtr node2, int Heu_type=2, bool tie_break=true);
		void AstarGetSucc(GridNodePtr currentPtr);		

		bool isFree(const int & idx_x, const int & idx_y, const int & idx_z) const;
		bool isFree(const Eigen::Vector3i & index) const;
		bool isOccupied(const Eigen::Vector2i index);

		Eigen::Vector3d gridIndex2coord(const Eigen::Vector3i & index);
		Eigen::Vector2d gridIndex2coord(const Eigen::Vector2i & index);

		Eigen::Vector2i coord2gridIndex(const Eigen::Vector2d & pt);
		double tetherLength(std::vector<Eigen::Vector2d>& contPts);
		unsigned int getIz(eu::ent_state_orig& entangle_state);
		unsigned int power_int(unsigned int base, unsigned int exponent);

	public:
		AstarPathFinder(){};
		~AstarPathFinder(){};
		void AstarGraphSearch(Eigen::Vector2d start_pt, Eigen::Vector2d end_pt, int& search_result);
		void resetGrid(GridNodePtr ptr);
		void resetUsedGrids();

		void init(double _resolution, double x_min, double x_max, double y_min, 
                  double y_max, double z_min, double z_max, double cableLength);
    	void setObs(const double coord_x, const double coord_y, const double coord_z);

		Eigen::Vector2d coordRounding(const Eigen::Vector2d & coord);
		std::vector<Eigen::Vector3d> getPath(double& length);
		std::vector<Eigen::Vector2d> getVisitedNodes();
		void setStaticObstVert(std::vector<mt::Polygon_Std>& inflatedHullOfStaticObs, 
							   std::vector<mt::Polygon_Std>& noInflatedHullOfStaticObs);
		void setHsig(eu::ent_state_orig& entangle_state_orig_);
		void setStaticObstRep(std::vector<Eigen::Matrix<double, 2, 2>>& staticObsRep);
		void getRuntime(double& runtime_this_round, double& time_spent_contact_pt, int& node_used_num);
		void getHsigTerminal(eu::ent_state_orig& entangle_state_orig_);


};

#endif