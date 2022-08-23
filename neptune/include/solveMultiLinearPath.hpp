#include <Eigen/Dense>
#include <stdio.h>
#include "cgal_utils.hpp"
#include "backward.hpp"
#include <unordered_map>
#include <queue>
#define inf 1<<20

using namespace Eigen;

class mtlp 
{
public:
	mtlp(double _resolution, double x_min, double x_max, double y_min, 
	                double y_max, double z_min, double z_max, double radius)
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
	        Vector3i tmpIdx(i,i,i);
	        Vector3d pos(0.0, 0.0, 0.0);
	        GridNodeMap[i] = new GridNode(tmpIdx, pos);
	    }
	    radius_squared = radius * radius;
	}	
	bool solveMultiTetherLinearPath(std::vector<Vector3d>& basepoints, std::vector<Matrix<double, 3, 2>>& startEndPts,
	std::vector<std::vector<Vector3d>>& outputPath);

private:
	struct GridNode;

	typedef GridNode* GridNodePtr;

	struct GridNode
	{     
	    int id;        // 1--> open set, -1 --> closed set
	    Eigen::Vector3d coord; 
	    Eigen::Vector3i dir;   // direction of expanding
	    Eigen::Vector3i index;
		
	    double gScore, fScore;
	    GridNodePtr cameFrom;
	    std::multimap<double, GridNodePtr>::iterator nodeMapIt;
		double lengthUsed;

	    GridNode(Eigen::Vector3i _index, Eigen::Vector3d _coord){  
			id = 0;
			index = _index;
			coord = _coord;
			dir   = Eigen::Vector3i::Zero();

			gScore = inf;
			fScore = inf;
			cameFrom = NULL;
	    }

	    GridNode(){};
	    ~GridNode(){};
	};	

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

	  void insert(Eigen::Vector3i idx, GridNodePtr node) {
	    data_3d_.insert(std::make_pair(idx, node));
	  }

	  GridNodePtr find(Eigen::Vector2i idx) {
	    auto iter = data_2d_.find(idx);
	    return iter == data_2d_.end() ? NULL : iter->second;
	  }

	  GridNodePtr find(Eigen::Vector3i idx) {
	    auto iter = data_3d_.find(idx);
	    return iter == data_3d_.end() ? NULL : iter->second;
	  }

	  void clear() {
	    data_3d_.clear();
	    data_2d_.clear();
	  }
	};

	std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, CompareCost> openList_;  //= OpenSet, = Q
	GridNodePtr * GridNodeMap;
	int GLX_SIZE, GLY_SIZE, GLZ_SIZE;
	int GLXYZ_SIZE, GLYZ_SIZE;	
	uint8_t * data;
	Eigen::Vector3i goalIdx;

	double resolution, inv_resolution;
	double gl_xl, gl_yl, gl_zl;
	double gl_xu, gl_yu, gl_zu;
 	GridNodePtr endPtr; 

	GridNodePtr terminatePtr;
	std::multimap<double, GridNodePtr, std::less<double>> openSet;
	GridNodeHashTable expanded_nodes_;
	int node_used_num_ = 0;
	double bias_ = 1.1;
	bool use_optimal_path_finding_ = true;
	std::vector<Matrix<double, 3, 2>> cableLines_;
	Vector3d basepoint_;
	double radius_squared = 0.0;
	double goal_height_ = 5.0;

	void AstarGraphSearch(Vector3d start_pt, Vector3d end_pt, 
                         std::vector<Matrix<double, 3, 2>>& cable_lines, 
                         Vector3d basepoint, int& search_result);
	void resetUsedGrids();
	void resetGrid(GridNodePtr ptr);

	Vector3i coord2gridIndex(const Vector3d & pt) 
	{
	    Vector3i idx;
	    idx <<  std::min( std::max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
	            std::min( std::max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
	            std::min( std::max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                        
	  
	    return idx;
	}

	Vector3d gridIndex2coord(const Vector3i & index) 
	{
	    Vector3d pt;
	    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
	    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
	    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

	    return pt;
	}

	double getHeu(GridNodePtr node1, GridNodePtr node2)
	{

	    //int Heu_type = 3; // 0-dijkstra, 1-manhattan, 2-euclidean, 3-diagonal
	    Vector3d diff = node1->index.cast<double>() -node2->index.cast<double>();
	    // std::cout<<"diff is "<<diff<<std::endl;
		double	heu_cost = diff.lpNorm<2>();
	    // ROS_WARN("heuristic cost is %f", heu_cost);
	    return heu_cost;
	}
	void AstarGetSucc(GridNodePtr currentPtr);
	std::vector<Vector3d> getPath();
	bool intersectSphereSegment(Sphere_3 sphe, Matrix<double, 3, 2> segment);


};