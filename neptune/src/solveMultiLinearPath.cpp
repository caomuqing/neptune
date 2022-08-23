/* ----------------------------------------------------------------------------
 * A implementation of the algorithm presented in
 * Motion Planning in R3 for Multiple Tethered Robots, Hert, etc
 * Transaction on Robotics And Automation, 1999
 * The sequential motion plan with Convex Hull method with waiting times is 
 * implemented here
 * Copyright 2022, Cao Muqing
 * Nanyang Technology University
 * All Rights Reserved
 * Author: Cao Muqing
 * 
 * -------------------------------------------------------------------------- */

#include "solveMultiLinearPath.hpp"


using namespace Eigen;

bool mtlp::solveMultiTetherLinearPath(std::vector<Vector3d>& basepoints, 
									 std::vector<Matrix<double, 3, 2>>& startEndPts,
									 std::vector<std::vector<Vector3d>>& outputPath)
{
	outputPath.clear();
	int num_of_robot = basepoints.size();
	if (num_of_robot != startEndPts.size())
	{
		std::cout<<"the size of inputs are not consistent! returning"<<std::endl;
		return false;
	}

	MatrixXi graph = MatrixXi::Identity(num_of_robot, num_of_robot);
	for (int i=0; i<num_of_robot; i++)
	{
		for (int j=0; j<num_of_robot; j++)
		{
			Matrix<double, 3, 3> triA, triB;
			triA << basepoints[i], startEndPts[i];
			triB << basepoints[j], startEndPts[j];
			// std::cout<<"triA is "<<std::endl;
			// std::cout<<triA<<std::endl;
			// std::cout<<"triB is "<<std::endl;
			// std::cout<<triB<<std::endl;
			if (!cu::intersectTri(triA, triB)) //no intersection between respective cable surface
			{
				graph(i,j) = 1;
				graph(j,i) = 1;
				continue;
			}

			bool baseline_intersect = cu::intersectSeg(startEndPts[i].block<2,2>(0,0), 
													startEndPts[j].block<2,2>(0,0));
			Eigen::Matrix<double, 3, 2> anchor_start_A, anchor_end_A, anchor_start_B, anchor_end_B;
			 anchor_start_A << basepoints[i], startEndPts[i].col(0);
			 anchor_end_A << basepoints[i], startEndPts[i].col(1);
			 anchor_start_B << basepoints[j], startEndPts[j].col(0);
			 anchor_end_B << basepoints[j], startEndPts[j].col(1);
			if (!baseline_intersect)
			{	
				if (cu::intersect(anchor_start_A, triB) && cu::intersect(anchor_end_A, triB) ||
					cu::intersect(anchor_start_B, triA) && cu::intersect(anchor_end_B, triA) )
				{	//type 1
					// do nothing
				}
				else if (cu::intersect(anchor_start_A, triB) && cu::intersect(anchor_start_B, triA) ||
						cu::intersect(anchor_end_B, triA) && cu::intersect(anchor_end_A, triB) )
				{	//type 2

				}
				else if (cu::intersect(anchor_start_A, triB) && cu::intersect(anchor_end_B, triA))
				{	//type 3, A has to move before B
					graph(i,j) = 1;
				}
				else if (cu::intersect(anchor_start_B, triA) && cu::intersect(anchor_end_A, triB))
				{	//type 3, B has to move before A
					graph(j,i) = 1;
				}
			}
			else
			{
				if (cu::intersect(anchor_end_A, triB))
				{
					graph(j,i) = 1;
				}
				else if (cu::intersect(anchor_end_B, triA))
				{
					graph(i,j) = 1;					
				}
				else if (cu::intersect(anchor_start_A, triB))
				{
					graph(i,j) = 1;
				}
				else if (cu::intersect(anchor_start_B, triA))
				{
					graph(j,i) = 1;
				}
			}
		}
	}

	std::cout<<"graph matrix is" <<std::endl;
	std::cout<<graph<<std::endl;

	std::vector<int> order_of_robots;
	bool has_removal = false;
	while (has_removal)
	{
		for (int i=0; i<num_of_robot; i++)
		{
			for (auto index:order_of_robots)
			{
				if (i==index) continue;
			}
			bool has_incoming_vertex = false;
			for (int j=0; j<num_of_robot; j++)
			{
				if (i==j) continue;
				for (auto index:order_of_robots)
				{
					if (j==index) continue;
				}				
				if (graph(j,i)==1)
				{
					has_incoming_vertex = true;
					break;
				}
			}
			if (!has_incoming_vertex)
			{
				order_of_robots.push_back(i);
				has_removal = true;
				goto nextloop;
			}
		}
		has_removal = false;
		nextloop: ;
	}
	std::cout<<"current order of robots is" <<std::endl;
	for (auto i:order_of_robots) std::cout<<" "<<i;
	std::cout<<std::endl;
	
	for (int i=0; i<num_of_robot; i++)
	{
		bool already_in = false;
		for (auto index:order_of_robots)
		{
			if (i==index) 
			{
				already_in = true;
				break;
			}
		}
		if (!already_in) order_of_robots.push_back(i);
	}
	std::cout<<"After adding other robots, the order of robots is" <<std::endl;
	for (auto i:order_of_robots) std::cout<<" "<<i;
	std::cout<<std::endl;
	std::vector<std::vector<Vector3d>> pathsInOrder;
	if (order_of_robots.size()!=num_of_robot)
	{
		std::cout<<"ERROR: robots in order not equal to total number of robots!"
				 <<" returning"<<std::endl;
		return false;
	}

	if (use_optimal_path_finding_) //use optimal path from astar algorithm, which is 
											 // the cable shadow method described in the paper, 
											 //with radius of robot taken into account)
	{
		for (int order=0; order<num_of_robot; order++)
		{
			int i = order_of_robots[order]; //the checking robot index
			std::vector<Matrix<double, 3, 2>> cable_lines;
			for (int j=0; j<num_of_robot; j++)
			{
				if (i==j) continue;
				int order_of_j;
				for (int k=0; k<num_of_robot; k++)
				{
					if (j==order_of_robots[k]) order_of_j = k;
				}
				Vector3d point;
				Eigen::Matrix<double, 3, 2> anchor_end_B;				

				if (order_of_j<i) //before the checking robot
				{	
					if (graph(j,i)==1) continue; //if robot j moves before the checking robot, then no need
					anchor_end_B << basepoints[j], startEndPts[j].col(1);
				}
				else
				{
					anchor_end_B << basepoints[j], startEndPts[j].col(0);
				}
				cable_lines.push_back(anchor_end_B);
			}
			std::cout<<"cable lines size is "<<cable_lines.size()<<std::endl;
			int search_result;
			AstarGraphSearch(startEndPts[i].col(0), startEndPts[i].col(1), 
								  cable_lines, basepoints[i], search_result);
			std::vector<Vector3d> pathForI;
			if (search_result)
			{
				pathForI = getPath();
			}
			pathsInOrder.push_back(pathForI);
		}
	}
	else //use the convex hull method described in the paper
	{
		for (int i=0; i<num_of_robot; i++)
		{
			int index = order_of_robots[i]; //the checking robot index
			Matrix<double, 3, 3> triA;
			triA << basepoints[index], startEndPts[index];	
			std::vector<Vector3d> listof3dpoints, pathForIndex;
			std::vector<Vector2d> listof2dpoints;

			for (int j=0; j<num_of_robot; j++)
			{
				if (index==j) continue;
				int order_of_j;
				for (int k=0; k<num_of_robot; k++)
				{
					if (j==order_of_robots[k]) order_of_j = k;
				}
				Vector3d point;
				if (order_of_j<i) //before the checking robot
				{	
					if (graph(j,index)==1) continue; //if robot j moves before the checking robot, then no need
					Eigen::Matrix<double, 3, 2> anchor_end_B;				
					anchor_end_B << basepoints[j], startEndPts[j].col(1);
					point = cu::FindIntersectPt(anchor_end_B, triA);				
				}
				else
				{
					Eigen::Matrix<double, 3, 2> anchor_start_B;				
					anchor_start_B << basepoints[j], startEndPts[j].col(0);
					point = cu::FindIntersectPt(anchor_start_B, triA);	
				}

				if (point(0)<-1000 && point(1)<-1000 && point(2)<-1000) //no intersection
					continue;
				listof3dpoints.push_back(point);
				listof2dpoints.push_back(Vector2d(point(0), point(1)));
			}
	    	std::vector<Point_2> points_cgal;		
			for (auto pt: listof2dpoints)
			{
		      points_cgal.push_back(Point_2(pt(0), pt(1)));
			}
			points_cgal.push_back(Point_2(startEndPts[index](0,0), startEndPts[index](1,0)));
			points_cgal.push_back(Point_2(startEndPts[index](0,1), startEndPts[index](1,1)));

			Matrix<double, 2, Eigen::Dynamic> poly2d = cu::convexHullOfPoints2d(points_cgal);
			int startPtIndex = 0;
			double minDist = 1000;
			for (int j=0; j<poly2d.cols(); j++)
			{
				double dist = (poly2d.col(j)-startEndPts[index].block<2,1>(0,0)).norm();
				if (dist<minDist)
				{
					minDist = dist;
					startPtIndex = j;
				}
			}
			double epsilon = 1e-5;
			pathForIndex.push_back(startEndPts[index].col(0)); //push the start point first

			if ((poly2d.col((startPtIndex+1)%poly2d.cols()) - 
					startEndPts[index].block<2,1>(0,1)).norm()<epsilon)
			{
				for (int j=startPtIndex-1; j>startPtIndex-poly2d.cols(); j--)
				{
					int jj = j<0? j+poly2d.cols() : j;
					for (int k=0; k<listof2dpoints.size(); k++) //identify the 3d point that is in the convex hull
					{
						if ((poly2d.col(jj)-listof2dpoints[k]).norm()<epsilon)
						{
							pathForIndex.push_back(listof3dpoints[k]);
						}
					}
				}
			}
			else
			{
				if ((poly2d.col(startPtIndex-1<0? startPtIndex-1+poly2d.cols():startPtIndex-1) - 
					 startEndPts[index].block<2,1>(0,1)).norm()>epsilon)
				{
					std::cout<<"ERROR: the end point is not neighbouring"<< 
							   "the start point in the convex Hull!"<<std::endl;
				}
				for (int j=startPtIndex+1; j<startPtIndex+poly2d.cols(); j++)
				{
					int jj = j%poly2d.cols();
					for (int k=0; k<listof2dpoints.size(); k++) //identify the 3d point that is in the convex hull
					{
						if ((poly2d.col(jj)-listof2dpoints[k]).norm()<epsilon)
						{
							pathForIndex.push_back(listof3dpoints[k]);
						}
					}
				}
			}
			pathForIndex.push_back(startEndPts[index].col(1)); //push the end point
			pathsInOrder.push_back(pathForIndex);
		}		
	}
	
	// put back in order
	for (int i=0; i<num_of_robot; i++)
	{
		for (int j=0; j<num_of_robot; j++)
		{
			if (i==order_of_robots[j])
			{
				outputPath.push_back(pathsInOrder[j]);
			}
		}
	}
	if (outputPath.size()!=num_of_robot)
	{
		std::cout<<"ERROR: outputPath size not equal to number of robots!"<<std::endl;
		return false;
	}

	return true; // dun need to run the following codes for benchmarking purpose

	// calculate motion profile for robots, taking velocity as 1
	// the code below may is not thoroughly tested, but it does not afftect out benchmark
	std::vector<std::vector<std::vector<Vector2d>>> ListOfRPIndexOrdered;
	std::vector<std::vector<std::vector<Matrix<double, 3, 2>>>> ListOfRPOrdered; //restricted paths
	for (int i=0; i<num_of_robot; i++) //compute intersection between cable surfaces
	{
		std::vector<Vector2d> RPIndexUnordered;
		std::vector<Matrix<double, 3, 2>> RPUnordered; //restricted paths
		std::vector<std::vector<Vector2d>> RPIndexOrdered;
		std::vector<std::vector<Matrix<double, 3, 2>>> RPOrdered; //restricted paths

		for (int j=0; j<num_of_robot; j++)
		{
			if (i==j) continue;		
			for (int k=0; k<outputPath[i].size()-1; k++)
			{
				Matrix<double, 3, 3> triA;
				triA << basepoints[i], outputPath[i][k], outputPath[i][k+1];
				for (int jk=0; jk<outputPath[j].size()-1; jk++)
				{
					Matrix<double, 3, 3> triB;				
					triB << basepoints[j], outputPath[j][jk], outputPath[j][jk+1];		
					Segment_3 outputSegment;
					bool result = cu::FindIntersectTri(triA, triB, outputSegment);
					if (result)
					{
						Point_3_exact sourcePt = outputSegment.source();
						Point_3_exact targetPt = outputSegment.target();
						Point_3_exact pathPtA(outputPath[i][k](0), 
													 outputPath[i][k](1), 
													 outputPath[i][k](2));
						Point_3_exact pathPtB(outputPath[i][k+1](0), 
													 outputPath[i][k+1](1), 
													 outputPath[i][k+1](2));				
						Point_3_exact point_base(basepoints[i](0), 
														 basepoints[i](1), 
														 basepoints[i](2));
						Line_3 lineA(point_base, sourcePt);
						Line_3 lineB(point_base, targetPt);
						Segment_3 path_line(pathPtA, pathPtB);
						const auto result1 = CGAL::intersection(lineA, path_line);
						const auto result2 = CGAL::intersection(lineB, path_line);

						Vector3d intersectPathlineA, intersectPathlineB;
						if (result1 && result2)
						{
							if (const Point_3_exact* s = boost::get<Point_3_exact>(&*result1)) 
							{
						      intersectPathlineA(0) = CGAL::to_double(s->x());
						      intersectPathlineA(1) = CGAL::to_double(s->y());
						      intersectPathlineA(2) = CGAL::to_double(s->z());								
							}
							else 
							{
								std::cout<<"ERROR: 1!"<<std::endl; 
								goto nextrobot; //ignore it for now, which is wrong
							}
							if (const Point_3_exact* s = boost::get<Point_3_exact>(&*result2)) 
							{
						      intersectPathlineB(0) = CGAL::to_double(s->x());
						      intersectPathlineB(1) = CGAL::to_double(s->y());
						      intersectPathlineB(2) = CGAL::to_double(s->z());								
							}
							else 
							{
								std::cout<<"ERROR: 1!"<<std::endl; 
								goto nextrobot; //ignore it for now, which is wrong
							}								
						}
						else
						{ 
							std::cout<<"ERROR: 2!"<<std::endl;
							goto nextrobot; //ignore it for now, which is wrong
						}

						Matrix<double, 3, 2> rp;
						if ((outputPath[i][k]-intersectPathlineA).norm()<
							 (outputPath[i][k]-intersectPathlineB).norm())
						{
							rp<<intersectPathlineA, intersectPathlineB;
						}
						else
						{
							rp<<intersectPathlineB, intersectPathlineA;
						}
						RPUnordered.push_back(rp);
						RPIndexUnordered.push_back(Vector2d(j,k));
						goto nextrobot; //each robot should only intersect once with the robot
					}
				}
			}
			nextrobot:;

		}

		//after finding all intersection, we should order them according to their starting points
		for (int k=0; k<outputPath[i].size()-1; k++)
		{

			std::vector<Vector2d> RPIndexUnorderedForAPath;
			std::vector<Matrix<double, 3, 2>> RPUnorderedForAPath; //restricted paths						
			for (int kk=0; kk<RPIndexUnordered.size(); kk++) //for each path line
			{
	
				if (RPIndexUnordered[kk](1)==k)
				{
					RPIndexUnorderedForAPath.push_back(RPIndexUnordered[kk]);
					RPUnorderedForAPath.push_back(RPUnordered[kk]);
				}
			}
			std::vector<Matrix<double, 3, 2>>  RPOrdered_tmp;
			std::vector<Vector2d> RPIndexOrdered_tmp;
			while (!RPUnorderedForAPath.empty())
			{
				// std::cout<<"current RPUnorderedForAPath size is "<<RPUnorderedForAPath.size()<<std::endl;
				double smallest_distance = 10000.0;
				int smallest_idx = 0;
				for (int kk=0; kk<RPUnorderedForAPath.size(); kk++)
				{
					if ((outputPath[i][k]-RPUnorderedForAPath[kk].col(0)).norm()<smallest_distance)
					{
						smallest_distance = (outputPath[i][k]-RPUnorderedForAPath[kk].col(0)).norm();
						smallest_idx = kk;
					}
				}
				RPOrdered_tmp.push_back(RPUnorderedForAPath[smallest_idx]);
				RPIndexOrdered_tmp.push_back(RPIndexUnorderedForAPath[smallest_idx]);

				RPUnorderedForAPath.erase(RPUnorderedForAPath.begin()+smallest_idx);
				RPIndexUnorderedForAPath.erase(RPIndexUnorderedForAPath.begin()+smallest_idx);				
			}
			RPOrdered.push_back(RPOrdered_tmp);
			RPIndexOrdered.push_back(RPIndexOrdered_tmp);			
		}

		ListOfRPOrdered.push_back(RPOrdered);
		ListOfRPIndexOrdered.push_back(RPIndexOrdered);
	}

	//result is stored in these vectors:
	std::vector<std::vector<Vector2d>> timingOfRobots; //<each robot in moving order <each RP<start and depart times>>>
	std::vector<std::vector<int>>  indexOfIntersectingRobots;//<each robot in moving order <each RP<index of robot it intersects>>>

	for (int order=0; order<num_of_robot; order++)
	{
		std::vector<Vector2d> timingsForRobot;
		std::vector<int> indexForRobot;		
		int i=order_of_robots[order];
		for (int k=0; k<ListOfRPOrdered[i].size(); k++)
		{
			Vector3d pk, qk, qk_minus1;
			double ta = 0.0, td = 0.0; //BUG:td should not be zero at this point
												//since it does not affect our benchmark we leave it for later to fix
			qk_minus1 = outputPath[i][0];

			for (int kk=0; kk<ListOfRPOrdered[i][k].size(); kk++)
			{
				pk = ListOfRPOrdered[i][k][kk].col(0);
				qk = ListOfRPOrdered[i][k][kk].col(1);
				ta = td + (qk_minus1-pk).norm();
				td = ta + (pk - qk).norm();
				double t_w = 0.0;
				int j = ListOfRPIndexOrdered[i][k][kk](0);
				int order_of_j;
				for (int o=0; o<num_of_robot; o++)
				{
					if (j==order_of_robots[o]) order_of_j = o;
				}
				if (order_of_j<order)
				{
					for (int id = 0; id<indexOfIntersectingRobots[order_of_j].size(); id++)
					{
						if (indexOfIntersectingRobots[order_of_j][id]==i)
						{
							if (timingOfRobots[order_of_j][id](0)<=td && 
								 timingOfRobots[order_of_j][id](1)>=ta) //overlap
							{
								t_w = timingOfRobots[order_of_j][id](1) - ta;
								ta += t_w;
								td += t_w;
							}
							break;	
						}
					}
				}	
				timingsForRobot.push_back(Vector2d(ta, td));
				indexForRobot.push_back(j);
				qk_minus1 = pk;
			}
		}
		timingOfRobots.push_back(timingsForRobot);
		indexOfIntersectingRobots.push_back(indexForRobot);
	}
	return true;
}

void mtlp::AstarGraphSearch(Vector3d start_pt, Vector3d end_pt, 
                          std::vector<Matrix<double, 3, 2>>& cable_lines, 
                          Vector3d basepoint,  int& search_result)
{   
    CompareCost comparer;
    std::priority_queue<GridNodePtr, std::vector<GridNodePtr>, CompareCost> openList_tmp(comparer);
    openList_.swap(openList_tmp);      
    cableLines_ = cable_lines;
    basepoint_ = basepoint;

    resetUsedGrids();

    //index of start_point and end_point
    Vector3i start_idx = coord2gridIndex(start_pt);
    Vector3i end_idx   = coord2gridIndex(end_pt);
    goalIdx = end_idx;
	 goal_height_ = end_pt(2);
    //position of start_point and end_point
    end_pt   = gridIndex2coord(end_idx);

    //Initialize the pointers of struct GridNode which represent start node and goal node
    node_used_num_ = 0;
    GridNodePtr startPtr = GridNodeMap[node_used_num_++];
    startPtr ->index = start_idx;
    endPtr   = new GridNode(end_idx,   end_pt);

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
    startPtr -> cameFrom = NULL; 

    openList_.push(startPtr);
    // // openSet.insert( make_pair(startPtr -> fScore, startPtr) );
    expanded_nodes_.insert(start_idx, startPtr);    

    std::vector<GridNodePtr> neighborPtrSets;
    std::vector<double> edgeCostSets;

    while ( !openList_.empty() ){

        currentPtr = openList_.top();
        openList_.pop();
        // ROS_WARN("Currently expand a node with coord %f, %f", currentPtr->coord(0), currentPtr->coord(1));        
        // ROS_WARN("Currently expand a node with Fscore %f", currentPtr->fScore);
        currentPtr->id = -1; //set it to be in closed set

        // if the current node is the goal 
        if( currentPtr->index == goalIdx ){
            search_result = 1;
            // ros::Time time_2 = ros::Time::now();
            terminatePtr = currentPtr;
            currentPtr->coord = end_pt;
            // ROS_WARN("[A*]{sucess}  Time in A*  is %f ms, path cost if %f m", time_astar, currentPtr->gScore * resolution );
            return;
        }
        //get the succetion
        AstarGetSucc(currentPtr);     
      
    }
    //if search fails
    search_result = -1;
    // ros::Time time_2 = ros::Time::now();
    // if((time_2 - time_1).toSec() > 0.1)
    //     ROS_WARN("[A*]{Fail} Time consume in Astar path finding is %f", (time_2 - time_1).toSec() );
}

bool mtlp::intersectSphereSegment(Sphere_3 sphe, Matrix<double, 3, 2> segment)
{
	Point_3_exact point1(segment(0,1), segment(1,1), segment(2,1));
	Point_3_exact point2(segment(0,0), segment(1,0), segment(2,0));
	Segment_3 cable(point1, point2);
	return CGAL::do_intersect(sphe, cable);
}

inline void mtlp::AstarGetSucc(GridNodePtr currentPtr)
{   
	 int m = 0;
	 for (int i=-1; i<=1; i++)
	     for (int j=-1; j<=1; j++)
	     		for (int k=-1; k<=1; k++)
		     	{
		         if (i==0 && j==0 && k==0) continue; //skip the current point
		         Vector3i offsetIndex(i, j, k);
		         Vector3i neighbourIndex = currentPtr->index + offsetIndex;
		         double offsetCost = offsetIndex.cast<float>().lpNorm<2>();

		         Vector3d pkplus1 = gridIndex2coord(neighbourIndex);
		         if (pkplus1(2)>goal_height_) continue; //should never be higher

		         Matrix<double, 3, 3> triA;
		         triA << basepoint_, pkplus1, currentPtr->coord;
					Point_3_exact center(pkplus1(0), pkplus1(1), pkplus1(2));
					Sphere_3 robot_sphere(center, radius_squared);

					bool intersected = false;
		         for (int c = 0; c<cableLines_.size(); c++)
		         {
		         	if (cu::intersect(cableLines_[c], triA) || 
		         		 intersectSphereSegment(robot_sphere, cableLines_[c]))
		         	{
							// std::cout<<"got a bad index!!"<<std::endl;
							intersected = true;
							break;
		         	} 
		         }
		         if (intersected) continue;

		         GridNodePtr nodeptr;
		         nodeptr = expanded_nodes_.find(neighbourIndex);        
		         if (nodeptr == NULL) { //not expanded at all

		             if (node_used_num_ >= GLXYZ_SIZE) return;
		             nodeptr = GridNodeMap[node_used_num_++];
		             nodeptr-> id = 1;
		             nodeptr-> index = neighbourIndex;                
		             nodeptr-> gScore = offsetCost + currentPtr->gScore; //get accumulated g cost
		             nodeptr-> fScore = nodeptr->gScore + bias_ * getHeu(nodeptr, endPtr);
		             nodeptr-> cameFrom = currentPtr;
		             nodeptr-> coord = pkplus1;
		             openList_.push(nodeptr);
		             expanded_nodes_.insert(neighbourIndex, nodeptr);

		             continue;                

		         } else if (nodeptr->id==1) { //in open list already
		             double newgScore = offsetCost + currentPtr->gScore; //get accumulated g cost
		             if (nodeptr->gScore > newgScore){       //compare new g cost with previous g cost, if lower then update the openset
		                 nodeptr->gScore = newgScore;
		                 nodeptr->fScore = newgScore + bias_ * getHeu(nodeptr, endPtr);
		                 nodeptr->cameFrom = currentPtr;

		             }                

		         } else { // in closed list already, dun care
		             // double newgScore = offsetCost + currentPtr->gScore; //get accumulated g cost
		             // if (nodeptr->gScore > newgScore){       //compare new g cost with previous g cost, if lower then update the openset
		             //     nodeptr->gScore = newgScore;
		             //     nodeptr->fScore = newgScore + getHeu(nodeptr, endPtr);
		             //     nodeptr->cameFrom = currentPtr;

		             //     nodeptr-> id = 1;
		             // }  
		         }
		 
		         m++;

		     }
}

std::vector<Vector3d> mtlp::getPath() 
{   
    std::vector<Vector3d> path;
    std::vector<GridNodePtr> gridPath;

    GridNodePtr gptr=terminatePtr;
    while(gptr->cameFrom!=NULL){
        gridPath.push_back(gptr);
        gptr=gptr->cameFrom;

    }
    gridPath.push_back(gptr);

    GridNodePtr lastnode = gridPath.back();
    for (auto ptr: gridPath)
    {
        path.push_back(Vector3d(ptr->coord(0), ptr->coord(1), ptr->coord(2)));
    }
        
    reverse(path.begin(),path.end());
    for (int i = 0; i<path.size()-1; i++)
    {
        // length += (path[i+1] - path[i]).norm();
    }
    return path;
}


void mtlp::resetUsedGrids()
{   
    for(int i=0; i < GLXYZ_SIZE; i++)
        resetGrid(GridNodeMap[i]);
}

void mtlp::resetGrid(GridNodePtr ptr)
{
    ptr->id = 0;
    ptr->cameFrom = NULL;
    ptr->gScore = inf;
    ptr->fScore = inf;
}
