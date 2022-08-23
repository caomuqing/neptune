/* ----------------------------------------------------------------------------
 * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Jesus Tordesillas, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include "cgal_utils.hpp"
#include <CGAL/convex_hull_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/intersections.h>

// Convert several polyhedra to a vector that contains all the edges of all these polyhedra
mt::Edges cu::vectorGCALPol2edges(const ConvexHullsOfCurves& convexHulls)
{
  // See example here:
  // http://cgal-discuss.949826.n4.nabble.com/Take-the-triangles-of-a-polyhedron-td4460275.html

  // Other related questions:
  // https://doc.cgal.org/5.0/Triangulation_3/Triangulation_3_2for_loop_8cpp-example.html
  // https://stackoverflow.com/questions/4837179/getting-a-vertex-handle-from-an-edge-iterator

  mt::Edges all_edges;

  for (int index_curve = 0; index_curve < convexHulls.size(); index_curve++)  // for each curve
  {
    for (int i = 0; i < convexHulls[index_curve].size(); i++)  // for each interval along the curve
    {
      CGAL_Polyhedron_3 poly = convexHulls[index_curve][i];

      for (CGAL_Polyhedron_3::Edge_iterator w = poly.edges_begin(); w != poly.edges_end();
           ++w)  // for all the edges of that polyhedron
      {
        // std::cout << "First Vertex of the edge" << w->opposite()->vertex()->point() << std::endl;
        // std::cout << "Second Vertex of the edge" << w->vertex()->point() << std::endl;

        Eigen::Vector3d vertex1(w->opposite()->vertex()->point().x(), w->opposite()->vertex()->point().y(),
                                w->opposite()->vertex()->point().z());

        Eigen::Vector3d vertex2(w->vertex()->point().x(), w->vertex()->point().y(), w->vertex()->point().z());

        std::pair<Eigen::Vector3d, Eigen::Vector3d> edge;
        edge.first = vertex1;
        edge.second = vertex2;

        all_edges.push_back(edge);
      }
    }
  }

  return all_edges;
}

mt::ConvexHullsOfCurves_Std cu::vectorGCALPol2vectorStdEigen(ConvexHullsOfCurves& convexHulls)
{
  mt::ConvexHullsOfCurves_Std convexHulls_of_curves_std;

  // std::cout << "convexHulls.size()= " << convexHulls.size() << std::endl;

  for (int index_curve = 0; index_curve < convexHulls.size(); index_curve++)  // for each curve
  {
    mt::ConvexHullsOfCurve_Std convexHulls_of_curve_std;

    for (int i = 0; i < convexHulls[index_curve].size(); i++)  // for each interval along the curve
    {
      CGAL_Polyhedron_3 poly = convexHulls[index_curve][i];

      mt::Polyhedron_Std convexHull_std(3,
                                        poly.size_of_vertices());  // poly.size_of_vertices() is the number of vertexes
      // std::vector<Eigen::Vector3d> convexHull_std;
      int j = 0;
      for (CGAL_Polyhedron_3::Vertex_iterator v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
      {
        Eigen::Vector3d vertex(v->point().x(), v->point().y(), v->point().z());
        convexHull_std.col(j) = vertex;
        // convexHull_std.push_back(vertex);
        j = j + 1;
        // std::cout << v->point() << std::endl;
      }

      convexHulls_of_curve_std.push_back(convexHull_std);
    }

    convexHulls_of_curves_std.push_back(convexHulls_of_curve_std);
  }
  // std::cout << "convexHulls_of_curves_std.size()= " << convexHulls_of_curves_std.size() << std::endl;

  return convexHulls_of_curves_std;
}

vec_E<Polyhedron<3>> cu::vectorGCALPol2vectorJPSPol(ConvexHullsOfCurves& convex_hulls_of_curves)
{
  vec_E<Polyhedron<3>> vector_of_polyhedron_jps;

  for (auto convex_hulls_of_curve : convex_hulls_of_curves)
  {
    for (auto polyhedron_i : convex_hulls_of_curve)
    {
      vec_E<Hyperplane<3>> hyperplanes;
      for (auto it = polyhedron_i.planes_begin(); it != polyhedron_i.planes_end(); it++)
      {
        Vector_3 n = it->orthogonal_vector();
        Point_3 p = it->point();
        hyperplanes.push_back(
            Hyperplane<3>(Eigen::Vector3d(p.x(), p.y(), p.z()), Eigen::Vector3d(n.x(), n.y(), n.z())));
      }
      Polyhedron<3> polyhedron_jps(hyperplanes);
      vector_of_polyhedron_jps.push_back(polyhedron_jps);
      // std::cout << red << bold << "Size=" << vector_of_polyhedron_jps.size() << reset << std::endl;
    }
  }
  return vector_of_polyhedron_jps;
}

CGAL_Polyhedron_3 cu::convexHullOfPoints(const std::vector<Point_3>& points)
{
  // generate 3 points randomly on a sphere of radius 1.0
  // and copy them to a vector
  /*  CGAL::Random_points_in_sphere_3<Point_3, PointCreator> gen(2.0);
    std::vector<Point_3> points;
    CGAL::copy_n(gen, 6, std::back_inserter(points));*/

  // define object to hold convex hull
  CGAL::Object ch_object;

  // compute convex hull
  // std::cout << "Computing the convex hull CGAL for these points:" << std::endl;
  // for (auto points_i : points)
  // {
  //   std::cout << points_i.x() << ", " << points_i.y() << ", " << points_i.z() << std::endl;
  // }
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object);
  // std::cout << "convexHullCgal Computed!" << std::endl;
  // determine what kind of object it is
  // if (CGAL::object_cast<Segment_3>(&ch_object))
  //   std::cout << "convex hull is a segment " << std::endl;
  // else if (CGAL::object_cast<CGAL_Polyhedron_3>(&ch_object))
  //   std::cout << "convex hull is a polyhedron " << std::endl;
  // else
  //   std::cout << "convex hull error!" << std::endl;

  CGAL_Polyhedron_3 poly = *CGAL::object_cast<CGAL_Polyhedron_3>(&ch_object);

  std::transform(poly.facets_begin(), poly.facets_end(), poly.planes_begin(),
                 cu::Plane_equation());  // Compute the planes

  return poly;

  /*
CGAL::set_pretty_mode(std::cout);
// std::copy(poly.planes_begin(), poly.planes_end(), std::ostream_iterator<Plane_3>(std::cout, "\n"));
*/
}

mt::Polygon_Std cu::convexHullOfPoints2d(const std::vector<Point_2>& points)
{
  // std::vector<Point_2>::iterator out;
  std::vector<Point_2> out;
  // std::iota(indices.begin(), indices.end(),0);
  // CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
  //                     Convex_hull_traits_2(CGAL::make_property_map(points)));
  CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(out));

  mt::Polygon_Std convexHull_std2d(2, out.size());
  int j=0;
  for( Point_2 pt : out){
    convexHull_std2d.col(j)= Eigen::Vector2d(pt.x(), pt.y());
    j++;
  }
  // std::cout<<"convex hull 2 returns points: "<<convexHull_std2d<<std::endl;
  return convexHull_std2d;
}

bool cu::intersect(const Eigen::Matrix<double, 2, Eigen::Dynamic>& pointsA,
               const Eigen::Matrix<double, 2, Eigen::Dynamic>& pointsB)
{
  std::vector<Point_2_exact> in1;
  std::vector<Point_2_exact> in2;
  for (int i = 0; i<pointsA.cols(); i++)
  {
    in1.push_back(Point_2_exact(pointsA(0,i), pointsA(1,i)));
  }  

  for (int i = 0; i<pointsB.cols(); i++)
  {
    in2.push_back(Point_2_exact(pointsB(0,i), pointsB(1,i)));
  }  
  // std::cout<<"solving for polygons"<<std::endl;

  Polygon_2 poly1(in1.begin(), in1.end());
  Polygon_2 poly2(in2.begin(), in2.end());
  // std::cout<<"solving for intersections"<<std::endl;

  return CGAL::do_intersect(poly1, poly2);
}

bool cu::intersectSeg(const Eigen::Matrix<double, 2, 2>& segA, 
                      const Eigen::Matrix<double, 2, 2>& segB)
{
  std::vector<Point_2_exact> in1;
  std::vector<Point_2_exact> in2;
  for (int i = 0; i<2; i++)
  {
    in1.push_back(Point_2_exact(segA(0,i), segA(1,i)));
  }  

  for (int i = 0; i<2; i++)
  {
    in2.push_back(Point_2_exact(segB(0,i), segB(1,i)));
  }  

  Segment_2 seg1(in1[0], in1[1]);
  Segment_2 seg2(in2[0], in2[1]);

  return CGAL::do_intersect(seg1, seg2);  
}

bool cu::intersectTri(const Eigen::Matrix<double, 3, 3>& triangleA,
                      const Eigen::Matrix<double, 3, 3>& triangleB)
{
  std::vector<Point_3_exact> in1;
  std::vector<Point_3_exact> in2;
  for (int i = 0; i<triangleA.cols(); i++)
  {
    in1.push_back(Point_3_exact(triangleA(0,i), triangleA(1,i), triangleA(2,i)));
  }

  for (int i = 0; i<triangleB.cols(); i++)
  {
    in2.push_back(Point_3_exact(triangleB(0,i), triangleB(1,i), triangleB(2,i)));
  }
  Triangle_3 triA(in1[0], in1[1], in1[2]);
  Triangle_3 triB(in2[0], in2[1], in2[2]);

  return CGAL::do_intersect(triA, triB);
}

bool cu::intersect(const Eigen::Matrix<double, 3, 2>& lineSeg,
               const Eigen::Matrix<double, 3, 3>& triangle)
{
  std::vector<Point_3_exact> in1;
  std::vector<Point_3_exact> in2;
  for (int i = 0; i<lineSeg.cols(); i++)
  {
    in1.push_back(Point_3_exact(lineSeg(0,i), lineSeg(1,i), lineSeg(2,i)));
  }

  for (int i = 0; i<triangle.cols(); i++)
  {
    in2.push_back(Point_3_exact(triangle(0,i), triangle(1,i), triangle(2,i)));
  }
  Segment_3 seg(in1[0], in1[1]);
  Triangle_3 tri(in2[0], in2[1], in2[2]);

  return CGAL::do_intersect(seg, tri);
}

Eigen::Vector3d cu::FindIntersectPt(const Eigen::Matrix<double, 3, 2>& lineSeg,
               const Eigen::Matrix<double, 3, 3>& triangle)
{
  std::vector<Point_3_exact> in1;
  std::vector<Point_3_exact> in2;
  for (int i = 0; i<lineSeg.cols(); i++)
  {
    in1.push_back(Point_3_exact(lineSeg(0,i), lineSeg(1,i), lineSeg(2,i)));
  }

  for (int i = 0; i<triangle.cols(); i++)
  {
    in2.push_back(Point_3_exact(triangle(0,i), triangle(1,i), triangle(2,i)));
  }
  Segment_3 seg(in1[0], in1[1]);
  Triangle_3 tri(in2[0], in2[1], in2[2]);

  const auto result = CGAL::intersection(seg, tri);
  Eigen::Vector3d outputPt(-10000, -10000, -10000);

  if (result) 
  {
    if (const Point_3_exact* s = boost::get<Point_3_exact>(&*result)) 
    {
      outputPt[0] = CGAL::to_double(s->x());
      outputPt[1] = CGAL::to_double(s->y());
      outputPt[2] = CGAL::to_double(s->z());
    } 
    else 
    {
      const Segment_3* p = boost::get<Segment_3 >(&*result);
      std::cout << "segment is on the triangle! wrong!" << std::endl;
    }
  }
  else
  {
    std::cout << "not getting intersection! wrong!" << std::endl;
  }  
  return outputPt;
}

bool cu::FindIntersectTri(const Eigen::Matrix<double, 3, 3>& triangleA,
                      const Eigen::Matrix<double, 3, 3>& triangleB,
                      Segment_3& outputSegment)
{
  std::vector<Point_3_exact> in1;
  std::vector<Point_3_exact> in2;
  for (int i = 0; i<triangleA.cols(); i++)
  {
    in1.push_back(Point_3_exact(triangleA(0,i), triangleA(1,i), triangleA(2,i)));
  }

  for (int i = 0; i<triangleB.cols(); i++)
  {
    in2.push_back(Point_3_exact(triangleB(0,i), triangleB(1,i), triangleB(2,i)));
  }
  Triangle_3 triA(in1[0], in1[1], in1[2]);
  Triangle_3 triB(in2[0], in2[1], in2[2]);

  const auto result = CGAL::intersection(triA, triB);

  if (result) 
  {
    if (const Point_3_exact* s = boost::get<Point_3_exact>(&*result)) 
    {
      return false;
    } 
    else if (const Segment_3* s = boost::get<Segment_3>(&*result)) 
    {
      const Segment_3* p = boost::get<Segment_3 >(&*result);
      outputSegment = *p;
      return true;
    }
    else
    {
      std::cout<<"FindIntersectTri: getting intersection other than point or segment, strange!"<<std::endl;
      return false;
    }
  }
  else
  {
    return false;
  }  

}

bool cu::check_inside(Point_2 pt, Point_2 *pgn_begin, Point_2 *pgn_end, K traits)
{
  switch(CGAL::bounded_side_2(pgn_begin, pgn_end,pt, traits)) {
    case CGAL::ON_BOUNDED_SIDE :
      return true;
      break;
    case CGAL::ON_BOUNDARY:
      return true;
      break;
    case CGAL::ON_UNBOUNDED_SIDE:
      return false;
      break;
  }
}