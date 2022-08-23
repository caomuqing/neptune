/* ----------------------------------------------------------------------------
 * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Jesus Tordesillas, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once

#include "mader_types.hpp"
#include <Eigen/Dense>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Convex_hull_traits_3.h>
#include <decomp_geometry/polyhedron.h>
#include <CGAL/convex_hull_2.h>
// #include <CGAL/Convex_hull_traits_adapter_2.h>
// #include <CGAL/property_map.h>
// #include <CGAL/ch_graham_andrew.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Sphere_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel K2;

typedef CGAL::Convex_hull_traits_3<K> Traits;
typedef Traits::Polyhedron_3 CGAL_Polyhedron_3;
//typedef K::Segment_3 Segment_3;
typedef K::Plane_3 Plane_3;
// define point creator
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef CGAL::Creator_uniform_3<double, Point_3> PointCreator;
typedef CGAL::Polygon_2<K2> Polygon_2;

typedef std::vector<CGAL_Polyhedron_3> ConvexHullsOfCurve;
typedef std::vector<ConvexHullsOfCurve> ConvexHullsOfCurves;

typedef K::Point_2 Point_2;
typedef K2::Point_2 Point_2_exact;
typedef K2::Point_3 Point_3_exact;
typedef CGAL::Segment_3<K2> Segment_3;
typedef CGAL::Triangle_3<K2> Triangle_3;
typedef CGAL::Segment_2<K2> Segment_2;
typedef CGAL::Line_3<K2> Line_3;
typedef CGAL::Sphere_3<K2> Sphere_3;

// typedef CGAL::Convex_hull_traits_adapter_2<K,
//           CGAL::Pointer_property_map<Point_2>::type > Convex_hull_traits_2;

namespace cu  // cgal utils
{
struct Plane_equation
{
  template <class Facet>
  typename Facet::Plane_3 operator()(Facet& f)
  {
    typename Facet::Halfedge_handle h = f.halfedge();
    typedef typename Facet::Plane_3 Plane;
    return Plane(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
  }
};

mt::ConvexHullsOfCurves_Std vectorGCALPol2vectorStdEigen(ConvexHullsOfCurves& convexHulls);

vec_E<Polyhedron<3>> vectorGCALPol2vectorJPSPol(ConvexHullsOfCurves& convex_hulls_of_curves);

CGAL_Polyhedron_3 convexHullOfPoints(const std::vector<Point_3>& points);

mt::Edges vectorGCALPol2edges(const ConvexHullsOfCurves& convexHulls);

mt::Polygon_Std convexHullOfPoints2d(const std::vector<Point_2>& points);

bool intersect(const Eigen::Matrix<double, 2, Eigen::Dynamic>& pointsA,
               const Eigen::Matrix<double, 2, Eigen::Dynamic>& pointsB);

bool intersectSeg(const Eigen::Matrix<double, 2, 2>& segA, 
                  const Eigen::Matrix<double, 2, 2>& segB);

bool intersectTri(const Eigen::Matrix<double, 3, 3>& triangleA,
                      const Eigen::Matrix<double, 3, 3>& triangleB);

bool intersect(const Eigen::Matrix<double, 3, 2>& lineSeg,
               const Eigen::Matrix<double, 3, 3>& triangle);

Eigen::Vector3d FindIntersectPt(const Eigen::Matrix<double, 3, 2>& lineSeg,
               const Eigen::Matrix<double, 3, 3>& triangle);

bool FindIntersectTri(const Eigen::Matrix<double, 3, 3>& triangleA,
                      const Eigen::Matrix<double, 3, 3>& triangleB, Segment_3& outputSegment);

bool check_inside(Point_2 pt, Point_2 *pgn_begin, Point_2 *pgn_end, K traits);

}  // namespace cu