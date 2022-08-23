//  Original C implementation created by Igor Kroitor on 29/12/15.
//  CPP implementation with Eigen library created Muqing Cao on aug 2021.

#pragma once

#include <Eigen/Dense>
#include <stdio.h>
//-----------------------------------------------------------------------------
// Gilbert-Johnson-Keerthi (GJK) collision detection algorithm in 2D
// http://www.dyn4j.org/2010/04/gjk-gilbert-johnson-keerthi/
// http://mollyrocket.com/849
//-----------------------------------------------------------------------------

using namespace Eigen;

namespace gjk  // cgal utils
{
//-----------------------------------------------------------------------------
// Basic vector arithmetic operations
Vector2d perpendicular (Vector2d v);

//-----------------------------------------------------------------------------
// Triple product expansion is used to calculate perpendicular normal vectors 
// which kinda 'prefer' pointing towards the Origin in Minkowski space

Vector2d tripleProduct (Vector2d a, Vector2d b, Vector2d c);
//-----------------------------------------------------------------------------
// This is to compute average center (roughly). It might be different from
// Center of Gravity, especially for bodies with nonuniform density,
// but this is ok as initial direction of simplex search in GJK.

Vector2d averagePoint(const Matrix<double, 2, Dynamic> &mat);

//-----------------------------------------------------------------------------
// Get furthest vertex along a certain direction

int indexOfFurthestPoint (const Matrix<double, 2, Dynamic>& vertices, Vector2d d);

//-----------------------------------------------------------------------------
// Minkowski sum support function for GJK

Vector2d support (const Matrix<double, 2, Dynamic>& vertices1, 
                  const Matrix<double, 2, Dynamic>& vertices2, Vector2d d);

//-----------------------------------------------------------------------------
// The GJK yes/no test


bool collision(const Matrix<double, 2, Dynamic>& vertices1,
          const Matrix<double, 2, Dynamic>& vertices2);

}