//  Original C implementation created by Igor Kroitor on 29/12/15.
//  CPP implementation with Eigen library created Muqing Cao on aug 2021.

#include <gjk.hpp>
#include <Eigen/Dense>
#include <stdio.h>
//-----------------------------------------------------------------------------
// Gilbert-Johnson-Keerthi (GJK) collision detection algorithm in 2D
// http://www.dyn4j.org/2010/04/gjk-gilbert-johnson-keerthi/
// http://mollyrocket.com/849
//-----------------------------------------------------------------------------

using namespace Eigen;

//-----------------------------------------------------------------------------
// Basic vector arithmetic operations
Vector2d gjk::perpendicular (Vector2d v) {
    return Vector2d(v(1), -v(0)); 
}

//-----------------------------------------------------------------------------
// Triple product expansion is used to calculate perpendicular normal vectors 
// which kinda 'prefer' pointing towards the Origin in Minkowski space

Vector2d gjk::tripleProduct (Vector2d a, Vector2d b, Vector2d c) {

    return b * a.dot(c) - a * b.dot(c);
}
//-----------------------------------------------------------------------------
// This is to compute average center (roughly). It might be different from
// Center of Gravity, especially for bodies with nonuniform density,
// but this is ok as initial direction of simplex search in GJK.

Vector2d gjk::averagePoint(const Matrix<double, 2, Dynamic> &mat){
    Vector2d avg = mat.rowwise().mean();
    return avg;
}

//-----------------------------------------------------------------------------
// Get furthest vertex along a certain direction

int gjk::indexOfFurthestPoint (const Matrix<double, 2, Dynamic>& vertices, Vector2d d) {
    
    double maxProduct = d.dot(vertices.col(0));
    int index = 0;
    for (int i = 1; i < vertices.cols(); i++) {
        double product = d.dot(vertices.col(i));
        if (product > maxProduct) {
            maxProduct = product;
            index = i;
        }
    }
    return index;
}

//-----------------------------------------------------------------------------
// Minkowski sum support function for GJK

Vector2d gjk::support (const Matrix<double, 2, Dynamic>& vertices1, 
                  const Matrix<double, 2, Dynamic>& vertices2, Vector2d d) {

    // get furthest point of first body along an arbitrary direction
    int i = indexOfFurthestPoint (vertices1, d);
    
    // get furthest point of second body along the opposite direction
    int j = indexOfFurthestPoint (vertices2, -d);

    // subtract (Minkowski sum) the two points to see if bodies 'overlap'
    return vertices1.col(i) - vertices2.col(j);
}

//-----------------------------------------------------------------------------
// The GJK yes/no test


bool gjk::collision(const Matrix<double, 2, Dynamic>& vertices1,
          const Matrix<double, 2, Dynamic>& vertices2) {
    
    int index = 0; // index of current vertex of simplex
    Vector2d a, b, c, d, ao, ab, ac, abperp, acperp;
    Matrix<double, 2, 3> simplex;
    
    Vector2d position1 = averagePoint (vertices1); // not a CoG but
    Vector2d position2 = averagePoint (vertices2); // it's ok for GJK )

    // initial direction from the center of 1st body to the center of 2nd body
    d = position1 - position2;
    
    // if initial direction is zero – set it to any arbitrary axis (we choose X)
    if ((d(0) == 0) && (d(1) == 0))
        d(0) = 1.0;
    
    // set the first support as initial point of the new simplex
    simplex.col(0) = support (vertices1, vertices2, d);
    a = simplex.col(0);
    
    if (a.dot(d) <= 0)
        return false; // no collision
    
    d = -a; // The next search direction is always towards the origin, so the next search direction is negate(a)
    
    while (1) {
        
        a = simplex.col(++index) = support (vertices1, vertices2, d);
        
        if (a.dot(d) <= 0)
            return false; // no collision
        
        ao = -a; // from point A to Origin is just negative A
        
        // simplex has 2 points (a line segment, not a triangle yet)
        if (index < 2) {
            b = simplex.col(0);
            ab = b- a; // from point A to B
            d = tripleProduct (ab, ao, ab); // normal to AB towards Origin
            if (d.norm() == 0)
                d = perpendicular (ab);
            continue; // skip to next iteration
        }
        
        b = simplex.col(1);
        c = simplex.col(0);
        ab = b - a; // from point A to B
        ac = c - a; // from point A to C
        
        acperp = tripleProduct (ab, ac, ac);
        
        if (acperp.dot(ao) >= 0) {
            
            d = acperp; // new direction is normal to AC towards Origin
            
        } else {
            
            abperp = tripleProduct (ac, ab, ab);
            
            if (abperp.dot(ao) < 0)
                return true; // collision
            
            simplex.col(0) = simplex.col(1); // swap first element (point C)

            d = abperp; // new direction is normal to AB towards Origin
        }
        
        simplex.col(1) = simplex.col(2); // swap element in the middle (point B)
        --index;
    }
    
    return false;
}