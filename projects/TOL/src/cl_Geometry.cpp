/*
 * cl_Geometry.cpp
 *
 *  Created on: Apr 10, 2017
 *      Author: gleim
 */

#include "cl_Geometry.hpp"
#include "algorithms.hpp"
#include "linalg.hpp"

moris::real
moris::Geometry::vol_tetrahedron(moris::Mat<moris::real>  & aCoord)
{
    //explanation: www.colorado.edu/engineering/Aerospace/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
    moris::Mat<moris::real> J(4,4,1.0); // 4 Coordinates are needed for the Jacobinan-matrix.
    // The order of the coordinates is irrelevant.

    J( 1, 0 ) = aCoord(0,0); J( 1, 1 ) = aCoord(1,0); J( 1, 2 ) = aCoord(2,0); J( 1, 3 ) = aCoord(3,0);
    J( 2, 0 ) = aCoord(0,1); J( 2, 1 ) = aCoord(1,1); J( 2, 2 ) = aCoord(2,1); J( 2, 3 ) = aCoord(3,1);
    J( 3, 0 ) = aCoord(0,2); J( 3, 1 ) = aCoord(1,2); J( 3, 2 ) = aCoord(2,2); J( 3, 3 ) = aCoord(3,2);

    moris::real volume = std::abs(moris::det( J )/6); // Volume = 1/6*det(J)

    return volume;
}
