/*
 * cl_Geometry.hpp
 *
 *  Created on: Apr 7, 2017
 *      Author: gleim
 */

#ifndef SRC_TOOLS_CL_GEOMETRY_HPP_
#define SRC_TOOLS_CL_GEOMETRY_HPP_

#include "cl_Mat.hpp" // LNA/src

namespace moris
{
    class Geometry
    {
    public:
        /*
         * Geometry for a general Tetrahedron
         *
         * @param[in] aCoord .... Volume of 4 points: 3 coloumns for the basis directions and 4 rows for the Coordinates
         * @param[out] vol   .... Volume of the tetrahedron
         */
        static moris::real
        vol_tetrahedron(moris::Mat<moris::real>  & aCoord);

    };

}

#endif /* SRC_TOOLS_CL_GEOMETRY_HPP_ */
