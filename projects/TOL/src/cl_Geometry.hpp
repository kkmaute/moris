/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Geometry.hpp
 *
 */

#ifndef SRC_TOOLS_CL_GEOMETRY_HPP_
#define SRC_TOOLS_CL_GEOMETRY_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

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
        static real
        vol_tetrahedron(Matrix< DDRMat >  & aCoord);

    };

}

#endif /* SRC_TOOLS_CL_GEOMETRY_HPP_ */

