/*
 * cl_MTK_Enums.hpp
 *
 *  Created on: Jul 25, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_ENUMS_HPP_
#define SRC_MESH_CL_MTK_ENUMS_HPP_

#include "assert.hpp"
#include "typedefs.hpp"

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

        enum class Geometry_Type
        {
            POINT, // point
            LINE,  // 1D line or curve
            QUAD,  // rectangle
            TRI,   // triangle
            HEX,   // quadrangle
            TET,   // tetrahedron
            PENTA,  // pentahedron
            UNDEFINED
        };

//------------------------------------------------------------------------------

        enum class Interpolation_Order
        {
            CONSTANT,
            LINEAR,
            QUADRATIC,
            SERENDIPITY,
            CUBIC,
            UNDEFINED
        };

        enum class Master_Slave
        {
            MASTER,
            SLAVE,
            UNDEFINED
        };

        enum class Primary_Void
        {
            PRIMARY,
            VOID,
            INTERP,
            UNDEFINED
        };


//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */


#endif /* SRC_MESH_CL_MTK_ENUMS_HPP_ */
