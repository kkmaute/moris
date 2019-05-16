/*
 * cl_MTK_HEX8_CONNECTIVITY.HPP
 *
 *  Created on: May 14, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_HEX8_CONNECTIVITY_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_HEX8_CONNECTIVITY_HPP_


#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"



namespace moris
{
    namespace Hex8
    {
    /*!
     * Node ordinal to face map
     */
    inline
    moris::Matrix<moris::IndexMat>
    get_node_to_face_map()
    {
        return {{0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {0,4,7,3},{0,3,2,1},{4,5,6,7}};
    }

    inline
    moris::Matrix<moris::IndexMat>
    get_node_to_face_map(moris::uint aSideOrdinal)
    {
        switch (aSideOrdinal)
        {
            case(0):
            {
                return {{0, 1, 5, 4}};
                break;
            }
            case(1):
            {
                return {{1, 2, 6, 5}};
                break;
            }
            case(2):
            {
                return {{2, 3, 7, 6}};
                break;
            }

            case(3):
            {
                return {{0, 4, 7, 3}};
                break;
            }
            case(4):
            {
                return {{0, 3, 2, 1}};
                break;
            }
            case(5):
            {
                return {{4, 5, 6, 7}};
                break;
            }

            default:
                MORIS_ERROR(0,"Invalid side ordinal specified");
                return moris::Matrix<moris::IndexMat>(0,0);
                break;
        }
    }


    }
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_HEX8_CONNECTIVITY_HPP_ */
