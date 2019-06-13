/*
 * cl_MTK_Hex64_CONNECTIVITY.HPP
 *
 *  Created on: May 14, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_HEX64_CONNECTIVITY_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_HEX64_CONNECTIVITY_HPP_


#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"



namespace moris
{
    namespace Hex64
    {
    /*!
     * Node ordinal to face map
     */
    inline
    moris::Matrix<moris::IndexMat>
    get_node_to_face_map()
    {
        return {{ 0, 1, 5, 4,  8,  9, 16, 17, 25, 24, 13, 12, 36, 37, 38, 39 },
                { 1, 2, 6, 5, 14, 15, 20, 21, 29, 28, 17, 16, 44, 45, 46, 47 },
                { 2, 3, 7, 6, 18, 19, 22, 23, 31, 30, 21, 20, 48, 49, 50, 51 },
                { 0, 4, 7, 3, 12, 13, 26, 27, 23, 22, 11, 10, 40, 41, 42, 43 },
                { 0, 3, 2, 1, 10, 11, 19, 18, 15, 14,  9,  8, 32, 33, 34, 35 },
                { 4, 5, 6, 7, 24, 25, 28, 29, 30, 31, 27, 26, 52, 53, 54, 55 }};
    }

    inline
    moris::Matrix<moris::IndexMat>
    get_node_to_face_map(moris::uint aSideOrdinal)
    {
        switch (aSideOrdinal)
        {
            case(0):
            {
                return {{ 0, 1, 5, 4,  8,  9, 16, 17, 25, 24, 13, 12, 36, 37, 38, 39 }};
                break;
            }
            case(1):
            {
                return {{ 1, 2, 6, 5, 14, 15, 20, 21, 29, 28, 17, 16, 44, 45, 46, 47 }};
                break;
            }
            case(2):
            {
                return {{ 2, 3, 7, 6, 18, 19, 22, 23, 31, 30, 21, 20, 48, 49, 50, 51 }};
                break;
            }

            case(3):
            {
                return {{ 0, 4, 7, 3, 12, 13, 26, 27, 23, 22, 11, 10, 40, 41, 42, 43 }};
                break;
            }
            case(4):
            {
                return {{ 0, 3, 2, 1, 10, 11, 19, 18, 15, 14,  9,  8, 32, 33, 34, 35 }};
                break;
            }
            case(5):
            {
                return {{ 4, 5, 6, 7, 24, 25, 28, 29, 30, 31, 27, 26, 52, 53, 54, 55 }};
                break;
            }

            default:
                MORIS_ERROR(0,"Invalid side ordinal specified");
                return moris::Matrix<moris::IndexMat>(0,0);
                break;
        }
    }

    /*!
     * Nodes on edge to cross for outward facing normals of each face ordinal
     * This is just one possible way to compute the correct cross product.
     *   Cell - Face ordinal
     *   Col 0 - Nodes on edge 0
     *   Col 1 - Nodes on edge 1
     * @param aSideOrdinal Side ordinal of element
     * @return Edge nodes to use for outward normal
     */

    inline
    moris::Matrix<moris::IndexMat>
    get_node_map_outward_normal(moris::uint aSideOrdinal)
    {
        switch (aSideOrdinal)
        {
            case(0):
            {
                return {{0,1},{1,5}};
                break;
            }
            case(1):
            {
                return {{1,2},{2,6}};
                break;
            }
            case(2):
            {
                return {{2,3},{3,7}};
                break;
            }

            case(3):
            {
                return {{7,3},{3,0}};
                break;
            }

            case(4):
            {
                return {{0,3},{3,2}};
                break;
            }
            case(5):
            {
                return {{4,5},{5,6}};
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

#endif /* PROJECTS_MTK_SRC_CL_MTK_Hex64_CONNECTIVITY_HPP_ */
