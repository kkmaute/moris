/*
 * cl_XTK_Tetra4_Connectivity.hpp
 *
 *  Created on: Jan 16, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_TETRA4_CONNECTIVITY_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_TETRA4_CONNECTIVITY_HPP_

#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace xtk
{
    namespace Tetra4_Connectivity
    {
    /*!
     * Node ordinal to face map
     */
    inline
    moris::Matrix<moris::IndexMat>
    get_node_to_face_map()
    {
        return {{0, 1, 3}, {2, 1, 3}, {0, 2, 3}, {0, 2, 1}};
    }

    inline
    moris::Matrix<moris::IndexMat>
    get_node_to_edge_map()
    {
        return {{0, 1},{1, 2},{0, 2},{0, 3},{1, 3},{2, 3}};
    }

    inline
    moris::Matrix<moris::IndexMat>
    get_node_to_face_map(moris::uint aSideOrdinal)
    {
        switch (aSideOrdinal)
        {
            case(0):
            {
                return {{0, 1, 3}};
                break;
            }
            case(1):
            {
                return {{2, 1, 3}};
                break;
            }
            case(2):
            {
                return {{0, 2, 3}};
                break;
            }

            case(3):
            {
                return {{0, 2, 1}};
                break;
            }
            default:
                MORIS_ASSERT(0,"Invalid side ordinal specified");
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
                return {{0,0},{1,3}};
                break;
            }
            case(1):
            {
                return {{1,1},{2,3}};
                break;
            }
            case(2):
            {
                return {{0,3},{3,2}};
                break;
            }

            case(3):
            {
                return {{0,2},{2,1}};
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


#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_TETRA4_CONNECTIVITY_HPP_ */
