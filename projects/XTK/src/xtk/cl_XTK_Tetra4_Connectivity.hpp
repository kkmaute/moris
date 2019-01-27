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
    class Tetra4_Connectivity
    {
    public:
        //TODO: [@keenan.doble@colorado.edu ] Figure out how make these constexp such that we know they
        // are resolved completely at compile time.

        /*
         * Face information
         */
        // Node ordinals associated with each face of the tetrahedral
        // FIXME: The ordering of face ordinal 1 seems to be CW rather than CCW
        const moris::Matrix<moris::IndexMat> mTetra4FaceMap = {{0, 1, 3}, {2, 1, 3}, {0, 2, 3}, {0, 2, 1}};

        // Nodes on edge to cross for outward facing normals of each face ordinal
        // This is just one possible way to compute the correct cross product.
        // Cell - Face ordinal
        // Col 0 - Nodes on edge 0
        // Col 1 - Nodes on edge 1
        const moris::Cell<moris::Matrix<moris::IndexMat>> mTetra4NodesForOutwardNormal = {{{0,0},{1,3}}, /*Nodes on edges to cross to get outward normal for face 0*/
                                                                                          {{1,1},{2,3}}, /*Nodes on edges to cross to get outward normal for face 1*/
                                                                                          {{0,3},{3,2}}, /*Nodes on edges to cross to get outward normal for face 2*/
                                                                                          {{0,2},{2,1}}};

        // Node ordinals on each edge
        const moris::Matrix<moris::IndexMat> mTetra4EdgeMap = {{0, 1},{1, 2},{0, 2},{0, 3},{1, 3},{2, 3}};


    };
}


#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_TETRA4_CONNECTIVITY_HPP_ */
