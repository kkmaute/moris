/*
 * cl_XTK_Enums.hpp
 *
 *  Created on: Nov 29, 2016
 *      Author: doble
 */

#ifndef SRC_XTK_CL_XTK_ENUMS_HPP_
#define SRC_XTK_CL_XTK_ENUMS_HPP_

enum class Decomposition
{
    // first decomposition(conformal) _ second decomposition(nonconformal)

    REGULAR_HIER, // regular subdivision of hex 8 followed by heirarchical subdivision
    TREE
};

enum class TemplateType
{
    REGULAR_SUBDIVISION_HEX8, // Topology created using a regularized subdivision (for generate_templated_mesh)
    QUAD_4,                    // Topology created using a Quad 4 template topology
    TET_4,
    HIERARCHY_TET4,
    HIERARCHY_TET4_3N,  // 3 node intersection pattern
    HIERARCHY_TET4_4Na, // 4 node intersection with high and low across from each other
    HIERARCHY_TET4_4Nb, // 4 node intersection pattern with high and MH acoss from each other
    HIERARCHY_TET4_4Nc, // 4 node intersection pattern with high and ML acoss from each other
};

enum class RequestType
{
    // Please order in the following manner NC, C, T then alphabetical
    // Because Nonconformal check happens first, Conformal check happens second, Tests happen last
    // NC - specifies a nonconformal request
    // C  - specifies a conformal request
    // T  - specifies
    NC_REGULAR_SUBDIVISION_HEX8,  // Nonconformal and a regular subdivision template will be used
    C_HIERARCHY_TET4,             // Conformal and a hierarchy template will be used
    T_QUAD4_EDGE_NODE             // Parallel test in test/src/cl_model.hpp edge requested on a node
};

#endif /* SRC_XTK_CL_XTK_ENUMS_HPP_ */
