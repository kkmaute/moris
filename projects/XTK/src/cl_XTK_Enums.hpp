/*
 * cl_XTK_Enums.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_ENUMS_HPP_
#define SRC_XTK_CL_XTK_ENUMS_HPP_

#include <string>
enum class TemplateType
{
    REGULAR_SUBDIVISION_HEX8, // Topology created using a regularized subdivision (for generate_templated_mesh)
    REGULAR_SUBDIVISION_QUAD4, // Topology created using a regular subdivision of Quad 4
    TRI_3,
    QUAD_4,                    // Topology created using a Quad 4 template topology
    TET_4,                    // Standard tet 4 topology
    HEX_8,
    HIERARCHY_TET4,
    HIERARCHY_TET4_3N,  // 3 node intersection pattern
    HIERARCHY_TET4_4Na, // 4 node intersection with high and low across from each other
    HIERARCHY_TET4_4Nb, // 4 node intersection pattern with high and MH across from each other
    HIERARCHY_TET4_4Nc, // 4 node intersection pattern with high and ML acoss from each other
    HIERARCHY_TET4_2,
    CONFORMAL_TRI3,
    BISECTED_TET4,      // A tet4 split into two
    INVALID_TEMPLATE_TYPE
};

enum class Subdivision_Method
{
    // Please order in the following manner NC, C, T then alphabetical
    // Because Nonconformal check happens first, Conformal check happens second, Tests happen last
    // NC - specifies a nonconformal request
    // C  - specifies a conformal request
    // T  - specifies a Test method
    NC_REGULAR_SUBDIVISION_HEX8,  // Nonconformal and a regular subdivision template will be used
    NC_REGULAR_SUBDIVISION_QUAD4,
    C_HIERARCHY_TET4,             // Conformal and a hierarchy template will be used
    C_TRI3,             // Conformal tri 3 mesh  will be constructed
    NO_METHOD
};

inline
const std::string get_enum_str(enum Subdivision_Method aSubdivisionEnum)
{
    switch (aSubdivisionEnum)
    {
       case Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8: return "NC_REGULAR_SUBDIVISION_HEX8";
       case Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4: return "NC_REGULAR_SUBDIVISION_QUAD4";
       case Subdivision_Method::C_HIERARCHY_TET4: return "C_HIERARCHY_TET4";
       case Subdivision_Method::C_TRI3: return "C_TRI3";
       case Subdivision_Method::NO_METHOD: return "NO_METHOD";
       default: return "invalid subdivision method";
    }
}

enum class Topology_Type
{
    EDGE, // Edge with 2 Node
    TRI_3,
    QUAD_4,
    TET_4,
    TET_10,
    HEXA_8, // hexahedron with 8 nodes topology

};

#endif /* SRC_XTK_CL_XTK_ENUMS_HPP_ */
