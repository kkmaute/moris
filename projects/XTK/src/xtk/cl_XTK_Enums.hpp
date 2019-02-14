/*
 * cl_XTK_Enums.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_ENUMS_HPP_
#define SRC_XTK_CL_XTK_ENUMS_HPP_


// Enums in this header
// 1.) TemplateType
// 2.) Subdivision_Method
// 3.) Topology_Type
// 4.) Phase_Table_Structure
// 5.) Enrichment_Method

enum class TemplateType
{
    REGULAR_SUBDIVISION_HEX8, // Topology created using a regularized subdivision (for generate_templated_mesh)
    QUAD_4,                    // Topology created using a Quad 4 template topology
    TET_4,                    // Standard tet 4 topology
    HEX_8,
    STACKED_2_TET4, // 2 Stacked tet 4s
    HIERARCHY_TET4,
    HIERARCHY_TET4_3N,  // 3 node intersection pattern
    HIERARCHY_TET4_4Na, // 4 node intersection with high and low across from each other
    HIERARCHY_TET4_4Nb, // 4 node intersection pattern with high and MH acoss from each other
    HIERARCHY_TET4_4Nc, // 4 node intersection pattern with high and ML acoss from each other
    HIERARCHY_TET4_2,
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
    C_HIERARCHY_TET4,             // Conformal and a hierarchy template will be used
    T_GENERATE_TET_MESH_FROM_HEX  // For test cases, intersects all elements provided
};

inline
const char* get_enum_str(enum Subdivision_Method aSubdivisionEnum)
{
    switch (aSubdivisionEnum)
    {
       case Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8: return "NC_REGULAR_SUBDIVISION_HEX8";
       case Subdivision_Method::C_HIERARCHY_TET4: return "C_HIERARCHY_TET4";
       case Subdivision_Method::T_GENERATE_TET_MESH_FROM_HEX: return "T_GENERATE_TET_MESH_FROM_HEX";
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

enum class Phase_Table_Structure
{
    EXP_BASE_2
};

enum class Matrix_Backend
{
    EIGEN,
    TEUCHOS
};




#endif /* SRC_XTK_CL_XTK_ENUMS_HPP_ */
