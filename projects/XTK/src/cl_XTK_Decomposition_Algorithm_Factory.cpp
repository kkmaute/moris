#include "cl_XTK_Decomposition_Algorithm_Factory.hpp"
#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Regular_Subdivision_Interface.hpp"
#include "cl_XTK_Node_Hierarchy_Interface.hpp"

    // free function
std::shared_ptr<xtk::Decomposition_Algorithm>
xtk::create_decomposition_algorithm(enum Subdivision_Method aSubdivisionMethod)
{
    switch (aSubdivisionMethod)
    {
    case Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8 :
        return std::make_shared<xtk::Regular_Subdivision_Interface>( );
        break;

    case Subdivision_Method::C_HIERARCHY_TET4 :
        return std::make_shared<xtk::Node_Hierarchy_Interface>( );
        break;
    
    default:
        MORIS_ERROR(0,"Decomposition algorithm corresponding to provided enum not implemented");
        return std::make_shared<xtk::Regular_Subdivision_Interface>( );
        break;
    }

}

