/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Decomposition_Algorithm_Factory.cpp
 *
 */

#include "cl_XTK_Decomposition_Algorithm_Factory.hpp"
#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Regular_Subdivision_Interface.hpp"
#include "cl_XTK_Node_Hierarchy_Interface.hpp"
#include "cl_XTK_Octree_Interface.hpp"
#include "cl_XTK_Elevate_Order_Interface.hpp"

// free function
std::shared_ptr<xtk::Decomposition_Algorithm>
xtk::create_decomposition_algorithm(
  enum Subdivision_Method aSubdivisionMethod,
  moris::ParameterList&   aParameterList)
{
  switch (aSubdivisionMethod)
  {
  case Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4:
    return std::make_shared<xtk::Regular_Subdivision_Interface>( aParameterList, mtk::CellTopology::QUAD4 );
    break;

  case Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8:
    return std::make_shared<xtk::Regular_Subdivision_Interface>( aParameterList, mtk::CellTopology::HEX8 );
    break;

  case Subdivision_Method::C_TRI3:
    return std::make_shared<xtk::Node_Hierarchy_Interface>( aParameterList );
    // return std::make_shared<xtk::Node_Hierarchy_Interface>( aParameterList, mtk::CellTopology::TRI3 );
    break;

  case Subdivision_Method::C_HIERARCHY_TET4:
    return std::make_shared<xtk::Node_Hierarchy_Interface>( aParameterList );
    // return std::make_shared<xtk::Node_Hierarchy_Interface>( aParameterList, mtk::CellTopology::TET4 );
    break;

  case Subdivision_Method::NC_OCTREE:
    return std::make_shared<xtk::Octree_Interface>( aParameterList );
    break;

  case Subdivision_Method::P_ELEVATE_ORDER_TRI3_TRI6:
    return std::make_shared<xtk::Elevate_Order_Interface>( aParameterList, Subdivision_Method::P_ELEVATE_ORDER_TRI3_TRI6  );
    break;

//  case Subdivision_Method::P_ELEVATE_ORDER_TRI3_TRI10:
//     return std::make_shared<xtk::Elevate_Order_Interface>( aParameterList, Subdivision_Method::P_ELEVATE_ORDER_TRI3_TRI10  );
//     break;

  case Subdivision_Method::P_ELEVATE_ORDER_TET4_TET10:
    return std::make_shared<xtk::Elevate_Order_Interface>( aParameterList, Subdivision_Method::P_ELEVATE_ORDER_TET4_TET10  );
    break;

//   case Subdivision_Method::P_ELEVATE_ORDER_TET4_TET20:
//     return std::make_shared<xtk::Elevate_Order_Interface>( aParameterList, Subdivision_Method::P_ELEVATE_ORDER_TET4_TET20  );
//     break;

  default:
    MORIS_ERROR( 0, "Decomposition algorithm corresponding to provided enum not implemented" );
    return std::make_shared<xtk::Regular_Subdivision_Interface>( aParameterList, mtk::CellTopology::HEX8 );
    break;
  }
}

