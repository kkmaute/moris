/**
 * fn_determine_cell_topology.hpp  
 * 
 *  Created on: Nov 16, 2021 
 *      Author: Nils Wunsch
 */
#ifndef SRC_fn_determine_cell_topology
#define SRC_fn_determine_cell_topology

#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_MTK_Interpolation_Enum_Int_Conversion.hpp"

namespace xtk
{
    //------------------------------------------------------------------------------

    enum CellTopology
    determine_cell_topology(
        uint aNumSpatialDims,
        uint aPolynomialOrder,
        enum CellShape aCellShape )
    {
        // --------------------------------------------------
        // 2D
        // --------------------------------------------------

        // get number of spatial dimensions and decide on cell topology
        if ( aNumSpatialDims == 2 )
        {
            if ( aCellShape == CellShape::SIMPLEX )
            {
                switch ( aPolynomialOrder )
                {
                    case 1:
                    {
                        return CellTopology::TRI3;
                        break;
                    }

                    case 2:
                    {
                        return CellTopology::TRI6;
                        break;
                    }

                    case 3:
                    {
                        return CellTopology::TRI10;
                        break;
                    }
                
                    default:
                    {
                        MORIS_ERROR( false, "determine_cell_topology() - cell topology can only be determined for polynomial orders 1,2, or 3" );
                        return CellTopology::INVALID;
                        break;
                    }
                }
            }
            else if ( aCellShape == CellShape::RECTANGULAR )
            {
                switch ( aPolynomialOrder )
                {
                    case 1:
                    {
                        return CellTopology::QUAD4;
                        break;
                    }

                    case 2:
                    {
                        return CellTopology::QUAD9;
                        break;
                    }

                    case 3:
                    {
                        return CellTopology::QUAD16;
                        break;
                    }
                
                    default:
                    {
                        MORIS_ERROR( false, "determine_cell_topology() - cell topology can only be determined for polynomial orders 1,2, or 3" );
                        return CellTopology::INVALID;
                        break;
                    }
                }
            }
            else
            {
                MORIS_ERROR( false, "determine_cell_topology() - cell topology can only be determined for rectangular or simplex elements" );
                return CellTopology::INVALID;
            }
        }

        // --------------------------------------------------
        // 3D
        // --------------------------------------------------

        else if ( aNumSpatialDims == 3 )
        {
            if ( aCellShape == CellShape::SIMPLEX )
            {
                switch ( aPolynomialOrder )
                {
                    case 1:
                    {
                        return CellTopology::TET4;
                        break;
                    }

                    case 2:
                    {
                        return CellTopology::TET10;
                        break;
                    }

                    // case 3:
                    // {
                    //     return CellTopology::TET20;
                    //     break;
                    // }
                
                    default:
                    {
                        MORIS_ERROR( false, "determine_cell_topology() - cell topology can only be determined for polynomial orders 1,2, or 3" );
                        return CellTopology::INVALID;
                        break;
                    }
                }
            }
            else if ( aCellShape == CellShape::RECTANGULAR )
            {
                switch ( aPolynomialOrder )
                {
                    case 1:
                    {
                        return CellTopology::HEX8;
                        break;
                    }

                    case 2:
                    {
                        return CellTopology::HEX27;
                        break;
                    }

                    case 3:
                    {
                        return CellTopology::HEX64;
                        break;
                    }
                
                    default:
                    {
                        MORIS_ERROR( false, "determine_cell_topology() - cell topology can only be determined for polynomial orders 1,2, or 3" );
                        return CellTopology::INVALID;
                        break;
                    }
                }
            }
            else
            {
                MORIS_ERROR( false, "determine_cell_topology() - cell topology can only be determined for rectangular or simplex elements" );
                return CellTopology::INVALID;
            }
        }

        // --------------------------------------------------

        else
        {
            MORIS_ERROR( false, "determine_cell_topology() - function only works for 2D or 3D" );
            return CellTopology::INVALID;
        }
    }

    //------------------------------------------------------------------------------

    enum CellTopology
    determine_cell_topology(
        uint aNumSpatialDims,
        enum mtk::Interpolation_Order aInterpolationOrder,
        enum CellShape aCellShape )
    {
        // call the above function with the enum replaced by integer
        return determine_cell_topology( aNumSpatialDims, mtk::ip_order_enum_to_uint( aInterpolationOrder ), aCellShape );
    }

    //------------------------------------------------------------------------------
}

#endif /* fn_determine_cell_topology.hpp */
