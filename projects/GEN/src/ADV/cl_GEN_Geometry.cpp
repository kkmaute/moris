/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry.cpp
 *
 */

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Derived_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Geometry::Geometry( Node_Manager& aNodeManager )
            : mNodeManager( aNodeManager )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Geometry::register_node_manager( Node_Manager& aNodeManager )
    {
        mNodeManager = aNodeManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Geometry::get_geometric_region(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        if ( aNodeIndex < mNodeManager.get_number_of_base_nodes() )
        {
            // This is a base node, so just pass on its index and coordinates
            return this->get_base_geometric_region( aNodeIndex, aNodeCoordinates );
        }
        else
        {
            // More information may be needed, so get the derived node from the node manager
            return this->get_derived_geometric_region( mNodeManager.get_derived_node( aNodeIndex ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Geometry::get_derived_geometric_region( Derived_Node* aDerivedNode )
    {
        return this->get_base_geometric_region( aDerivedNode->get_index(), aDerivedNode->get_coordinates() );
    }

    //--------------------------------------------------------------------------------------------------------------

}
