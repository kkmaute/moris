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

}
