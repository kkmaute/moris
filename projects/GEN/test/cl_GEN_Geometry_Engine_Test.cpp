/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry_Engine_Test.cpp
 *
 */

#include "cl_GEN_Geometry_Engine_Test.hpp"

#include <utility>

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    Geometry_Engine_Test::Geometry_Engine_Test(
            mtk::Interpolation_Mesh*   aMesh,
            Geometry_Engine_Parameters aParameters )
            : Geometry_Engine( aMesh, std::move( aParameters ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Node_Manager& Geometry_Engine_Test::get_node_manager()
    {
        return mNodeManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< Geometry > Geometry_Engine_Test::get_geometry(uint aGeometryIndex)
    {
        return mGeometries(aGeometryIndex);
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr<Property> Geometry_Engine_Test::get_property(uint aPropertyIndex)
    {
        return mProperties(aPropertyIndex);
    }

    //--------------------------------------------------------------------------------------------------------------

}
