/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_GEN_create_field.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_GEN_ADV_Manager.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Library_IO.hpp"
#include "cl_GEN_Node_Manager.hpp"

namespace moris::gen
{
    /**
     * Creates an instance of the specified Field class and returns a shared pointer to it.
     *
     * @param aFieldParameterList Parameter list for creating a field class
     * @param aADVManager ADV manager, constaining an ADV vector and bounds
     * @param aFieldDependencies Other fields that the created field depends on
     * @param aLibrary Pointer to library for loading user-defined functions
     * @param aMTKMesh MTK mesh
     * @return Pointer to a Field object
     */
    std::shared_ptr< Field > create_field(
            const Parameter_List&                aFieldParameterList,
            ADV_Manager&                         aADVManager,
            Vector< std::shared_ptr< Field > >   aFieldDependencies = {},
            const std::shared_ptr< Library_IO >& aLibrary           = nullptr,
            mtk::Mesh*                           aMTKMesh           = nullptr );

    /**
     * Gets the number of ADVs that will be created based on the given parameter list.
     *
     * @param aFieldParameterList Field parameter list
     * @return Number of active ADVs
     */
    Vector< char > get_active_parameter_ids( const Parameter_List& aFieldParameterList );
}
