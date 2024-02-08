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

#include "cl_GEN_Field_Analytic.hpp"
#include "cl_Param_List.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Library_IO.hpp"
#include "cl_GEN_Node_Manager.hpp"

namespace moris::gen
{
    /**
     * Creates an instance of the specified Field class and returns a shared pointer to it.
     *
     * @param aFieldParameterList Parameter list for creating a field class
     * @param aADVs ADV Vector
     * @param aLibrary Pointer to library for loading user-defined functions
     * @return Pointer to a Field object
     */
    std::shared_ptr< Field > create_field(
            ParameterList                         aFieldParameterList,
            Matrix< DDRMat >&                     aADVs,
            Vector< std::shared_ptr< Field > >      aFieldDependencies = {},
            std::shared_ptr<Library_IO>           aLibrary = nullptr,
            mtk::Mesh*                            aMTKMesh = nullptr,
            uint                                  aIndex = 0 );
}
