/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Surface_Mesh.cpp
 *
 */

#include "cl_GEN_Surface_Mesh.hpp"

namespace moris::ge {
    Surface_Mesh::Surface_Mesh( 
        const std::string&          aFileName,
        const Matrix< DDRMat >&     aOffsets,
        Geometry_Field_Parameters aParameters )
    : Object( aFileName, aOffsets )
    {
    }


}