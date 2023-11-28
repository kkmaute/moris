/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Surface_Mesh.hpp
 *
 */

#pragma once

#include "cl_SDF_Object.hpp"
#include "cl_GEN_Geometry.hpp"

namespace moris::ge
{
    class Surface_Mesh : public sdf::Object
    {

        // functions
        public:
            /**
             * Constructor. Creates vertices and facets from the designated input file
             * 
             * @param aFilePath .stl or .obj file that contains the surface mesh vertices and edges
             */
            Surface_Mesh( 
                const std::string& aFilePath, 
                const Matrix< DDRMat >& aOffsets,
                Geometry_Field_Parameters aParameters );

            /**
             * Gets potential ancestor nodes that are near the parent nodes.
             * 
             * @param aFirstParentGlobalCoordinates Global coordinates of the first parent node queued for intersection
             * @param aSecondParentGlobalCoordinates Global coordinates of the second parent node queued for intersection
             * @return Cell where each entry is a matrix with an ancestor node coordinate
             */
            Cell< Matrix< DDRMat > >
            get_candidate_intersection_nodes( 
                const Matrix< DDRMat >&         aFirstParentGlobalCoordinates,
                const Matrix< DDRMat >&         aSecondParentGlobalCoordinates );
    };
} // namespace moris::ge

