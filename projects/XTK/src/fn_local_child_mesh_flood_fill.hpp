/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_local_child_mesh_flood_fill.hpp
 *
 */

#ifndef SRC_XTK_FN_LOCAL_CHILD_MESH_FLOOD_FILL_HPP_
#define SRC_XTK_FN_LOCAL_CHILD_MESH_FLOOD_FILL_HPP_

// XTKL: Linalg Includes
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_Matrix.hpp"

// MTKL : Mesh Includes
#include "fn_MTK_mesh_flood_fill.hpp"

namespace moris::xtk
{
    /*
     * Performs local child mesh flood fill operation and returns the elemental subphases
     * see test case xtk/fn_flood_fill.cpp
     */
    inline Matrix< IndexMat >
    local_child_mesh_flood_fill( Child_Mesh& aChildMesh )
    {
        // Get number of elements in the child mesh
        moris::size_t tNumElements = aChildMesh.get_num_entities( mtk::EntityRank::ELEMENT );

        // Specify dummy value as maximum moris::size_t val
        moris::size_t tMax = std::numeric_limits< moris::moris_index >::max();

        // Allocate space for active elements
        Vector< moris_index > tActiveElements( tNumElements );

        for ( moris::size_t iE = 0; iE < tNumElements; iE++ )
        {
            // Add element index to active element list
            ( tActiveElements )( iE ) = iE;
        }

        // Mark all elements as included
        Vector< moris_index > tIncludedElementMarker( tNumElements, 1 );

        // maximum subphase
        moris_index tMaxSubphase = 0;

        // Run flood fill Algorithm
        Matrix< IndexMat > tElementSubphase = mtk::flood_fill( 
                aChildMesh.get_element_to_element(),
                aChildMesh.get_element_phase_indices(),
                tActiveElements,
                tIncludedElementMarker,
                tMax,
                tMaxSubphase,
                true );

        return tElementSubphase;
    }
}    // namespace moris::xtk

#endif /* SRC_XTK_FN_LOCAL_CHILD_MESH_FLOOD_FILL_HPP_ */
