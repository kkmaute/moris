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

// XTKL: XTK Includes
#include "fn_mesh_flood_fill.hpp"

namespace xtk
{
    /*
     * Performs local child mesh flood fill operation and returns the elemental subphases
     * see test case xtk/fn_flood_fill.cpp
     */
    inline moris::Matrix< moris::IndexMat >
    local_child_mesh_flood_fill( Child_Mesh& aChildMesh )
    {
        // Get number of elements in the child mesh
        moris::size_t tNumElements = aChildMesh.get_num_entities( EntityRank::ELEMENT );

        // Specify dummy value as maximum moris::size_t val
        moris::size_t tMax = std::numeric_limits< moris::moris_index >::max();

        // Maximum number of phases
        moris::size_t tNumPhases = 2;

        // Allocate space for active elements
        moris::Matrix< moris::IndexMat > tActiveElements( 1, tNumElements );

        for ( moris::size_t iE = 0; iE < tNumElements; iE++ )
        {
            // Add element index to active element list
            ( tActiveElements )( 0, iE ) = iE;
        }

        // Mark all elements as included
        moris::Matrix< moris::IndexMat > tIncludedElementMarker( 1, tNumElements, 1 );

        // maximum subphase
        moris_index tMaxSubphase = 0;

        // Run flood fill Algorithm
        moris::Matrix< moris::IndexMat > tElementSubphase = flood_fill( aChildMesh.get_element_to_element(),
                aChildMesh.get_element_phase_indices(),
                tActiveElements,
                tIncludedElementMarker,
                tNumPhases,
                tMax,
                tMaxSubphase,
                true );

        return tElementSubphase;
    }
}    // namespace xtk

#endif /* SRC_XTK_FN_LOCAL_CHILD_MESH_FLOOD_FILL_HPP_ */
