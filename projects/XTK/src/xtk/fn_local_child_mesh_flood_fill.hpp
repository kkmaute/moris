/*
 * fn_local_child_mesh_flood_fill.hpp
 *
 *  Created on: May 16, 2018
 *      Author: doble
 */

#ifndef SRC_XTK_FN_LOCAL_CHILD_MESH_FLOOD_FILL_HPP_
#define SRC_XTK_FN_LOCAL_CHILD_MESH_FLOOD_FILL_HPP_


// XTKL: Linalg Includes
#include "cl_XTK_Child_Mesh.hpp"
#include "linalg/cl_XTK_Matrix.hpp"


// XTKL: XTK Includes
#include "xtk/fn_mesh_flood_fill.hpp"

namespace xtk
{
/*
 * Performs local child mesh flood fill operation and returns the elemental subphases
 * see test case xtk/fn_flood_fill.cpp
 */
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
Mat<Integer,Integer_Matrix>
local_child_mesh_flood_fill(Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> & aChildMesh)
{
    // Get number of elements in the child mesh
    Integer tNumElements = aChildMesh.get_num_entities(EntityRank::ELEMENT);

    // Specify dummy value as maximum integer val
    Integer tMax = std::numeric_limits<Integer>::max();

    // Maximum number of element neighbors
    Integer tMaxNeighbors = 4;

    // Maximum number of phases
    Integer tNumPhases = 2;

    // Allocate space for active elements
    Mat<Integer, Integer_Matrix> tActiveElements(1,tNumElements);

    for(Integer iE = 0; iE<tNumElements; iE++)
    {
        // Add element index to active element list
        (tActiveElements)(0,iE) = iE;
    }

    // Mark all elements as included
    Mat<Integer, Integer_Matrix> tIncludedElementMarker(1,tNumElements,1);

    // Run flood fill Algorithm
    Mat<Integer, Integer_Matrix> tElementSubphase = flood_fill( aChildMesh.get_element_to_element(),
                                                                aChildMesh.get_element_phase_indices(),
                                                                tActiveElements,
                                                                tIncludedElementMarker,
                                                                tNumPhases,
                                                                tMax,
                                                                true);

    return tElementSubphase;
}
}


#endif /* SRC_XTK_FN_LOCAL_CHILD_MESH_FLOOD_FILL_HPP_ */
