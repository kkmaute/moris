/*
 * cl_HMR_Implementation.cpp
 *
 *  Created on: Feb 9, 2018
 *      Author: gleim
 */

#include "cl_HMR_Implementation.hpp"
using namespace moris;

uint
HMR_Implementation::get_elem_topology_num_edges(
        uint aElemId ) const
{
    MORIS_ASSERT( mHMRMesh.give_modeldim() >= 2 && mHMRMesh.give_modeldim() <=3, " Edges are only in 2D and 3D available" );
    uint tNumberElementEdges = 0;
    if( mHMRMesh.give_modeldim() == 2 )
        tNumberElementEdges = 4;
    else if( mHMRMesh.give_modeldim() == 3 )
        tNumberElementEdges = 12;
    return tNumberElementEdges;
}

//--------------------------------------------------------------------------------

uint
HMR_Implementation::get_face_topology_num_edges(
        uint aFaceId ) const
{
    MORIS_ASSERT( mHMRMesh.give_modeldim() >= 2 && mHMRMesh.give_modeldim() <=3, " Edges are only in 2D and 3D available" );
    uint tNumberFaceEdges = 0;
    if( mHMRMesh.give_modeldim() == 2 )
    {
        tNumberFaceEdges = 2;
    }
    else if( mHMRMesh.give_modeldim() == 3 )
    {
        tNumberFaceEdges = 4;
    }
    return tNumberFaceEdges;
}

//--------------------------------------------------------------------------------

Mat<uint>
HMR_Implementation::get_edges_connected_to_element(
        uint const aElementId ) const
{
    //Determine edges of element
    Mat<uint> tEdgesOfElement = mHMRMesh.give_element_edges( aElementId );
    //Reordering (2D: 0 = bottom, 1 = right, 2 = top, 3 = left, 3D: 0-3 same and 4 = back, 5 = front)
    Mat<uint> tEdgesReordered;
    //Numbering comes from cubit side numbering
    if( mHMRMesh.give_modeldim() == 2 )
    {
        tEdgesReordered.set_size(4,1,0);
        tEdgesReordered( 0 ) = tEdgesOfElement(0);
        tEdgesReordered( 1 ) = tEdgesOfElement(3);
        tEdgesReordered( 2 ) = tEdgesOfElement(1);
        tEdgesReordered( 3 ) = tEdgesOfElement(2);
    }
    else if( mHMRMesh.give_modeldim() == 3 )
    {
        tEdgesReordered.set_size(12,1,0);
        tEdgesReordered( 0  ) = tEdgesOfElement(0);
        tEdgesReordered( 1  ) = tEdgesOfElement(5);
        tEdgesReordered( 2  ) = tEdgesOfElement(1);
        tEdgesReordered( 3  ) = tEdgesOfElement(4);
        tEdgesReordered( 4  ) = tEdgesOfElement(8);
        tEdgesReordered( 5  ) = tEdgesOfElement(9);
        tEdgesReordered( 6  ) = tEdgesOfElement(11);
        tEdgesReordered( 7  ) = tEdgesOfElement(10);
        tEdgesReordered( 8  ) = tEdgesOfElement(2);
        tEdgesReordered( 9  ) = tEdgesOfElement(7);
        tEdgesReordered( 10 ) = tEdgesOfElement(3);
        tEdgesReordered( 11 ) = tEdgesOfElement(6);
    }
    return tEdgesReordered;
}

Mat<uint>
HMR_Implementation::get_faces_connected_to_element(
        uint const aElementId ) const
{
    //Determine faces of element
    Mat<uint> tFacesOfElement = mHMRMesh.give_element_faces( aElementId );
    //Reordering (2D: 0 = bottom, 1 = right, 2 = top, 3 = left, 3D: 0-3 same and 4 = back, 5 = front)
    Mat<uint> tFacesReordered;
    if( mHMRMesh.give_modeldim() == 2 )
    {
        tFacesReordered.set_size(4,1,0);
        tFacesReordered( 0 ) = tFacesOfElement(2);
        tFacesReordered( 1 ) = tFacesOfElement(1);
        tFacesReordered( 2 ) = tFacesOfElement(3);
        tFacesReordered( 3 ) = tFacesOfElement(0);
    }
    else if( mHMRMesh.give_modeldim() == 3 )
    {
        tFacesReordered.set_size(6,1,0);
        tFacesReordered( 0 ) = tFacesOfElement(2);
        tFacesReordered( 1 ) = tFacesOfElement(1);
        tFacesReordered( 2 ) = tFacesOfElement(3);
        tFacesReordered( 3 ) = tFacesOfElement(0);
        tFacesReordered( 4 ) = tFacesOfElement(4);
        tFacesReordered( 5 ) = tFacesOfElement(5);
    }
    return tFacesReordered;
}
