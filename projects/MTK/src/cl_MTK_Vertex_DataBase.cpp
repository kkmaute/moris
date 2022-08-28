/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_DataBase.cpp
 *
 */

#include "cl_MTK_Vertex_DataBase.hpp"
#include <typeinfo>
#include "fn_TOL_Capacities.hpp"

//#include "cl_MTK_Interpolation_Mesh_Analysis.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    Vertex_DataBase::Vertex_DataBase( moris_index aVertexIndex,
        mtk::Mesh* const&                         aMesh ) :
        mVertexIndex( aVertexIndex ),
        mMesh( aMesh )
    {
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Vertex_DataBase::get_coords() const
    {
        real* tCoordsPointer =  mMesh->get_vertex_coords_ptr( mVertexIndex );
        // use advanced constructor
        Matrix< DDRMat > tCoordMatrix(tCoordsPointer, 1, mMesh->get_spatial_dim(), false, true );

        return tCoordMatrix;
    }

    //------------------------------------------------------------------------------

    moris_id
    Vertex_DataBase::get_id() const
    {
        // ask mesh for the id
        return mMesh->get_entity_id( EntityRank::NODE, mVertexIndex );
    }

    //------------------------------------------------------------------------------

    moris_index
    Vertex_DataBase::get_index() const
    {
        return mVertexIndex;
    }

    //------------------------------------------------------------------------------

    moris_index
    Vertex_DataBase::get_owner() const
    {
        return mMesh->get_entity_owner( EntityRank::NODE, mVertexIndex );
    }

    //------------------------------------------------------------------------------

    Vertex_Interpolation*
    Vertex_DataBase::get_interpolation( const uint aBSplineMeshIndex )
    {
        // get the first interpolation vertex
        Vertex_Interpolation** tFirstIPVertex = mMesh->get_vertex_interpolation( mVertexIndex );
        uint tLocalOrder = mMesh->get_local_mesh_index( aBSplineMeshIndex );

        return *( tFirstIPVertex + tLocalOrder );
    }

    //------------------------------------------------------------------------------

    const Vertex_Interpolation*
    Vertex_DataBase::get_interpolation( const uint aBSplineMeshIndex ) const
    {
        // get the first interpolation vertex
        Vertex_Interpolation** tFirstIPVertex = mMesh->get_vertex_interpolation( mVertexIndex );
        uint tLocalOrder = mMesh->get_local_mesh_index( aBSplineMeshIndex );

        return *( tFirstIPVertex + tLocalOrder );
    }

    //------------------------------------------------------------------------------

    bool
    Vertex_DataBase::has_interpolation( const uint aBSplineMeshIndex )
    {
        // FIXME: new method needed for this function
        MORIS_ERROR( 0, "Vertex_DataBase::has_interpolation() - no implemneted in the Vertex_DataBase class" );
        return false;
    }

    //------------------------------------------------------------------------------

    size_t
    Vertex_DataBase::capacity()
    {
        size_t tCapacity = 0;

        tCapacity += sizeof( mVertexIndex );
        tCapacity += sizeof( mMesh );

        return tCapacity;
    }

    //------------------------------------------------------------------------------
}// namespace moris::mtk
