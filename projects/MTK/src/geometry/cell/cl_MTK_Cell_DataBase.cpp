/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_DataBase.cpp
 *
 */

#include "cl_MTK_Cell_DataBase.hpp"

#include "cl_MTK_Mesh_Core.hpp"
#include "fn_TOL_Capacities.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    Cell_DataBase::Cell_DataBase( mtk::Cell&         aCell,
            std::shared_ptr< moris::mtk::Cell_Info > aCellInfo,
            moris_index                              aCellIndex2,
            mtk::Mesh*                               aMesh )
            : Cell( aCell.get_id(), aCell.get_index(), aCell.get_owner(), aCellInfo )
            , mBaseCell( aCell.get_base_cell() )
            , mCellIndex2( aCellIndex2 )
            , mMesh( aMesh )
    {
    }

    //------------------------------------------------------------------------------

    Cell_DataBase::Cell_DataBase( mtk::Cell& aCell,
            moris_index                      aCellIndex2,
            mtk::Mesh*                       aMesh )
            : Cell( aCell.get_id(), aCell.get_index(), aCell.get_owner(), aCell.get_cell_info_sp() )
            , mCellIndex2( aCellIndex2 )
            , mMesh( aMesh )
    {
    }

    //------------------------------------------------------------------------------

    Cell_DataBase::Cell_DataBase( moris_index aCellIndex2,
            mtk::Mesh*                        aMesh )
            : mCellIndex2( aCellIndex2 )
            , mMesh( aMesh )
    {
    }

    //------------------------------------------------------------------------------

    moris::Vector< Vertex* >
    Cell_DataBase::get_vertex_pointers() const
    {
        // get number of vertices to obtain how much need to be marched in the coordinate list
        uint tNumVertices = this->get_number_of_vertices();

        // initialize the vertices cell
        moris::Vector< mtk::Vertex* > tVertexPtrs;

        Vertex** tVertices = mMesh->get_cell_vertices( mCellIndex2 );
        // std::copy( mVertices, mVertices + tNumVertices, ( tVertexPtrs.data() ).data() );

        tVertexPtrs.insert( 0, tVertices, tVertices + tNumVertices );

        return tVertexPtrs;
    }

    //------------------------------------------------------------------------------

    // TODO MESH-CLEANUP
    void
    Cell_DataBase::remove_vertex_pointer( moris_index aIndex )
    {
        std::cout << "In Cell Database" << std::endl;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cell_DataBase::get_vertex_coords() const
    {
        // get dimension of the matrix
        uint tNumVertices = this->get_number_of_vertices();
        uint tDim         = mMesh->get_spatial_dim();

        // set the dimension of the problem
        Matrix< DDRMat > tVertexCoords( tNumVertices, tDim );

        Vertex** tVertices = mMesh->get_cell_vertices( mCellIndex2 );

        // loop over vertices to set the coords of the individual vertices
        for ( uint i = 0; i < tNumVertices; i++ )
        {
            Matrix< DDRMat > tVertCoord = tVertices[ i ]->get_coords();

            tVertexCoords.set_row( i, tVertCoord );
        }

        // return the output
        return tVertexCoords;
    }

    //------------------------------------------------------------------------------
    uint
    Cell_DataBase::get_level() const
    {
        return mBaseCell->get_level();
    }

    //------------------------------------------------------------------------------

    mtk::Cell const *
    Cell_DataBase::get_base_cell() const
    {
        return mBaseCell;
    }

    //------------------------------------------------------------------------------

    mtk::Cell*
    Cell_DataBase::get_base_cell()
    {
        return mBaseCell;
    }

    //------------------------------------------------------------------------------

    const luint*
    Cell_DataBase::get_ijk() const
    {
        return mBaseCell->get_ijk();
    }

    //------------------------------------------------------------------------------

    size_t
    Cell_DataBase::capacity()
    {
        size_t tTotalSize = 0;

        tTotalSize += sizeof( mBaseCell );
        tTotalSize += sizeof( mCellIndex2 );
        tTotalSize += sizeof( mMesh );

        tTotalSize += sizeof( mCellId );
        tTotalSize += sizeof( mCellIndex );
        tTotalSize += sizeof( mCellOwner );
        tTotalSize += sizeof( mCellInfo );

        return tTotalSize;
    }

    //------------------------------------------------------------------------------

    Cell_Info const *
    Cell_DataBase::get_cell_info() const
    {
        return this->get_cell_info_sp().get();
    }
    //------------------------------------------------------------------------------

    std::shared_ptr< mtk::Cell_Info >
    Cell_DataBase::get_cell_info_sp() const
    {
        return mMesh->get_cell_info_sp( mCellIndex2 );
    }

    //------------------------------------------------------------------------------

    moris_id
    Cell_DataBase::get_id() const
    {
        return mMesh->get_entity_id( EntityRank::ELEMENT, mCellIndex2 );
    }

    //------------------------------------------------------------------------------

    moris_index
    Cell_DataBase::get_index() const
    {
        return mCellIndex2;
    }

    //------------------------------------------------------------------------------

    moris_id
    Cell_DataBase::get_owner() const
    {
        return mMesh->get_entity_owner( EntityRank::ELEMENT, mCellIndex2 );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::mtk
