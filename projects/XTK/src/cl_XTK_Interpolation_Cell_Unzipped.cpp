/*
 * cl_XTK_Interpolation_Cell_Unzipped.cpp
 *
 *  Created on: Jul 10, 2019
 *      Author: doble
 */


#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
namespace xtk
{

//------------------------------------------------------------------------------

Interpolation_Cell_Unzipped::Interpolation_Cell_Unzipped(
        moris::mtk::Cell*                        aBaseCell,
        moris_index                              aSubphaseIndex,
        moris_index                              aBulkPhaseIndex,
        moris_id                                 aCellId,
        moris_index                              aCellIndex,
        moris_id                                 aCellOwner,
        std::shared_ptr< moris::mtk::Cell_Info > aConnectivity )
        : Interpolation_Cell( aCellId, aCellIndex, aCellOwner, aConnectivity )
        , mBaseCell( aBaseCell )
        , mSubPhaseIndex( aSubphaseIndex )
        , mBulkPhaseIndex( aBulkPhaseIndex )
{
}

//------------------------------------------------------------------------------

Interpolation_Cell_Unzipped::Interpolation_Cell_Unzipped(
        moris::mtk::Cell*                        aBaseCell,
        moris_index                              aPrimarySubPhaseIndex,
        moris_index                              aPrimaryBulkPhaseIndex,
        moris_index                              aEnrCellIndex,
        moris_id                                 aEnrCellOwner,
        moris_index                              aNumMeshIndices,
        std::shared_ptr< moris::mtk::Cell_Info > aConnectivity,
        bool                                     aIsSpgBasedConstruction )
        : Interpolation_Cell( MORIS_ID_MAX, aEnrCellIndex, aEnrCellOwner, aConnectivity )
        , mBaseCell( aBaseCell )
        , mSubPhaseIndex( aPrimarySubPhaseIndex )
        , mBulkPhaseIndex( aPrimaryBulkPhaseIndex )
{
    MORIS_ASSERT( aIsSpgBasedConstruction, 
        "Interpolation_Cell_Unzipped::Interpolation_Cell_Unzipped() - second constructor can only be used with SPG-based enrichment" );

    // size lists
    mSpgIndices.resize( aNumMeshIndices, -1 );
    mBulkPhaseIndices.resize( aNumMeshIndices, -1 );
}

//------------------------------------------------------------------------------

uint
Interpolation_Cell_Unzipped::get_level() const
{
    return mBaseCell->get_level();
}

//------------------------------------------------------------------------------

uint
Interpolation_Cell_Unzipped::get_number_of_vertices() const
{
    return mVertices.size();
}

//------------------------------------------------------------------------------

moris::Cell< mtk::Vertex* >
Interpolation_Cell_Unzipped::get_vertex_pointers() const
{

    uint tNumVerts = this->get_number_of_vertices();
    moris::Cell< mtk::Vertex* > tVerts(tNumVerts);

    for(uint i = 0; i < tNumVerts; i++)
    {
        tVerts(i) = mVertices(i);
    }

    return tVerts;
}

//------------------------------------------------------------------------------

Matrix< DDRMat >
Interpolation_Cell_Unzipped::get_vertex_coords() const
{
    return mBaseCell->get_vertex_coords();
}

//------------------------------------------------------------------------------

void
Interpolation_Cell_Unzipped::set_vertices(moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & aVertexPointers)
{
    mVertices = aVertexPointers;
}

//------------------------------------------------------------------------------

moris::mtk::Cell const*
Interpolation_Cell_Unzipped::get_base_cell() const
{
    return mBaseCell;
}

//------------------------------------------------------------------------------

moris::mtk::Cell*
Interpolation_Cell_Unzipped::get_base_cell()
{
    return mBaseCell;
}

//------------------------------------------------------------------------------

moris_index
Interpolation_Cell_Unzipped::get_subphase_index() const
{
    return mSubPhaseIndex;
}

//------------------------------------------------------------------------------

moris_index
Interpolation_Cell_Unzipped::get_bulkphase_index() const
{
    return mBulkPhaseIndex;
}

//------------------------------------------------------------------------------

void
Interpolation_Cell_Unzipped::set_SPG_and_BP_indices_for_DM_list_index(
        const moris_index aMeshListIndex,
        const moris_index aSpgIndex,
        const moris_index aBulkPhaseIndex )
{
    mSpgIndices( aMeshListIndex ) = aSpgIndex;
    mBulkPhaseIndices( aMeshListIndex ) = aBulkPhaseIndex;
}

//------------------------------------------------------------------------------

moris_index
Interpolation_Cell_Unzipped::get_SPG_index_for_DM_list_index( const moris_index aMeshListIndex ) const
{
    return mSpgIndices( aMeshListIndex );
}

//------------------------------------------------------------------------------

moris_index
Interpolation_Cell_Unzipped::get_bulkphase_index_for_DM_list_index( const moris_index aMeshListIndex ) const
{
    return mBulkPhaseIndices( aMeshListIndex );
}

//------------------------------------------------------------------------------

moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const &
Interpolation_Cell_Unzipped::get_xtk_interpolation_vertices() const
{
    return mVertices;
}
//------------------------------------------------------------------------------

moris::Cell< xtk::Interpolation_Vertex_Unzipped* > &
Interpolation_Cell_Unzipped::get_xtk_interpolation_vertices()
{
    return mVertices;
}
//------------------------------------------------------------------------------
size_t
Interpolation_Cell_Unzipped::capacity()
{
    size_t tTotal = 0;
    tTotal += sizeof(mBaseCell);
    tTotal += sizeof(mSubPhaseIndex);
    tTotal += sizeof(mSubPhaseIndex);

    // FIXME: add for SPGs
    // tTotal += sizeof(mSubphaseGroupIndex);
    // tTotal += sizeof(mSubphaseGroupIndex);
    tTotal += sizeof(mBulkPhaseIndex);
    tTotal += mVertices.capacity();
    return tTotal;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

}

