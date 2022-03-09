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
        , mLocalSpgIndex( -1 )
        , mBulkPhaseIndex( aBulkPhaseIndex )
{
}

//------------------------------------------------------------------------------

Interpolation_Cell_Unzipped::Interpolation_Cell_Unzipped(
        moris::mtk::Cell*                        aBaseCell,
        moris_index                              aLocalSpgIndex,
        moris_id                                 aEnrCellId,
        moris_index                              aEnrCellIndex,
        moris_id                                 aEnrCellOwner,
        std::shared_ptr< moris::mtk::Cell_Info > aConnectivity,
        bool                                     aIsSpgBasedConstruction )
        : Interpolation_Cell( aEnrCellId, aEnrCellIndex, aEnrCellOwner, aConnectivity )
        , mBaseCell( aBaseCell )
        , mSubPhaseIndex( -1 )
        , mLocalSpgIndex( aLocalSpgIndex )
        , mBulkPhaseIndex( -1 )
{
    MORIS_ASSERT( aIsSpgBasedConstruction, 
        "Interpolation_Cell_Unzipped::Interpolation_Cell_Unzipped() - second constructor can only be used with SPG-based enrichment" );
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
    MORIS_ASSERT( mSubPhaseIndex != -1, "Interpolation_Cell_Unzipped::get_subphase_index() - UIPV has been constructed using SPGs, no SP index available." );
    return mSubPhaseIndex;
}

//------------------------------------------------------------------------------

moris_index
Interpolation_Cell_Unzipped::get_local_SPG_index() const
{
    MORIS_ASSERT( mLocalSpgIndex != -1, "Interpolation_Cell_Unzipped::get_local_SPG_index() - UIPV has been constructed using SPs, no SPG index available." );
    return mLocalSpgIndex;
}

//------------------------------------------------------------------------------

moris_index
Interpolation_Cell_Unzipped::get_bulkphase_index() const
{
    MORIS_ASSERT( mBulkPhaseIndex != -1, "Interpolation_Cell_Unzipped::get_bulkphase_index() - UIPV has been constructed using SPGs, bulk-phase index not well-defined." );
    return mBulkPhaseIndex;
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

