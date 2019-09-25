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
Interpolation_Cell_Unzipped::Interpolation_Cell_Unzipped(moris::mtk::Cell*         aBaseCell,
                                                         moris_index               aSubphaseIndex,
                                                         moris_index               aBulkPhaseIndex,
                                                         moris_id                  aCellId,
                                                         moris_index               aCellIndex,
                                                         moris_id                  aCellOwner,
                                                         moris::mtk::Cell_Info* aConnectivity):
                Interpolation_Cell(aCellId,aCellIndex,aCellOwner,aConnectivity),
                mBaseCell(aBaseCell),
                mSubPhaseIndex(aSubphaseIndex),
                mBulkPhaseIndex(aBulkPhaseIndex)
{
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
moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const &
Interpolation_Cell_Unzipped::get_xtk_interpolation_vertices() const
{
    return mVertices;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

}

