/*
 * cl_XTK_Interpolation_Vertex_Unzipped.cpp
 *
 *  Created on: Jul 10, 2019
 *      Author: doble
 */


#include "cl_XTK_Interpolation_Vertex_Unzipped.hpp"

#include "cl_MTK_Vertex.hpp"

namespace xtk
{
//------------------------------------------------------------------------------

Interpolation_Vertex_Unzipped::Interpolation_Vertex_Unzipped(mtk::Vertex*               aBaseInterpVertex,
                                                             moris_id                   aVertexId,
                                                             moris_index                aVertexIndex,
                                                             moris_index                aVertexOwner,
                                                             uint                       aInterpolationOrder,
                                                             mtk::Vertex_Interpolation* aVertexInterp):
                                                                     mBaseInterpVertex(aBaseInterpVertex),
                                                                     mVertexId(aVertexId),
                                                                     mVertexIndex(aVertexIndex),
                                                                     mVertexOwner(aVertexOwner),
                                                                     mOrder(aInterpolationOrder),
                                                                     mInterpolation(aVertexInterp)
{

}
//------------------------------------------------------------------------------
Matrix< DDRMat >
Interpolation_Vertex_Unzipped::get_coords() const
{
    return mBaseInterpVertex->get_coords();
}
//------------------------------------------------------------------------------
moris_id
Interpolation_Vertex_Unzipped::get_id() const
{
    return mVertexId;
}
//------------------------------------------------------------------------------
moris_index
Interpolation_Vertex_Unzipped::get_index() const
{
    return mVertexIndex;
}
//------------------------------------------------------------------------------
moris_index
Interpolation_Vertex_Unzipped::get_owner() const
{
    return mVertexOwner;
}
//------------------------------------------------------------------------------
mtk::Vertex_Interpolation *
Interpolation_Vertex_Unzipped::get_interpolation( const uint aOrder )
{
    MORIS_ERROR(mOrder == aOrder,"Invalid interpolation order specified");
    return mInterpolation;
}
//------------------------------------------------------------------------------
const mtk::Vertex_Interpolation *
Interpolation_Vertex_Unzipped::get_interpolation( const uint aOrder ) const
{
    MORIS_ERROR(mOrder == aOrder,"Invalid interpolation order specified");
    return mInterpolation;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

}
