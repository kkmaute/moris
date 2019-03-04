/*
 * cl_MTK_Mesh_XTK_Impl.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: doble
 */

#include "cl_MTK_Mesh_XTK_Impl.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enrichment.hpp"

namespace moris
{
namespace mtk
{
// ----------------------------------------------------------------------------------
// Constructor/Deconstructor Source code
// ----------------------------------------------------------------------------------

XTK_Impl::XTK_Impl(xtk::Model* aModelPtr):
    mXTKModelPtr(aModelPtr),
    mOutputMeshPtr(aModelPtr->get_output_mesh()){}

// ----------------------------------------------------------------------------------
// access connectivity
// ----------------------------------------------------------------------------------

Matrix<IndexMat>
XTK_Impl::get_entity_connected_to_entity_loc_inds(moris_index     aEntityIndex,
                                                  enum EntityRank aInputEntityRank,
                                                  enum EntityRank aOutputEntityRank) const
{
    if(aInputEntityRank == EntityRank::ELEMENT  && aOutputEntityRank == EntityRank::NODE)
    {
        xtk::Background_Mesh & tBM = mXTKModelPtr->get_background_mesh();

        return tBM.get_mtk_cell(aEntityIndex).get_vertex_inds();
    }

    else
    {
        MORIS_ASSERT(0," Specified connectivity not implemented: aInputEntityRank = %u and aOutputEntityRank = %u",(moris::uint)aInputEntityRank,(moris::uint)aOutputEntityRank);
        return Matrix<IndexMat>(0,0);
    }
}

// ----------------------------------------------------------------------------------
// T matrix things
// ----------------------------------------------------------------------------------
const Matrix< DDRMat > &
XTK_Impl::get_t_matrix_of_node_loc_ind(
        const moris_index aNodeIndex,
        const EntityRank  aBSplineRank )
{
    xtk::Vertex_Enrichment const & tVertexEnrichment = mXTKModelPtr->get_basis_enrichment().get_vertex_enrichment(aNodeIndex);
    return tVertexEnrichment.get_basis_weights();
}

//------------------------------------------------------------------------------

Matrix< IndexMat >
XTK_Impl::get_bspline_inds_of_node_loc_ind(
        const moris_index aNodeIndex,
        const EntityRank  aBSplineRank )
{
    xtk::Vertex_Enrichment const & tVertexEnrichment = mXTKModelPtr->get_basis_enrichment().get_vertex_enrichment(aNodeIndex);
    return tVertexEnrichment.get_basis_basis_indices();
}

// ----------------------------------------------------------------------------------
// Index and ID allocation
// ----------------------------------------------------------------------------------
moris_id
XTK_Impl::get_max_entity_id( enum EntityRank aEntityRank ) const
{
    MORIS_ERROR(0,"Entered  function in Mesh base class, (function is not implemented)");
    return mXTKModelPtr->get_background_mesh().allocate_entity_ids(1,aEntityRank);
}
}
}

