/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_Interpolation_XTK_Impl.cpp
 *
 */

#include "cl_MTK_Vertex_Interpolation_XTK_Impl.hpp"

#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Vertex_Enrichment.hpp"
namespace moris
{
    namespace mtk
    {
        void
        Vertex_Interpolation_XTK::set_vertex_enrichment( xtk::Vertex_Enrichment *aVertexEnrichment )
        {
            MORIS_ASSERT( aVertexEnrichment != nullptr, "Null pointer provided to XTK vertex interpolation" );
            mVertexEnrichment = aVertexEnrichment;
        }

        Matrix< IndexMat >
        Vertex_Interpolation_XTK::get_indices() const
        {
            MORIS_ERROR( mVertexEnrichment != nullptr, "mVertexEnrichment not set in XTK vertex interpolation" );
            return mVertexEnrichment->get_basis_indices();
        }

        const Matrix< DDRMat > *
        Vertex_Interpolation_XTK::get_weights() const
        {
            MORIS_ERROR( mVertexEnrichment != nullptr, "mVertexEnrichment not set in XTK vertex interpolation" );
            return &mVertexEnrichment->get_basis_weights();
        }

    }    // namespace mtk
}    // namespace moris
