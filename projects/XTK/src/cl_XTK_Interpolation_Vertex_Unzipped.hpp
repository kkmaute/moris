/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Interpolation_Vertex_Unzipped.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_VERTEX_UNZIPPED_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_VERTEX_UNZIPPED_HPP_

#include "moris_typedefs.hpp"
#include "cl_MTK_Vertex.hpp"

#include "cl_XTK_Vertex_Enrichment.hpp"

namespace moris::xtk
{
    class Vertex_Enrichment;
}

using namespace moris;

namespace moris::xtk
{
    class Interpolation_Vertex_Unzipped : public mtk::Vertex
    {
      private:
        // the parent vertex
        mtk::Vertex* mBaseInterpVertex;

        moris_id    mVertexId;
        moris_index mVertexIndex;
        moris_index mVertexOwner;

        // interpolation of vertex
        Vector< Vertex_Enrichment* > mInterpolation;    // Vertex Enrichment owned by ip mesh

      public:
        Interpolation_Vertex_Unzipped(){};

        Interpolation_Vertex_Unzipped(
                mtk::Vertex*       aBaseInterpVertex,
                moris_id           aVertexId,
                moris_index        aVertexIndex,
                moris_index        aVertexOwner,
                uint               aInterpolationOrder,
                Vertex_Enrichment* aVertexInterp,
                uint               aMaxIpOrder );

        //------------------------------------------------------------------------------

        Matrix< DDRMat >                 get_coords() const;
        moris_id                         get_id() const;
        moris_index                      get_index() const;
        mtk::Vertex const *              get_base_vertex() const;
        moris_index                      get_owner() const;
        mtk::Vertex_Interpolation*       get_interpolation( const uint aOrder );
        const mtk::Vertex_Interpolation* get_interpolation( const uint aOrder ) const;
        uint                             get_num_vertex_interpolations() const;
        bool                             has_interpolation( const uint aBSplineMeshIndex );

        //------------------------------------------------------------------------------

        Vertex_Enrichment*        get_xtk_interpolation( const uint aOrder );
        Vertex_Enrichment const * get_xtk_interpolation( const uint aOrder ) const;

        //------------------------------------------------------------------------------

        void
        add_vertex_interpolation(
                const uint         aOrder,
                Vertex_Enrichment* aVertexInterp );

        //------------------------------------------------------------------------------

        void
        print() const;

        //------------------------------------------------------------------------------

        // memory
        size_t
        capacity();

        friend class Ghost_Stabilization;
        friend class Enriched_Interpolation_Mesh;

        mtk::Vertex* get_base_vertex();

      protected:
        void
        set_vertex_id( moris_index const & aId );
    };

    inline std::ostream&
    operator<<( std::ostream& os, const xtk::Interpolation_Vertex_Unzipped& dt )
    {
        os << "Vertex Id: " << dt.get_id()
           << " | Vertex Index: " << dt.get_index()
           << " | Base Vertex Id: " << dt.get_base_vertex()->get_id()
           << " | Vertex Interpolation: " << *( dt.get_xtk_interpolation( 0 ) );
        return os;
    }
}    // namespace moris::xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_VERTEX_UNZIPPED_HPP_ */
