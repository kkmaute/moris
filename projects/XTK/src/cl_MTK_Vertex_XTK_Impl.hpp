/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_XTK_Impl.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_MTK_VERTEX_XTK_IMPL_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_MTK_VERTEX_XTK_IMPL_HPP_

#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Vertex_Interpolation_XTK_Impl.hpp"

//------------------------------------------------------------------------------

namespace moris::mtk
{
    //------------------------------------------------------------------------------
    class Vertex_XTK : public Vertex
    {
        //------------------------------------------------------------------------------

      private:
        moris::moris_id                     mVertexId;
        moris::moris_index                  mVertexIndex;
        moris::moris_index                  mVertexOwner;
        Vertex_Interpolation_XTK*           mVertexInterpolation = nullptr;    // (basis weights and basis identity)
        std::shared_ptr< Matrix< DDRMat > > mCoordinates;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        /**
         * trivial constructor
         */
        Vertex_XTK();

        Vertex_XTK( moris::moris_id                 aVertexId,
                moris::moris_index                  aVertexIndex,
                moris::moris_index                  aOwner,
                std::shared_ptr< Matrix< DDRMat > > aCoordinates );
        //------------------------------------------------------------------------------

        /**
         * Destructor
         */
        ~Vertex_XTK() override;

        //------------------------------------------------------------------------------
        void
        set_vertex_interpolation( Vertex_Interpolation_XTK* aVertexInterpolation );

        //------------------------------------------------------------------------------

        /**
         * returns a Matrix with node coordinates
         */
        Matrix< DDRMat >
        get_coords() const override;

        //------------------------------------------------------------------------------

        /**
         * returns the domain wide id of this vertex
         */
        moris_id
        get_id() const override;

        //------------------------------------------------------------------------------

        /**
         * returns the processor unique index of this vertex
         */
        moris_index
        get_index() const override;

        //------------------------------------------------------------------------------

        void
        set_index( const moris_index aIndex );

        //------------------------------------------------------------------------------

        // fixme: change this into moris_id
        // FIXME: add owner function in background mesh
        moris_index
        get_owner() const override;

        //------------------------------------------------------------------------------

        Vertex_Interpolation*
        get_interpolation( const uint aOrder ) override;

        //------------------------------------------------------------------------------

        const Vertex_Interpolation*
        get_interpolation( const uint aOrder ) const override;

        //------------------------------------------------------------------------------

        uint
        get_level() const override;

        //------------------------------------------------------------------------------

        size_t
        capacity();

        //------------------------------------------------------------------------------

    };    // class Vertex_XTK

    //------------------------------------------------------------------------------

}    // namespace moris::mtk

//------------------------------------------------------------------------------

#endif /* PROJECTS_XTK_SRC_XTK_CL_MTK_VERTEX_XTK_IMPL_HPP_ */
