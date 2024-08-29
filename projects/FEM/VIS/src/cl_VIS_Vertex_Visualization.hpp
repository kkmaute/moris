/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Vertex_Visualization.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_VERTEX_VISUALIZATION_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_VERTEX_VISUALIZATION_STK_HPP_

#include "cl_MTK_Vertex.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris::vis
{
    //------------------------------------------------------------------------------

    class Vertex_Visualization : public mtk::Vertex
    {
        //------------------------------------------------------------------------------

      private:
        moris_id     mVertexId;
        moris_index  mVertexInd;
        mtk::Vertex *mIntegrationVertex = nullptr;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        /**
         *  constructor
         */
        Vertex_Visualization( moris_id aVertexId,
                moris_index            aVertexInd,
                mtk::Vertex           *aIntegrationVertex )
                : mVertexId( aVertexId )
                , mVertexInd( aVertexInd )
                , mIntegrationVertex( aIntegrationVertex ){
                    // do nothing else
                };

        //------------------------------------------------------------------------------

        /**
         * trivial constructor
         */
        Vertex_Visualization(){};

        //------------------------------------------------------------------------------

        /**
         * destructor
         */
        ~Vertex_Visualization() override{};

        //------------------------------------------------------------------------------

        /**
         * returns a moris::Matrix with node coordinates
         */
        Matrix< DDRMat >
        get_coords() const override
        {
            MORIS_ASSERT( mIntegrationVertex != nullptr, "get_coords(), Integration vertex is nullptr" );

            return mIntegrationVertex->get_coords();
        }

        //------------------------------------------------------------------------------

        /**
         * returns the vis mesh domain wide id of this vertex
         */
        moris_id
        get_id() const override
        {
            return mVertexId;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the vis mesh domain wide id of this vertex
         */
        moris_index
        get_index() const override
        {
            return mVertexInd;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the id used in the integration mesh
         */
        moris_id
        get_integration_id() const
        {
            return mIntegrationVertex->get_id();
        }

        //------------------------------------------------------------------------------

        /**
         * returns the index used in the integration mesh
         */
        moris_index
        get_integration_index() const
        {
            return mIntegrationVertex->get_index();
        }

        //------------------------------------------------------------------------------

        /**
         * @brief returns
         *
         * @return mtk::Vertex const*
         */
        mtk::Vertex const *
        get_base_vertex() const override
        {
            return mIntegrationVertex;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the id of the proc that owns this vertex
         */
        moris_id
        get_owner() const override
        {
            return mIntegrationVertex->get_owner();
        }

        //------------------------------------------------------------------------------

        mtk::Vertex_Interpolation *
        get_interpolation( const uint aBSplineMeshIndex ) override
        {
            MORIS_ERROR( false, "get_interpolation(), not implemented for visualization vertex" );
            return nullptr;
        }

        //------------------------------------------------------------------------------

        const mtk::Vertex_Interpolation *
        get_interpolation( const uint aBSplineMeshIndex ) const override
        {
            MORIS_ERROR( false, "get_interpolation(), not implemented for visualization vertex" );
            return nullptr;
        }

        //------------------------------------------------------------------------------

    };    // class vis::Vertex_Visualization

    //------------------------------------------------------------------------------

}    // namespace moris::vis

//------------------------------------------------------------------------------

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_VERTEX_VISUALIZATION_STK_HPP_ */
