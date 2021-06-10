/*
 * cl_VIS_Vertex_Visualization.hpp
 *
 *  Created on: Dec 121, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_VERTEX_VISUALIZATION_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_VERTEX_VISUALIZATION_STK_HPP_

#include "cl_MTK_Vertex.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"

namespace moris
{
    namespace vis
    {
//------------------------------------------------------------------------------

        class Vertex_Visualization: public mtk::Vertex
        {
        private:

            moris_id       mVertexId;
            moris_index    mVertexInd;
            mtk::Vertex *  mIntergrationVertex = nullptr;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             *  constructor
             */
            Vertex_Visualization( moris_id      aVertexId,
                                  moris_index   aVertexInd,
                                  mtk::Vertex * aIntegrationVertex ) : mVertexId( aVertexId ),
                                                                        mVertexInd( aVertexInd ),
                                                                        mIntergrationVertex( aIntegrationVertex )

            {};

//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Vertex_Visualization(){};

//------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Vertex_Visualization(){};

//------------------------------------------------------------------------------

            /**
             * returns a moris::Matrix with node coordinates
             */
            Matrix< DDRMat > get_coords() const
            {
                MORIS_ASSERT( mIntergrationVertex!=nullptr, "get_coords(), Integration vertex is nullptr" );

                return mIntergrationVertex->get_coords();
            }

//------------------------------------------------------------------------------

            /**
             * returns the vis mesh domain wide id of this vertex
             */
            moris_id get_id() const
            {
                return mVertexId;
            }


//------------------------------------------------------------------------------

            /**
             * returns the vis mesh domain wide id of this vertex
             */
            moris_index get_index() const
            {
                return mVertexInd;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the id used in the integration mesh
             */
            moris_id get_integration_id() const
            {
                return mIntergrationVertex->get_id();
            }


            //------------------------------------------------------------------------------

            /**
             * returns the index used in the integration mesh
             */
            moris_index get_integration_index() const
            {
                return mIntergrationVertex->get_index();
            }

//------------------------------------------------------------------------------

            /**
             * returns the id of the proc that owns this vertex
             */
            moris_id get_owner() const
            {
                return mIntergrationVertex->get_owner();
            }

//------------------------------------------------------------------------------

            mtk::Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex )
            {
                MORIS_ERROR( false,"get_interpolation(), not implemented for visualization vertex");
                return nullptr;
            }

            const mtk::Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex ) const
            {
                MORIS_ERROR( false,"get_interpolation(), not implemented for visualization vertex");
                return nullptr;
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace vis */
} /* namespace moris */
    //------------------------------------------------------------------------------



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_VERTEX_VISUALIZATION_STK_HPP_ */
