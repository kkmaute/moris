/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Node.hpp
 *
 */

#ifndef PROJECTS_FEM_SRC_CL_FEM_NODE_HPP_
#define PROJECTS_FEM_SRC_CL_FEM_NODE_HPP_

#include "typedefs.hpp"            //MRS/COR/src
#include "cl_MTK_Vertex.hpp"       //MTK/src
#include "cl_FEM_Node_Base.hpp"    //MTK/src

#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class Node : public Node_Base
        {
            //------------------------------------------------------------------------------

          private:
            // mtk vertex pointer
            const mtk::Vertex *mVertex;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * constructor
             */
            Node( const mtk::Vertex *aVertex )
                    : mVertex( aVertex )
            {
                //                mOwner = mVertex->get_owner();
            }

            //------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Node(){};

            //------------------------------------------------------------------------------

            /**
             * returns the owner of this node
             */
            auto
            get_owner() const -> decltype( mVertex->get_owner() )
            {
                return mVertex->get_owner();
            }

            //------------------------------------------------------------------------------

            /**
             * returns the T-Matrix of this node
             */
            const Matrix< DDRMat > *
            get_t_matrix( const uint aBSplineMeshIndex ) const
            {
                return mVertex->get_interpolation( aBSplineMeshIndex )->get_weights();
            }

            //------------------------------------------------------------------------------

            /**
             * returns the B-Spline IDs of this node
             */
            Matrix< IdMat >
            get_adof_ids( const uint aBSplineMeshIndex ) const
            {
                return mVertex->get_interpolation( aBSplineMeshIndex )->get_ids();
            }

            //------------------------------------------------------------------------------

            /**
             * returns the B-Spline IDs of this node
             */
            Matrix< IndexMat >
            get_adof_indices( const uint aBSplineMeshIndex ) const
            {
                return mVertex->get_interpolation( aBSplineMeshIndex )->get_indices();
            }

            //------------------------------------------------------------------------------

            /**
             * returns the proc owners of the IDs of this node
             */
            Matrix< IdMat >
            get_adof_owners( const uint aBSplineMeshIndex ) const
            {
                return mVertex->get_interpolation( aBSplineMeshIndex )->get_owners();
            }

            //------------------------------------------------------------------------------

            /**
             * get the ID of this node
             * @param[ in ] id id for this node
             */
            moris_id
            get_id() const
            {
                return mVertex->get_id();
            }

            //------------------------------------------------------------------------------

            /**
             * get the Index of this node
             * @param[ out ] index index for this node
             */
            moris_index
            get_index() const
            {
                return mVertex->get_index();
            }

            //------------------------------------------------------------------------------

            /**
             * is owned
             * @param[ out ] bool true if node is owned by processor
             */
            bool
            id_owned()
            {
                bool tOwned = true;

                if ( mVertex->get_owner() != par_rank() )
                {
                    tOwned = false;
                }

                return tOwned;
            }

            //------------------------------------------------------------------------------

            /**
             * get vertex coordinates
             * @param[ in ] aVertexCoords matrix to fill with vertex coordinates
             */
            void
            get_vertex_coords( Matrix< DDRMat > &aVertexCoords )
            {
                aVertexCoords = mVertex->get_coords();
            }

            //------------------------------------------------------------------------------

            const mtk::Vertex *
            get_geometric_vertex() const
            {
                return mVertex;
            }

            //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    }    // namespace fem
} /* namespace moris */

#endif /* PROJECTS_FEM_SRC_CL_NODE_HPP_ */
