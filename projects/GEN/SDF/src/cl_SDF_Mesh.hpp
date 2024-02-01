/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Mesh.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_MESH_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_MESH_HPP_

#include <memory>
#include "moris_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Mesh_Core.hpp"

#include "cl_SDF_Vertex.hpp"
#include "cl_SDF_Cell.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        /**
         * Wrapper around an MTK mesh
         */
        class Mesh
        {
            //! pointer to underlying mesh
            mtk::Mesh * mMesh;

            //! vector with SDF Vertices
            moris::Cell< Vertex * > mVertices;

            //! vector with SDF Cells
            moris::Cell< Cell * > mCells;

            bool mVerbose;

            Matrix< F31RMat > mMinCoord;
            Matrix< F31RMat > mMaxCoord;

            Matrix< IdMat > mNodeIDs;

            // interpolation order
            uint mOrder;

//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            /**
             * constructor
             */
            Mesh( std::shared_ptr< mtk::Mesh > aMesh, bool aVerbose = false );

//-------------------------------------------------------------------------------

            Mesh( mtk::Mesh * aMesh , bool aVerbose = false );

//-------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Mesh();

//-------------------------------------------------------------------------------

            /**
             * expose mesh pointer
             */
            mtk::Mesh *
            get_mtk_mesh()
            {
                return mMesh;
            }

//-------------------------------------------------------------------------------
            uint
            get_num_nodes() const
            {
                return mMesh->get_num_nodes();
            }

//-------------------------------------------------------------------------------

            uint
            get_num_elems() const
            {
                return mMesh->get_num_elems();
            }

//-------------------------------------------------------------------------------

            const Matrix< F31RMat > &
            get_node_coordinate( const moris_index aIndex ) const
            {
                return mVertices( aIndex )->get_coords();
            }

//-------------------------------------------------------------------------------

            Vertex *
            get_vertex( const uint & aIndex )
            {
                return mVertices( aIndex );
            }

//-------------------------------------------------------------------------------

            Cell *
            get_cell( const uint & aIndex )
            {
                return mCells( aIndex );
            }

//-------------------------------------------------------------------------------

            real
            get_min_coord( const uint aIndex ) const
            {
                return mMinCoord( aIndex );
            }

//-------------------------------------------------------------------------------

            real
            get_max_coord( const uint aIndex ) const
            {
                return mMaxCoord( aIndex );
            }

//-------------------------------------------------------------------------------

            bool
            is_verbose() const
            {
                return mVerbose;
            }

//-------------------------------------------------------------------------------

            const Matrix< IdMat > &
            get_node_ids() const
            {
                return mNodeIDs;
            }

//-------------------------------------------------------------------------------

            /**
             * return the interpolation order of the mesh.
             * Needer for HDF5 output.
             * ( taken from first element on mesh, assuming that all elements
             *   are of the same order )
             */
            uint
            get_order() const
            {
                return mOrder;
            }

//-------------------------------------------------------------------------------
        private:
//-------------------------------------------------------------------------------

            void
            link_vertex_cells();

//-------------------------------------------------------------------------------
            void
            link_vertex_neighbors();

        };
//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_MESH_HPP_ */

