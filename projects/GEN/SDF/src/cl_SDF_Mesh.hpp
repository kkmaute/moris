/*
 * cl_SDF_Mesh.hpp
 *
 *  Created on: Oct 3, 2018
 *      Author: messe
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_MESH_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_MESH_HPP_

#include <memory>
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Mesh.hpp"

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

//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            /**
             * constructor
             */
            Mesh( std::shared_ptr< mtk::Mesh > aMesh, bool aVerbose = false );

            Mesh( mtk::Mesh * aMesh , bool aVerbose = false );

//-------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Mesh();

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
        private:
//-------------------------------------------------------------------------------

            void
            link_vertices_with_neighbors();

//-------------------------------------------------------------------------------
        };
//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_MESH_HPP_ */
