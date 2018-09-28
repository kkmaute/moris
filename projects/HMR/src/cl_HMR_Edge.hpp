/*
 * cl_Lagrange_Edge.hpp
 *
 *  Created on: Sep 26, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_EDGE_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_EDGE_HPP_

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Edge.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_HMR_Mesh_Base.hpp"
#include "cl_HMR_Element.hpp"
namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------

        class Edge : public mtk::Edge
        {
            //! id of this edge
            moris_id    mID;

            //! index of this edge
            moris_index mIndex;

            //! pointer with indices in elements
            Matrix< DDUMat > mIndicesInElements;

// ----------------------------------------------------------------------------
        protected:
// ----------------------------------------------------------------------------

            //! index on master
            uint      mIndexOfMaster;

            //! pointer with elements
            moris::Cell< Element * > mElements;

// ----------------------------------------------------------------------------
        public:
// ----------------------------------------------------------------------------

            /**
             * constructor
             */
            Edge( Mesh_Base       * aMesh,
                  Background_Edge * aBackgroundEdge );

// ----------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            virtual ~Edge(){};

// ----------------------------------------------------------------------------

            /**
             * returns the domain wide id of the edge
             *
             * @return moris_id ID
             */
            moris_id
            get_id() const;

// ----------------------------------------------------------------------------

            /**
             * returns the local index of the edge
             *
             * @return moris_index ID
             */
            moris_index
            get_index() const;

// ----------------------------------------------------------------------------

            /**
             * tells how many vertices are connected to this edge
             */
            virtual uint
            get_number_of_vertices() const = 0;

// ----------------------------------------------------------------------------

            /**
             * returns the proc id of the owner of this edge
             * ( this information is needed for STK )
             */
            moris_id
            get_owner() const;

// ----------------------------------------------------------------------------

            /**
             * fills a moris::cell with pointers to connected vertices
             */
            moris::Cell< mtk::Vertex* >
            get_vertex_pointers() const;

// ----------------------------------------------------------------------------

            /**
             * returns a Mat with IDs of connected vertices
             */
            Matrix< IdMat >
            get_vertex_ids() const;

// ----------------------------------------------------------------------------

            /**
             * returns a Mat with indices of connected vertices
             */
            virtual Matrix< IndexMat >
            get_vertex_inds() const;

// ----------------------------------------------------------------------------

            /**
             * returns a Mat of dimension
             * < number of vertices * number of dimensions >
             */
            Matrix< DDRMat >
            get_vertex_coords() const;

// ----------------------------------------------------------------------------

            /**
             * an edge is always a line
             */
            mtk::Geometry_Type
            get_geometry_type() const;

// ----------------------------------------------------------------------------

            virtual
            mtk::Interpolation_Order
            get_interpolation_order() const = 0;

// ----------------------------------------------------------------------------

            void
            set_index( const moris_index & aIndex );

// ----------------------------------------------------------------------------

            void
            set_id( const moris_id & aID );

// ----------------------------------------------------------------------------

            uint
            get_number_of_elements() const;

// ----------------------------------------------------------------------------

            Element *
            get_element( const uint & aIndex );

// ----------------------------------------------------------------------------

            uint
            get_index_on_element(  const uint & aIndex ) const;

// ----------------------------------------------------------------------------

            Element *
            get_hmr_master();

// ----------------------------------------------------------------------------

            uint
            get_index_on_master() const;

// ----------------------------------------------------------------------------

            virtual const Basis *
            get_basis( const uint aIndex ) const = 0;

// ----------------------------------------------------------------------------

            virtual Basis *
            get_basis( const uint aIndex )  = 0;

// ----------------------------------------------------------------------------

            bool
            is_active() const;

// ----------------------------------------------------------------------------
        private:
// ----------------------------------------------------------------------------

            void
            find_master(
                    Mesh_Base       * aMesh,
                    Background_Edge * aBackgroundEdge );

// ----------------------------------------------------------------------------

            virtual void
            copy_vertex_pointers() = 0;

// ----------------------------------------------------------------------------
        };

// ----------------------------------------------------------------------------
    }
}

#endif /* PROJECTS_HMR_SRC_CL_HMR_EDGE_HPP_ */
