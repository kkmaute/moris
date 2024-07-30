/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell.hpp
 *
 */

#ifndef SRC_MESH_CL_MTK_CELL_HPP_
#define SRC_MESH_CL_MTK_CELL_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src
#include "cl_Matrix.hpp"
#include "fn_isrow.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Vertex.hpp"    //MTK/src
#include "cl_MTK_Enums.hpp"     //MTK/src

//------------------------------------------------------------------------------
// Forward Declarations
namespace moris
{
    namespace mtk
    {
        class Cell_Info;
    }
}    // namespace moris
//------------------------------------------------------------------------------

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------
        /**
         * \brief the mtk::Cell class provides the cell information that is
         * provided by the mesh.
         */

        class Cell
        {
          protected:    // protected so the derived classes can set
            std::shared_ptr< moris::mtk::Cell_Info > mCellInfo;

            moris_id mCellId;
            moris_id mCellIndex;
            moris_id mCellOwner;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Cell(){};

            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Cell(
                    moris_id                                 aCellId,
                    moris_index                              aCellIndex,
                    moris_id                                 aCellOwner,
                    std::shared_ptr< moris::mtk::Cell_Info > aCellInfo );

            //------------------------------------------------------------------------------

            /**
             * Destructor
             */
            virtual ~Cell(){};

            //------------------------------------------------------------------------------
            virtual mtk::Cell_Info const *
            get_cell_info() const;

            //------------------------------------------------------------------------------

            virtual std::shared_ptr< mtk::Cell_Info >
            get_cell_info_sp() const;

            //------------------------------------------------------------------------------

            void
            set_mtk_cell_info( std::shared_ptr< moris::mtk::Cell_Info > aCellInfo );

            //------------------------------------------------------------------------------

            /**
             * returns the domain wide id of the cell
             *
             * @return moris_id ID
             */
            virtual moris_id
            get_id() const;

            //------------------------------------------------------------------------------
            /**
             * set the cell id
             */
            virtual void
            set_id( moris_id aId );

            //------------------------------------------------------------------------------
            /**
             * set the cell id
             */
            virtual void
            set_owner( moris_id aOwner );

            /**
             * set the cell id
             */
            virtual void
            set_index( moris_id aIndex );

            //------------------------------------------------------------------------------

            /**
             * returns the local index of the cell
             *
             * @return moris_index ID
             */
            virtual moris_index
            get_index() const;

            //------------------------------------------------------------------------------

            /**
             * returns the proc id of the owner of this cell
             * ( this information is needed for STK )
             */
            virtual moris_id
            get_owner() const;

            //------------------------------------------------------------------------------

            /*!
             * Returns the level that this cell is on. For most meshes this returns 0. However,
             * for HMR this is not trivial
             */
            virtual uint
            get_level() const;

            //------------------------------------------------------------------------------

            /**
             * @return Ptrs of vertices connected to this cell
             */
            virtual Vector< Vertex * >
            get_vertex_pointers() const = 0;

            //-----------------------------------------------------------------------------

            /**
             * @brief
             *
             * @param aIndex
             */

            virtual void
            remove_vertex_pointer( moris_index aIndex ) = 0;

            //-----------------------------------------------------------------------------

            virtual void
            remove_vertex( moris_index aIndex );

            //------------------------------------------------------------------------------

            /**
             *  @return number of vertices on this cell
             */
            virtual uint
            get_number_of_vertices() const;

            //------------------------------------------------------------------------------
            /**
             *  @return number of facets on this cell
             */
            virtual uint
            get_number_of_facets() const;

            //------------------------------------------------------------------------------
            /**
             *  @return number of edges on this cell
             */
            virtual uint
            get_number_of_edges() const;

            //------------------------------------------------------------------------------

            /**
             * returns a Mat with IDs of connected vertices
             */
            virtual Matrix< IdMat >
            get_vertex_ids() const;

            //------------------------------------------------------------------------------

            /**
             * returns a Mat with indices of connected vertices
             */
            virtual Matrix< IndexMat >
            get_vertex_inds() const;

            //------------------------------------------------------------------------------

            // TODO MESHCLEANUP
            // virtual void
            // set_vertex_pointers( mtk::Vertex* aVertex, moris_index aIndex );

            //------------------------------------------------------------------------------

            virtual bool
            check_unique_vertex_inds() const;

            //------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_vertex_owners() const;

            //------------------------------------------------------------------------------

            /**
             * returns a Mat of dimension
             * < number of vertices * number of dimensions >
             */
            virtual Matrix< DDRMat >
            get_vertex_coords() const = 0;

            //------------------------------------------------------------------------------

            virtual Vector< mtk::Vertex_Interpolation * >
            get_vertex_interpolations( const uint aOrder ) const;

            //------------------------------------------------------------------------------

            /*!
             * get vertices on side ordinal.
             * This functions is needed for side clustering
             */
            virtual Vector< moris::mtk::Vertex const * >
            get_vertices_on_side_ordinal( moris::moris_index aSideOrdinal ) const;

            //------------------------------------------------------------------------------

            /*!
             * get vertices on side ordinal that define the geometry (i.e. the corner nodes)
             * This functions is needed for side clustering
             */
            virtual Vector< moris::mtk::Vertex const * >
            get_geometric_vertices_on_side_ordinal( moris::moris_index aSideOrdinal ) const;

            //------------------------------------------------------------------------------

            /*!
             * Get vertex coordinates on side ordinal
             */

            virtual moris::Matrix< moris::DDRMat >
            get_cell_physical_coords_on_side_ordinal( moris::moris_index aSideOrdinal ) const;

            //------------------------------------------------------------------------------

            /*!
             * get vertices on side ordinal.
             * This functions is needed for side clustering
             */
            moris::Matrix< IndexMat >
            get_vertices_ind_on_side_ordinal( moris::moris_index aSideOrdinal ) const;

            //------------------------------------------------------------------------------

            /*!
             * @param[in] aVertexIndex - Proc local vertex index
             * @return the vertex ordinal wrt the mtk cell
             */
            virtual moris_index
            get_vertex_ordinal_wrt_cell( moris_index const &aVertexIndex ) const;

            //------------------------------------------------------------------------------

            /**
             * returns an enum that defines the geometry type of the element
             */
            virtual Geometry_Type
            get_geometry_type() const;

            //------------------------------------------------------------------------------

            /*!
             * Compute facet normal
             */
            virtual moris::Matrix< moris::DDRMat >
            compute_outward_side_normal( moris::moris_index aSideOrdinal ) const;

            //------------------------------------------------------------------------------

            /*
             * Volume in 3D, Surface Area in 2D
             */
            virtual moris::real
            compute_cell_measure() const;

            //------------------------------------------------------------------------------

            /*
             * Derivative of the Volume in 3D, Surface Area in 2D wrt to a single dof
             * @param[in] aLocalVertexID  Local ID of vertex to use (0, 1, 2, 3, etc).
             * @param[in] aDirection      Direction to take derivative (0,1, or 2).
             */
            virtual moris::real
            compute_cell_measure_deriv(
                    uint aLocalVertexID,
                    uint aDirection ) const;

            //------------------------------------------------------------------------------

            /*
             * Surface Area on side of cell in 3D, line length on side in 2D
             */
            virtual moris::real
            compute_cell_side_measure( moris_index const &aCellSideOrd ) const;

            //------------------------------------------------------------------------------

            /*
             * Derivative of the surface Area on side of cell in 3D, line length on side in 2D
             * @param[in] aCellSideOrd    Side Ordinal for the considered cell side
             * @param[in] aLocalVertexID  Local ID of vertex to use (0, 1, 2, 3, etc).
             * @param[in] aDirection      Direction to take derivative (0,1, or 2).
             */
            virtual moris::real
            compute_cell_side_measure_deriv(
                    moris_index const &aCellSideOrd,
                    uint               aLocalVertexID,
                    uint               aDirection ) const;

            //------------------------------------------------------------------------------

            /*
             * Cell centroid
             */
            virtual moris::Matrix< moris::DDRMat >
            compute_cell_centroid() const;

            //------------------------------------------------------------------------------

            /**
             * returns the order of the element
             */
            virtual Interpolation_Order
            get_interpolation_order() const;

            //------------------------------------------------------------------------------

            /**
             * returns the integration order of the element
             */
            virtual Integration_Order
            get_integration_order() const;

            //------------------------------------------------------------------------------

            virtual moris::mtk::Cell const *
            get_base_cell() const
            {
                return this;
            }

            //------------------------------------------------------------------------------

            virtual moris::mtk::Cell *
            get_base_cell()
            {
                return this;
            }

            //------------------------------------------------------------------------------

            virtual const luint *
            get_ijk() const
            {
                MORIS_ASSERT( false, "get_ijk( ) not implemented for base class" );
                return nullptr;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the order of the element
             */
            virtual size_t
            capacity()
            {
                MORIS_ERROR( 0, "No default implementation" );
                return 0;
            };

            //------------------------------------------------------------------------------

            /*!
             * Get geometric vertex coordinates on side ordinal
             * @param[ in ] aSideOrdinal Side ordinal of the cell
             */

            virtual moris::Matrix< moris::DDRMat >
            get_cell_geometric_coords_on_side_ordinal( moris::moris_index aSideOrdinal ) const;

            //------------------------------------------------------------------------------

        };    // class mtk::Cell

        //------------------------------------------------------------------------------

        // operators for printing
        inline std::ostream &
        operator<<( std::ostream &os, const mtk::Cell &dt )
        {
            os << "Cell Id: " << dt.get_id() << " | Cell Index: " << dt.get_index();

            return os;
        }

        //------------------------------------------------------------------------------

        inline std::ostream &
        operator<<( std::ostream &os, mtk::Cell const *const &dt )
        {
            os << "Cell Id: " << dt->get_id() << " | Cell Index: " << dt->get_index();

            return os;
        }

        //------------------------------------------------------------------------------

    } /* namespace mtk */
} /* namespace moris */

//------------------------------------------------------------------------------

#endif /* SRC_MESH_CL_MTK_CELL_HPP_ */
