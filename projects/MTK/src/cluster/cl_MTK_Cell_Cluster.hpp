/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Cluster.hpp
 *
 */
#ifndef PROJECTS_MTK_SRC_CL_MTK_CELL_CLUSTER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_CELL_CLUSTER_HPP_

#include "cl_Vector.hpp"
#include "cl_MTK_Cell.hpp"

#include "cl_MTK_Cluster.hpp"
#include "fn_norm.hpp"

namespace moris::mtk
{
    class Cluster_Group;

    class Cell_Cluster : public Cluster
    {
      public:
        Cell_Cluster() {};

        ~Cell_Cluster() override {};

        // ##############################################
        //  Characteristic functions
        // ##############################################

        bool
        is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override = 0;

        // ##############################################
        //  Cell/Vertex Access
        //  (Pure Virtual)
        // ##############################################

        //------------------------------------------------------------------------------

        Vector< moris::mtk::Cell const * > const &
        get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override = 0;

        //------------------------------------------------------------------------------

        Vector< moris::mtk::Cell const * > const &
        get_void_cells_in_cluster() const override = 0;

        //------------------------------------------------------------------------------

        moris::mtk::Cell const &
        get_interpolation_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override = 0;

        //------------------------------------------------------------------------------

        Vector< moris::mtk::Vertex const * >
        get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override = 0;

        //------------------------------------------------------------------------------

        void
        set_interpolation_cell( moris::mtk::Cell const *aInterpCell ) override
        {
            MORIS_ERROR( false, "set_interpolation_cell(), not implemented for this class" );
        };

        //------------------------------------------------------------------------------

        void
        add_primary_integration_cell( Vector< moris::mtk::Cell const * > const &aIntegrationCell ) override
        {
            MORIS_ERROR( false, "add_primary_integration_cell(), not implemented for this class" );
        };

        //------------------------------------------------------------------------------

        void
        add_void_integration_cell( Vector< moris::mtk::Cell const * > const &aIntegrationCell ) override
        {
            MORIS_ERROR( false, "add_primary_integration_cell(), not implemented for this class" );
        };

        //------------------------------------------------------------------------------

        void
        mark_as_nontrivial() override
        {
            MORIS_ERROR( false, "mark_as_nontrivial(), not implemented for this class" );
        };

        //------------------------------------------------------------------------------

        // ##############################################
        //  Local Coordinate Access
        //  (Pure Virtual)
        // ##############################################

        moris::Matrix< moris::DDRMat >
        get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override = 0;

        //------------------------------------------------------------------------------

        /*
         * Access a single local coordinate of a vertex
         */
        moris::Matrix< moris::DDRMat >
        get_vertex_local_coordinate_wrt_interp_cell(
                moris::mtk::Vertex const  *aVertex,
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override = 0;

        //------------------------------------------------------------------------------

        // ##############################################
        //  Size Access
        //  (Pure Virtual)
        // ##############################################

        //------------------------------------------------------------------------------

        /*!
         * Size of the xsi vector in this side cluster
         */
        moris_index
        get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override = 0;

        //------------------------------------------------------------------------------

        // ---------------------------------------------
        // EVERYTHING BELOW THIS LINE HAS A DEFAULT
        // IMPLEMENTATION
        // ---------------------------------------------

        // ##############################################
        //  Cell/Vertex Index Access
        // ##############################################

        //------------------------------------------------------------------------------

        moris::Matrix< moris::IndexMat >
        get_primary_cell_indices_in_cluster() const override
        {
            // number of cells in cluster
            moris::uint tNumCells = this->get_num_primary_cells();

            // cell access
            Vector< moris::mtk::Cell const * > const &tCells = this->get_primary_cells_in_cluster();

            // initialize output
            moris::Matrix< moris::IndexMat > tCellIndices( 1, tNumCells );

            // get cell indices and store
            for ( moris::uint i = 0; i < tNumCells; i++ )
            {
                tCellIndices( i ) = tCells( i )->get_index();
            }

            return tCellIndices;
        }

        //------------------------------------------------------------------------------

        moris::Matrix< moris::IndexMat >
        get_void_cell_indices_in_cluster() const override
        {
            MORIS_ERROR( !this->is_trivial(), "get_void_cell_indices_in_cluster on trivial cluster is not allowed" );

            // number of cells in cluster
            moris::uint tNumCells = this->get_num_void_cells();

            // cell access
            Vector< moris::mtk::Cell const * > const &tCells = this->get_void_cells_in_cluster();

            // initialize output
            moris::Matrix< moris::IndexMat > tCellIndices( 1, tNumCells );

            // get cell indices and store
            for ( moris::uint i = 0; i < tNumCells; i++ )
            {
                tCellIndices( i ) = tCells( i )->get_index();
            }

            return tCellIndices;
        }

        //----------------------------------------------------------------

        bool
        has_cluster_group( const moris_index aDiscretizationMeshIndex ) const override
        {
            MORIS_ERROR( false, "mtk::Cell_Cluster::has_cluster_group() - not implemented for this class" );
            return false;
        }

        //----------------------------------------------------------------

        std::shared_ptr< Cluster_Group >
        get_cluster_group( const moris_index aDiscretizationMeshIndex ) const override
        {
            MORIS_ERROR( false, "mtk::Cell_Cluster::get_cluster_group() - not implemented for this class" );
            return nullptr;
        }

        //----------------------------------------------------------------

        void
        set_cluster_group(
                const moris_index                            aDiscretizationMeshIndex,
                const std::shared_ptr< mtk::Cluster_Group > &aClusterGroupPtr ) override
        {
            MORIS_ERROR( false, "mtk::Cell_Cluster::set_cluster_group() - only implemented for child xtk::Cell_Cluster" );
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_cell_measure(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            moris::real tVolume = 0.0;

            Vector< moris::mtk::Cell const * > const *tCells = nullptr;

            if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY )
            {
                tCells = &this->get_primary_cells_in_cluster();
            }
            else
            {
                tCells = &this->get_void_cells_in_cluster();
            }

            for ( auto iC = tCells->cbegin(); iC < tCells->cend(); iC++ )
            {
                tVolume += ( *iC )->compute_cell_measure();
            }

            return tVolume;
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_measure(
                const moris_index          aDiscretizationMeshIndex,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            MORIS_ERROR( false, "mtk::Cell_Cluster::compute_cluster_group_cell_measure() - Only implemented in child xtk::Cell_Cluster" );
            return 0.0;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat >
        compute_cluster_ig_cell_measures(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            Vector< moris::mtk::Cell const * > const *tCells = nullptr;

            if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY )
            {
                tCells = &this->get_primary_cells_in_cluster();
            }
            else
            {
                tCells = &this->get_void_cells_in_cluster();
            }

            Matrix< DDRMat > tMeasureVec( tCells->size(), 1 );

            for ( uint iC = 0; iC < tCells->size(); iC++ )
            {
                tMeasureVec( iC ) = ( *tCells )( iC )->compute_cell_measure();
            }

            return tMeasureVec;
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_cell_measure_derivative(
                const Matrix< DDRMat >    &aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            moris::real tDerivative = 0.0;

            Vector< moris::mtk::Cell const * > const *tCells = nullptr;

            if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY )
            {
                tCells = &this->get_primary_cells_in_cluster();
            }
            else
            {
                tCells = &this->get_void_cells_in_cluster();
            }

            // loop over the ig cells in cluster
            for ( auto iC = tCells->cbegin(); iC < tCells->cend(); iC++ )
            {
                // get the cell coordinates
                Matrix< DDRMat > tCellCoords = ( *iC )->get_vertex_coords();

                // get number of nodes in this cell
                uint tNumNodesInCell = tCellCoords.n_rows();    // number of nodes in this cell

                // FIXME could be done with node index
                // check if this cell in cluster is affected by the perturbed node
                // flag true if cell is affected by perturbed cluster node
                bool tIsAffected = false;

                // init cell local node index
                uint tLocalVertexID = UINT_MAX;

                // loop over the nodes of the cell
                for ( uint iCellNode = 0; iCellNode < tNumNodesInCell; iCellNode++ )
                {
                    // check if perturbed cluster node affects this cell by using the distance between two nodes
                    tIsAffected = tIsAffected || ( moris::norm( aPerturbedVertexCoords - tCellCoords.get_row( iCellNode ) ) < 1e-12 );

                    // if the cell is affected by perturbed cluster node
                    if ( tIsAffected == true )
                    {
                        // get cell local node index
                        tLocalVertexID = iCellNode;

                        // break the loop as node correspondence was found
                        break;
                    }
                }

                // if the cell is affected by perturbed cluster node
                if ( tIsAffected == true )
                {
                    // add contribution from the cell to the cluster measure
                    tDerivative += ( *iC )->compute_cell_measure_deriv( tLocalVertexID, aDirection );
                }
            }
            return tDerivative;
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_measure_derivative(
                const moris_index          aDiscretizationMeshIndex,
                const Matrix< DDRMat >    &aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            MORIS_ERROR( false, "mtk::Cell_Cluster::compute_cluster_group_cell_measure_derivative() - Only implemented in child xtk::Cell_Cluster" );
            return 0.0;
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_cell_side_measure(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            MORIS_ERROR( 0, "compute_cluster_cell_side_measure only valid on side clusters" );
            return 0;
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_side_measure(
                const moris_index          aDiscretizationMeshIndex,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            MORIS_ERROR( false, "mtk::Cell_Cluster::compute_cluster_group_cell_side_measure() - Only valid on xtk::Side_Cluster" );
            return 0.0;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat >
        compute_cluster_ig_cell_side_measures(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            MORIS_ERROR( 0, "compute_cluster_ig_cell_side_measures only valid on side clusters" );
            return { { 0 } };
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_cell_side_measure_derivative(
                const Matrix< DDRMat >    &aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid,
                const mtk::Leader_Follower aIsLeader ) const override
        {
            MORIS_ERROR( 0, "compute_cluster_cell_side_measure_derivative only valid on side clusters" );
            return 0;
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_side_measure_derivative(
                const moris_index          aDiscretizationMeshIndex,
                const Matrix< DDRMat >    &aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid,
                const mtk::Leader_Follower aIsLeader ) const override
        {
            MORIS_ERROR( false, "mtk::Cell_Cluster::compute_cluster_group_cell_side_measure_derivative() - Only valid on xtk::Side_Cluster" );
            return 0.0;
        }

        //------------------------------------------------------------------------------

        moris::moris_index
        get_interpolation_cell_index() const override
        {
            return get_interpolation_cell().get_index();
        }

        //------------------------------------------------------------------------------
        // test
        moris::Matrix< moris::IndexMat >
        get_vertex_indices_in_cluster( mtk::Leader_Follower aLeaderFollower = mtk::Leader_Follower::LEADER ) const override
        {
            // number of cells in cluster
            moris::uint tNumVertices = this->get_num_vertices_in_cluster();

            // cell access
            Vector< moris::mtk::Vertex const * > tVertices = this->get_vertices_in_cluster();

            // initialize output
            moris::Matrix< moris::IndexMat > tVertexIndices( 1, tNumVertices );

            // get cell indices and store
            for ( moris::uint i = 0; i < tNumVertices; i++ )
            {
                tVertexIndices( i ) = tVertices( i )->get_index();
            }

            return tVertexIndices;
        }

        Matrix< DDRMat >
        get_quadrature_weights(  ) const
        {
            Matrix< DDRMat > tDummyWeights ;
            return tDummyWeights;
        }

        Matrix< DDRMat >
        get_quadrature_points(  ) const 
        {
            Matrix< DDRMat > tDummyPoints ;
            return tDummyPoints;
        }

        // ##############################################
        //  Cell/Vertex Id Access
        // ##############################################

        moris::Matrix< moris::IdMat >
        get_primary_cell_ids_in_cluster() const override
        {
            // number of cells in cluster
            moris::uint tNumCells = this->get_num_primary_cells();

            // cell access
            Vector< moris::mtk::Cell const * > const &tCells = this->get_primary_cells_in_cluster();

            // initialize output
            moris::Matrix< moris::IdMat > tCellIds( 1, tNumCells );

            // get cell indices and store
            for ( moris::uint i = 0; i < tNumCells; i++ )
            {
                tCellIds( i ) = tCells( i )->get_id();
            }

            return tCellIds;
        }

        //------------------------------------------------------------------------------

        moris::Matrix< moris::IdMat >
        get_void_cell_ids_in_cluster() const override
        {
            MORIS_ERROR( !this->is_trivial(), "get_void_cell_ids_in_cluster on trivial cluster is not allowed" );

            // number of cells in cluster
            moris::uint tNumCells = this->get_num_void_cells();

            // cell access
            Vector< moris::mtk::Cell const * > const &tCells = this->get_void_cells_in_cluster();

            // initialize output
            moris::Matrix< moris::IdMat > tCellIds( 1, tNumCells );

            // get cell indices and store
            for ( moris::uint i = 0; i < tNumCells; i++ )
            {
                tCellIds( i ) = tCells( i )->get_id();
            }

            return tCellIds;
        }

        //------------------------------------------------------------------------------

        moris::moris_id
        get_interpolation_cell_id() const override
        {
            return get_interpolation_cell().get_id();
        }

        //------------------------------------------------------------------------------

        moris::Matrix< moris::IdMat >
        get_vertex_ids_in_cluster() const override
        {
            // number of cells in cluster
            moris::uint tNumVertices = this->get_num_vertices_in_cluster();

            // cell access
            Vector< moris::mtk::Vertex const * > tVertices = this->get_vertices_in_cluster();

            // initialize output
            moris::Matrix< moris::IdMat > tVertexIds( 1, tNumVertices );

            // get cell indices and store
            for ( moris::uint i = 0; i < tNumVertices; i++ )
            {
                tVertexIds( i ) = tVertices( i )->get_id();
            }

            return tVertexIds;
        }

        // ##############################################
        //  Local Coordinate access
        // ##############################################

        /*!
         * Access a primary integration cells parametric coordinates relative to the interpolation cell
         * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
         */
        moris::Matrix< moris::DDRMat >
        get_primary_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aPrimaryCellClusterIndex ) const override
        {
            MORIS_ASSERT( aPrimaryCellClusterIndex < (moris_index)this->get_num_primary_cells(),
                    "Integration Cell Cluster index out of bounds" );

            // get the integration cell of interest
            moris::mtk::Cell const *tIntegrationCell = this->get_primary_cells_in_cluster()( aPrimaryCellClusterIndex );

            // get the vertex pointers on the side
            Vector< moris::mtk::Vertex * > tVerticesOnCell = tIntegrationCell->get_vertex_pointers();

            // allocate output (nnode x dim_xsi)
            moris::Matrix< moris::DDRMat > tVertexParamCoords( tVerticesOnCell.size(), this->get_dim_of_param_coord() );

            // iterate through vertices and collect local coordinates
            for ( moris::uint i = 0; i < tVerticesOnCell.size(); i++ )
            {
                tVertexParamCoords.get_row( i ) = this->get_vertex_local_coordinate_wrt_interp_cell( tVerticesOnCell( i ) ).get_row( 0 );
            }

            return tVertexParamCoords;
        }

        //------------------------------------------------------------------------------

        /*!
         * Access a void integration cells parametric coordinates relative to the interpolation cell
         * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
         */
        moris::Matrix< moris::DDRMat >
        get_void_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aVoidCellClusterIndex ) const override
        {
            MORIS_ERROR( !this->is_trivial(),
                    "get_void_cell_local_coords_on_side_wrt_interp_cell on trivial cluster is not allowed" );

            MORIS_ASSERT( aVoidCellClusterIndex < (moris_index)this->get_num_void_cells(),
                    "Integration Cell Cluster index out of bounds" );

            // get the integration cell of interest
            moris::mtk::Cell const *tIntegrationCell = this->get_void_cells_in_cluster()( aVoidCellClusterIndex );

            // get the vertex pointers on the side
            Vector< moris::mtk::Vertex * > tVerticesOnCell = tIntegrationCell->get_vertex_pointers();

            // allocate output (nnode x dim_xsi)
            moris::Matrix< moris::DDRMat > tVertexParamCoords( tVerticesOnCell.size(), this->get_dim_of_param_coord() );

            // iterate through vertices and collect local coordinates
            for ( moris::uint i = 0; i < tVerticesOnCell.size(); i++ )
            {
                tVertexParamCoords.get_row( i ) = this->get_vertex_local_coordinate_wrt_interp_cell( tVerticesOnCell( i ) ).get_row( 0 );
            }

            return tVertexParamCoords;
        }

        // ##############################################
        //  Size Access
        // ##############################################

        moris::uint
        get_num_primary_cells() const override
        {
            return this->get_primary_cells_in_cluster().size();
        }

        //------------------------------------------------------------------------------

        moris::uint
        get_num_void_cells() const override
        {
            return this->get_void_cells_in_cluster().size();
        }

        //------------------------------------------------------------------------------

        moris::uint
        get_num_vertices_in_cluster() const override
        {
            return this->get_vertices_in_cluster().size();
        }

    };    // class mtk::Cell_Cluster

}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_CL_MTK_CELL_CLUSTER_HPP_ */
