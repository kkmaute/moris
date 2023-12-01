/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Cluster.hpp
 *
 */
#ifndef PROJECTS_MTK_SRC_CL_MTK_CLUSTER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_CLUSTER_HPP_

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"

namespace moris
{
    namespace mtk
    {
        class Cluster_Group;

        class Cluster
        {
          private:
            moris::Cell< moris::mtk::Cell const * > mDummCellCell;

          public:
            Cluster(){};

            virtual ~Cluster(){};

            // ##############################################
            //  Characteristic functions
            // ##############################################

            virtual bool
            is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const = 0;

            // ##############################################
            //  Cell/Vertex Access
            // ##############################################

            virtual moris_index
            get_vertex_cluster_index(
                    const Vertex              *aVertex,
                    const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
            {
                MORIS_ERROR( false, "get_vertex_cluster_index(): not implemented for this cluster type" );
                return 0;
            }

            virtual moris::moris_index
            get_vertex_ordinal_on_facet(
                    moris_index               aCellIndexInCluster,
                    moris::mtk::Vertex const *aVertex ) const
            {
                MORIS_ERROR( false, "get_vertex_ordinal_on_facet(): not implemented for this cluster type" );
                return 0;
            }

            virtual moris::Cell< moris::mtk::Cell const * > const &
            get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const = 0;

            virtual void
            set_interpolation_cell( moris::mtk::Cell const *aInterpCell )
            {
                MORIS_ERROR( false, "set_interpolation_cell(), not implemented for this class" );
            };

            virtual void
            add_primary_integration_cell( moris::Cell< moris::mtk::Cell const * > const &aIntegrationCell )
            {
                MORIS_ERROR( false, "add_primary_integration_cell(), not implemented for this class" );
            };

            virtual void
            add_void_integration_cell( moris::Cell< moris::mtk::Cell const * > const &aIntegrationCell )
            {
                MORIS_ERROR( false, "add_primary_integration_cell(), not implemented for this class" );
            };

            virtual void
            mark_as_nontrivial()
            {
                MORIS_ERROR( false, "mark_as_nontrivial(), not implemented for this class" );
            };

            virtual moris::Cell< moris::mtk::Cell const * > const &
            get_void_cells_in_cluster() const
            {
                MORIS_ERROR( false, "get_void_cells_in_cluster(): not implemented for this cluster type" );
                return mDummCellCell;
            }

            virtual moris::mtk::Cell const &
            get_interpolation_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const = 0;

            virtual moris::Matrix< moris::IndexMat >
            get_cell_side_ordinals( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
            {
                MORIS_ERROR( false, "get_interpolation_cell(): not implemented for this cluster type" );
                return moris::Matrix< moris::IndexMat >( 0, 0 );
            };

            virtual moris_index
            get_cell_side_ordinal(
                    moris::moris_index         aCellIndexInCluster,
                    const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
            {
                MORIS_ERROR( false, "get_cell_side_ordinal(): not implemented for this cluster type" );
                return 0;
            };

            virtual moris::Cell< moris::mtk::Vertex const * >
            get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const = 0;

            virtual moris::mtk::Vertex const *
            get_leader_vertex_pair( moris::mtk::Vertex const *aLeaderVertex ) const
            {
                MORIS_ERROR( false, "get_leader_vertex_pair(): not implemented for this cluster type" );
                return nullptr;
            }

            /**
             * @brief Get the leader vertex pairs, used for double sided cluster
             * @return moris::Cell< moris::mtk::Vertex const *> const& cell of vertices on the leader side
             */
            virtual moris::Cell< moris::mtk::Vertex const * > const &
            get_leader_vertex_pairs() const
            {
                MORIS_ERROR( false, "get_leader_vertex_pair(): not implemented for this cluster type" );
                moris::Cell< moris::mtk::Vertex const * > *tDummyCell = new moris::Cell< moris::mtk::Vertex const * >( 0 );
                return *tDummyCell;
            }

            /**
             * @brief Get the leader side cluster of a double sided cluster
             * @return moris::mtk::Cluster const& leader side cluster
             */
            virtual moris::mtk::Cluster const &
            get_leader_side_cluster() const
            {
                MORIS_ERROR( false, "get_leader_side_cluster(): not implemented for this cluster type" );
                return *this;
            }

            /**
             * @brief Get the follower side cluster of a double sided cluster
             *
             * @return moris::mtk::Cluster const& follower side cluster
             */

            virtual moris::mtk::Cluster const &
            get_follower_side_cluster() const
            {
                MORIS_ERROR( false, "get_follower_side_cluster(): not implemented for this cluster type" );
                return *this;
            }

            virtual moris_index
            get_follower_vertex_ord_on_facet(
                    moris_index               aCellClusterIndex,
                    moris::mtk::Vertex const *aFollowerVertex ) const
            {
                MORIS_ERROR( false, "get_follower_vertex_ord_on_facet(): not implemented for this cluster type" );
                return 0;
            }

            virtual bool
            has_cluster_group( const moris_index aDiscretizationMeshIndex ) const
            {
                MORIS_ERROR( false, "mtk::Cluster::has_cluster_group() - not implemented for this class" );
                return false;
            }

            virtual std::shared_ptr< Cluster_Group >
            get_cluster_group( const moris_index aDiscretizationMeshIndex ) const
            {
                MORIS_ERROR( false, "mtk::Cluster::get_cluster_group() - not implemented for this class" );
                return nullptr;
            }

            virtual void
            set_cluster_group(
                    const moris_index                aDiscretizationMeshIndex,
                    std::shared_ptr< Cluster_Group > aClusterGroupPtr )
            {
                MORIS_ERROR( false, "mtk::Cluster::set_cluster_group() - not implemented for this class" );
            }

            // ##############################################
            //  Local Coordinate Access
            //  (Pure Virtual)
            // ##############################################
            virtual moris::Matrix< moris::DDRMat >
            get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const = 0;

            /*
             * Access a single local coordinate of a vertex
             */
            virtual moris::Matrix< moris::DDRMat >
            get_vertex_local_coordinate_wrt_interp_cell(
                    moris::mtk::Vertex const  *aVertex,
                    const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const = 0;

            virtual moris::Matrix< moris::DDRMat >
            get_cell_local_coords_on_side_wrt_interp_cell(
                    moris::moris_index         aLeaderClusterLocalIndex,
                    const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
            {
                MORIS_ERROR( false, "get_cell_local_coords_on_side_wrt_interp_cell(): not implemented for this cluster type" );
                return moris::Matrix< moris::DDRMat >( 0, 0 );
            }

            // ##############################################
            //  Size Access
            //  (Pure Virtual)
            // ##############################################
            /*!
             * Size of the xsi vector in this side cluster
             */
            virtual moris_index
            get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const = 0;

            /*
             * Compute the measure (volume 3d or area 2d) of the cells in the void or primary phase
             */
            virtual moris::real
            compute_cluster_cell_measure(
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            virtual moris::real
            compute_cluster_group_cell_measure(
                    const moris_index          aDiscretizationMeshIndex,
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;


            /*
             * Compute the measure of individual IG cells (volume 3d or area 2d) in the void or primary phase
             */
            virtual Matrix< DDRMat >
            compute_cluster_ig_cell_measures(
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            /*
             * Compute the derivative of measure (volume 3d or area 2d) of the cells in the void or primary phase
             * wrt a given vertex and space direction
             */
            virtual moris::real
            compute_cluster_cell_measure_derivative(
                    const Matrix< DDRMat >    &aPerturbedVertexCoords,
                    uint                       aDirection,
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            virtual moris::real
            compute_cluster_group_cell_measure_derivative(
                    const moris_index          aDiscretizationMeshIndex,
                    const Matrix< DDRMat >    &aPerturbedVertexCoords,
                    uint                       aDirection,
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            /*
             * Compute the side measure (surface area 3d or length 2d) of the cells in the void or primary phase on the side set.
             * Only valid on side cluster type mtk clusters
             */
            virtual moris::real
            compute_cluster_cell_side_measure(
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            virtual moris::real
            compute_cluster_group_cell_side_measure(
                    const moris_index          aDiscretizationMeshIndex,
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            /*
             * Compute vector of individual side measures for each IG cell
             */
            virtual Matrix< DDRMat >
            compute_cluster_ig_cell_side_measures(
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            /*
             * Compute the derivative side measure (surface area 3d or length 2d) of the cells in the void or primary phase on the side set.
             * wrt a given vertex and space direction
             * Only valid on side cluster type mtk clusters
             */
            virtual moris::real
            compute_cluster_cell_side_measure_derivative(
                    const Matrix< DDRMat >    &aPerturbedVertexCoords,
                    uint                       aDirection,
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            virtual moris::real
            compute_cluster_group_cell_side_measure_derivative(
                    const moris_index          aDiscretizationMeshIndex,
                    const Matrix< DDRMat >    &aPerturbedVertexCoords,
                    uint                       aDirection,
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            // ---------------------------------------------
            // EVERYTHING BELOW THIS LINE HAS A DEFAULT
            // IMPLEMENTATION
            // ---------------------------------------------

            // ##############################################
            //  Cell/Vertex Index Access
            // ##############################################

            virtual moris::Matrix< moris::IndexMat >
            get_primary_cell_indices_in_cluster() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_primary_cell_indices_in_cluster() - not implemented for this cluster type" );
                return moris::Matrix< moris::IndexMat >( 0, 0 );
            }

            virtual moris::Matrix< moris::IndexMat >
            get_void_cell_indices_in_cluster() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_void_cell_indices_in_cluster() - not implemented for this cluster type" );
                return moris::Matrix< moris::IndexMat >( 0, 0 );
            }

            virtual moris::moris_index
            get_interpolation_cell_index() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_interpolation_cell_index() - not implemented for this cluster type" );
                return moris::moris_index( 0 );
            }

            virtual moris::Matrix< moris::IndexMat >
            get_vertex_indices_in_cluster( mtk::Leader_Follower aLeaderFollower = mtk::Leader_Follower::LEADER ) const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_vertex_indices_in_cluster() - not implemented for this cluster type" );
                return moris::Matrix< moris::IndexMat >( 0, 0 );
            }

            // ##############################################
            //  Cell/Vertex Id Access
            // ##############################################

            virtual moris::Matrix< moris::IdMat >
            get_primary_cell_ids_in_cluster() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_primary_cell_ids_in_cluster() - not implemented for this cluster type" );
                return moris::Matrix< moris::IndexMat >( 0, 0 );
            }

            virtual moris::Matrix< moris::IdMat >
            get_void_cell_ids_in_cluster() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_void_cell_ids_in_cluster() - not implemented for this cluster type" );
                return moris::Matrix< moris::IdMat >( 0, 0 );
            }

            virtual moris::moris_id
            get_interpolation_cell_id() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_interpolation_cell_id() - not implemented for this cluster type" );
                return 0;
            }

            virtual moris::Matrix< moris::IdMat >
            get_vertex_ids_in_cluster() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_vertex_ids_in_cluster() - not implemented for this cluster type" );
                return moris::Matrix< moris::IdMat >( 0, 0 );
            }

            virtual moris::Matrix< moris::DDRMat >
            get_vertex_coords_in_cluster() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_vertex_ids_in_cluster() - not implemented for this cluster type" );
                return moris::Matrix< moris::DDRMat >( 0, 0 );
            }

            // ##############################################
            //  Local Coordinate access
            // ##############################################

            /*!
             * Access a primary integration cells parametric coordinates relative to the interpolation cell
             * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
             */
            virtual moris::Matrix< moris::DDRMat >
            get_primary_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aPrimaryCellClusterIndex ) const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell() - not implemented for this cluster type." );
                return moris::Matrix< moris::DDRMat >( 0, 0 );
            }

            /*!
             * Access a void integration cells parametric coordinates relative to the interpolation cell
             * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
             */
            virtual moris::Matrix< moris::DDRMat >
            get_void_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aVoidCellClusterIndex ) const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_void_cell_local_coords_on_side_wrt_interp_cell() - not implemented for this cluster type." );
                return moris::Matrix< moris::DDRMat >( 0, 0 );
            }

            // ##############################################
            //  Size Access
            // ##############################################

            virtual moris::uint
            get_num_primary_cells() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_num_primary_cells() - not implemented for this cluster type" );
                return 0;
            }

            virtual moris::uint
            get_num_void_cells() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_num_void_cells() - not implemented for this cluster type" );
                return 0;
            }

            virtual moris::uint
            get_num_vertices_in_cluster() const
            {
                MORIS_ERROR( false, "MTK::Cluster::get_num_vertices_in_cluster() - not implemented for this cluster type" );
                return 0;
            }

            moris::Cell< moris::mtk::Vertex * >
            get_primary_vertices_in_cluster() const
            {
                // initialize output array
                moris::Cell< moris::mtk::Vertex * > tPrimaryVertices;

                // get the primary cell
                moris::Cell< moris::mtk::Cell const * > const &tPrimaryCells = this->get_primary_cells_in_cluster();

                // initialize map marking which vertices found to prevent listing them multiple times
                map< moris_index, moris_index > tPrimaryVertexIndices;
                moris_index                     tPrimaryVertexCounter = 0;

                // collect the vertices from all primary cells
                for ( uint iPrimaryCell = 0; iPrimaryCell < tPrimaryCells.size(); iPrimaryCell++ )
                {
                    // get the vertices from the cells
                    moris::Cell< mtk::Vertex * > tVerticesOnCell = tPrimaryCells( iPrimaryCell )->get_vertex_pointers();

                    for ( uint iVertOnCell = 0; iVertOnCell < tVerticesOnCell.size(); iVertOnCell++ )
                    {
                        // get the index of the current vertex
                        moris_index tVertIndex = tVerticesOnCell( iVertOnCell )->get_index();

                        // check if this vertex has already been listed, if not add this vertex to the list
                        if ( !tPrimaryVertexIndices.key_exists( tVertIndex ) )
                        {
                            tPrimaryVertices.push_back( tVerticesOnCell( iVertOnCell ) );
                            tPrimaryVertexIndices[ tVertIndex ] = tPrimaryVertexCounter++;
                        }
                    }
                }

                // return list of primary vertices
                return tPrimaryVertices;
            }

            virtual Matrix< IndexMat >
            get_hanging_nodes() const
            {
                return Matrix< IndexMat >( 0, 0 );
            };

        };    // class MTK_Cluster

    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_CL_MTK_CLUSTER_HPP_ */
