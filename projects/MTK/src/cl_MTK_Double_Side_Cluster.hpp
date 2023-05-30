/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Double_Side_Cluster.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_DOUBLE_SIDE_CLUSTER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_DOUBLE_SIDE_CLUSTER_HPP_

#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Side_Cluster.hpp"

#include "cl_MTK_Cluster.hpp"

namespace moris
{
    namespace mtk
    {
        class Double_Side_Cluster : public Cluster
        {
            //----------------------------------------------------------------
            
            /*
             * An assumption here is made that the first facet in the list of leader facets is paired
             * with the first facet appearing in list of follower facets.
             */

            // Leader side cluster
            moris::mtk::Cluster const *mLeaderSideCluster;

            // Follower side cluster
            moris::mtk::Cluster const *mFollowerSideCluster;

            /*!
             * A one way pairing from leader vertices to follower vertices
             */
            moris::Cell< moris::mtk::Vertex const * > mLeaderToFollowerVertexPairs;

          public:
            // ----------------------------------------------------------------------------------

            /**
             * Default constructor, both leader and follower side clusters are initialized
             * as nullptr.
             */
            Double_Side_Cluster();

            // ----------------------------------------------------------------------------------


            /**
             * Default virtual destructor
             */
            virtual ~Double_Side_Cluster();

            // ----------------------------------------------------------------------------------

            /**
             * Constructor which should be used if to construct operational double side cluster
             * @param[in] aLeaderSideCluster Leader side cluster
             * @param[in] aFollowerSideCluster Follower side cluster
             * @param[in] aLeaderToFollowerVertexPair) Vertices on the sides of the clusters
             */
            Double_Side_Cluster(
                    moris::mtk::Cluster const                       *aLeaderSideCluster,
                    moris::mtk::Cluster const                       *aFollowerSideCluster,
                    moris::Cell< moris::mtk::Vertex const * > const &aLeaderToFollowerVertexPair );

            //----------------------------------------------------------------

            // ##############################################
            //  Side Cluster traits access
            // ##############################################

            //----------------------------------------------------------------

            /*!
             * Ask a cluster in the double side cluster whether it is trivial or not. A trivial
             * cluster has a 1-to-1 relationship between integration entities and interpolation entities.
             * This does not mean there is the same id for each entity.
             * @param[in] aIsLeader enum specifying which cluster the question is for
             * @return True if trivial cluster
             */
            // virtual
            bool
            is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            /*!
             * @return true if leader cluster is trivial
             */
            bool
            is_leader_trivial() const;

            //----------------------------------------------------------------

            /*!
             * @return true if follower cluster is trivial
             */
            bool
            is_follower_trivial() const;

            //----------------------------------------------------------------

            // ##############################################
            //  Single Side Cluster Access
            // ##############################################

            //----------------------------------------------------------------

            /*!
             * @return const leader side cluster
             */
            moris::mtk::Cluster const &
            get_leader_side_cluster() const;

            //----------------------------------------------------------------

            /*!
             * @return const follower side cluster
             */
            moris::mtk::Cluster const &
            get_follower_side_cluster() const;

            //----------------------------------------------------------------

            moris::mtk::Cluster const &
            get_cluster(
                    const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            // ##############################################
            //  Vertex Pair Access
            // ##############################################

            //----------------------------------------------------------------

            /*!
             * @param[in] aLeaderVertex A vertex on the leader side of double side cluster
             * @return Corresponding follower vertex
             */
            moris::mtk::Vertex const *
            get_leader_vertex_pair( moris::mtk::Vertex const *aLeaderVertex ) const;


            //----------------------------------------------------------------

            // ##############################################
            //  Cell Side Ordinals/Vertex Access
            // ##############################################

            //----------------------------------------------------------------

            /*
             * @param[in] aVertex Vertex in mesh
             * @param[in] aIsLeader Leader/Follower selector enum
             * @return Vertex local cluster index wrt leader/follower
             */
            moris_index
            get_vertex_cluster_index(
                    const Vertex           *aVertex,
                    const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            /*
             * @param[in] aIsLeader Leader/Follower selector enum
             * @return Interpolation cell of the side cluster.
             */
            moris::mtk::Cell const &
            get_interpolation_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            /*!
             * @return Interpolation cell of leader side cluster
             */
            moris::mtk::Cell const &
            get_leader_interpolation_cell() const;

            //----------------------------------------------------------------

            /*!
             * @return Interpolation cell of follower side cluster
             */
            moris::mtk::Cell const &
            get_follower_interpolation_cell() const;

            //----------------------------------------------------------------

            /*!
             * @param[in] aIsLeader Leader/Follower selector enum
             * @return Primary integration cells in the side cluster
             */
            moris::Cell< mtk::Cell const * > const &
            get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader ) const;

            //----------------------------------------------------------------

            /*!
             * @return Primary integration cells in the leader side cluster
             */
            moris::Cell< mtk::Cell const * > const &
            get_leader_integration_cells() const;

            //----------------------------------------------------------------

            /*!
             * @return Primary integration cells in the leader side cluster
             */
            moris::Cell< mtk::Cell const * > const &
            get_follower_integration_cells() const;

            //----------------------------------------------------------------
            /*
             * @param[in] aIsLeader Leader/Follower selector enum
             * @return All integration cell side ordinals
             */
            moris::Matrix< moris::IndexMat >
            get_cell_side_ordinals( const mtk::Leader_Follower aIsLeader ) const;

            //----------------------------------------------------------------
            /*
             *
             */
            moris_index
            get_cell_side_ordinal(
                    moris::moris_index      aCellIndexInCluster,
                    const mtk::Leader_Follower aIsLeader ) const;

            /*!
             * @return all integration cell side ordinals on leader side of the
             * double sided side cluster
             */
            moris::Matrix< moris::IndexMat >
            get_leader_integration_cell_side_ordinals() const;

            //----------------------------------------------------------------

            /*!
             * @return all integration cell side ordinals on leader side of the
             * double sided side cluster
             */
            moris_index
            get_leader_cell_side_ordinal( moris::moris_index aLeaderCellIndexInCluster ) const;
            //----------------------------------------------------------------

            /*!
             * @return all integration cell side ordinals on follower side of the
             * double sided side cluster
             */
            moris::Matrix< moris::IndexMat >
            get_follower_integration_cell_side_ordinals() const;

            //----------------------------------------------------------------

            /*!
             * Single side ordinal version of above
             */

            moris_index
            get_follower_cell_side_ordinal( moris::moris_index aFollowerCellIndexInCluster ) const;

            //----------------------------------------------------------------

            /*!
             * @return all the leader vertices in this cluster
             */
            moris::Cell< moris::mtk::Vertex const * >
            get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader ) const;


            moris::Cell< moris::mtk::Vertex const * >
            get_leader_vertices_in_cluster() const;

            //----------------------------------------------------------------

            /*!
             * Returns all the follower vertices in this cluster
             */

            moris::Cell< moris::mtk::Vertex const * >
            get_follower_vertices_in_cluster() const;

            //----------------------------------------------------------------

            moris_index
            get_follower_vertex_ord_on_facet( moris_index aCellClusterIndex,
                    moris::mtk::Vertex const          *aFollowerVertex ) const;

            //----------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_leader_vertex_indices_in_cluster() const;

            //----------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_follower_vertex_indices_in_cluster() const;

            //----------------------------------------------------------------

            // virtual
            moris::Matrix< moris::IndexMat >
            get_vertex_indices_in_cluster() const;

            //----------------------------------------------------------------

            // ##############################################
            //  Local Coordinate Access
            // ##############################################

            //----------------------------------------------------------------

            /*
             * Access the full array of local coordinates on the leader
             */

            moris::Matrix< moris::DDRMat >
            get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader ) const;

            //----------------------------------------------------------------

            moris::Matrix< moris::DDRMat >
            get_leader_vertices_local_coordinates_wrt_interp_cell() const;

            //----------------------------------------------------------------

            /*
             * Access the full array of local coordinates on the follower
             */
            moris::Matrix< moris::DDRMat >
            get_follower_vertices_local_coordinates_wrt_interp_cell() const;

            //----------------------------------------------------------------

            moris::Matrix< moris::DDRMat >
            get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const *aVertex,
                    const mtk::Leader_Follower                                        aIsLeader ) const;
            //----------------------------------------------------------------

            /*
             * Access a single local coordinate of a vertex on the leader
             */
            moris::Matrix< moris::DDRMat >
            get_leader_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const *aVertex ) const;

            //----------------------------------------------------------------

            /*
             * Access a single local coordinate of a vertex on the follower
             */
            moris::Matrix< moris::DDRMat >
            get_follower_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const *aVertex ) const;

            //----------------------------------------------------------------

            moris::Matrix< moris::DDRMat >
            get_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aClusterLocalIndex,
                    const mtk::Leader_Follower                                   aIsLeader ) const;

            //----------------------------------------------------------------

            /*!
             * Access an integration cells parametric coordinates on a side leader side
             * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
             */
            moris::Matrix< moris::DDRMat >
            get_leader_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aLeaderClusterLocalIndex ) const;

            //----------------------------------------------------------------

            /*!
             * Access an integration cells parametric coordinates on a side follower side
             * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
             */
            moris::Matrix< moris::DDRMat >
            get_follower_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aFollowerClusterLocalIndex ) const;

            //----------------------------------------------------------------

            // ##############################################
            //  Size Access
            // ##############################################

            //----------------------------------------------------------------

            /*!
             * Size of the xsi vector of leader
             */
            moris_index
            get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::real
            compute_cluster_cell_measure(
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::real
            compute_cluster_group_cell_measure(
                    const moris_index       aDiscretizationMeshIndex,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            Matrix< DDRMat >
            compute_cluster_ig_cell_measures(
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::real
            compute_cluster_cell_measure_derivative(
                    const Matrix< DDRMat > &aPerturbedVertexCoords,
                    uint                    aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::real
            compute_cluster_group_cell_measure_derivative(
                    const moris_index       aDiscretizationMeshIndex,
                    const Matrix< DDRMat > &aPerturbedVertexCoords,
                    uint                    aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::real
            compute_cluster_cell_side_measure(
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::real
            compute_cluster_group_cell_side_measure(
                    const moris_index       aDiscretizationMeshIndex,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            // ----------------------------------------------------------------------------------

            Matrix< DDRMat >
            compute_cluster_ig_cell_side_measures(
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::real
            compute_cluster_cell_side_measure_derivative(
                    const Matrix< DDRMat > &aPerturbedVertexCoords,
                    uint                    aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::real
            compute_cluster_group_cell_side_measure_derivative(
                    const moris_index       aDiscretizationMeshIndex,
                    const Matrix< DDRMat > &aPerturbedVertexCoords,
                    uint                    aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris_index
            get_leader_dim_of_param_coord() const;

            //----------------------------------------------------------------

            /*!
             * Size of the xsi vector of follower
             */
            moris_index
            get_follower_dim_of_param_coord() const;
            //----------------------------------------------------------------

            moris::uint
            get_leader_num_vertices_in_cluster() const;

            //----------------------------------------------------------------

            moris::uint
            get_follower_num_vertices_in_cluster() const;

            //----------------------------------------------------------------

            /**
             * @brief Get the leader vertex pairs cell
             *
             * @return moris::Cell<moris::mtk::Vertex const *> const& a cell of the leader side vertices
             */
            moris::Cell< moris::mtk::Vertex const * > const &
            get_leader_vertex_pairs() const;

            //----------------------------------------------------------------

            /**
             * @brief memory usage of the double sided cluster
             *
             * @return size_t
             */

            size_t
            capacity();

            //----------------------------------------------------------------

        };    // class Double_Side_Cluster

        //----------------------------------------------------------------

        inline std::ostream &
        operator<<( std::ostream &os, const Double_Side_Cluster &dt )
        {
            os << "\n  Leader Interpolation Cell: " << std::setw( 9 ) << dt.get_leader_side_cluster().get_interpolation_cell().get_id();
            os << " | Follower Interpolation Cell: " << std::setw( 9 ) << dt.get_follower_side_cluster().get_interpolation_cell().get_id();


            moris::Cell< mtk::Cell const * > const &tLeaderIGCells        = dt.get_leader_side_cluster().get_primary_cells_in_cluster();
            moris::Matrix< moris::IndexMat >        tLeaderIGCellSideOrds = dt.get_leader_side_cluster().get_cell_side_ordinals();
            moris::Cell< mtk::Cell const * > const &tFollowerIGCells         = dt.get_follower_side_cluster().get_primary_cells_in_cluster();
            moris::Matrix< moris::IndexMat >        tFollowerIGCellSideOrds  = dt.get_follower_side_cluster().get_cell_side_ordinals();

            os << "       Cell Pairs: " << std::endl;

            for ( moris::uint i = 0; i < tLeaderIGCells.size(); i++ )
            {
                std::cout << "  Leader Cell ID/Ord: " << std::setw( 9 ) << tLeaderIGCells( i )->get_id() << std::setw( 9 ) << tLeaderIGCellSideOrds( i );
                std::cout << "  Follower Cell ID/Ord: " << std::setw( 9 ) << tFollowerIGCells( i )->get_id() << std::setw( 9 ) << tFollowerIGCellSideOrds( i ) << std::endl;
            }

            moris::print( dt.get_vertices_local_coordinates_wrt_interp_cell( mtk::Leader_Follower::LEADER ), "Leader Local Coords" );
            moris::print( dt.get_vertices_local_coordinates_wrt_interp_cell( mtk::Leader_Follower::FOLLOWER ), "Follower Local Coords" );

            return os;
        }

        //----------------------------------------------------------------

        inline std::ostream &
        operator<<( std::ostream &os, Double_Side_Cluster const *const &dt )
        {
            os << *dt;

            return os;
        }

        //----------------------------------------------------------------

    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_CL_MTK_DOUBLE_SIDE_CLUSTER_HPP_ */
