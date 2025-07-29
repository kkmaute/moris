/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Double_Side_Cluster.cpp
 *
 */

#include "cl_MTK_Double_Side_Cluster.hpp"
#include "fn_TOL_Capacities.hpp"

//----------------------------------------------------------------------------------

namespace moris::mtk
{
    Double_Side_Cluster::Double_Side_Cluster()
            : mLeaderSideCluster( nullptr )
            , mFollowerSideCluster( nullptr )
    {
    }

    //----------------------------------------------------------------------------------

    Double_Side_Cluster::~Double_Side_Cluster() {}

    //----------------------------------------------------------------------------------

    Double_Side_Cluster::Double_Side_Cluster(
            moris::mtk::Cluster const                  *aLeaderSideCluster,
            moris::mtk::Cluster const                  *aFollowerSideCluster,
            Vector< moris::mtk::Vertex const * > const &aLeftToRightVertexPair )
            : mLeaderSideCluster( aLeaderSideCluster )
            , mFollowerSideCluster( aFollowerSideCluster )
    {
        mLeaderToFollowerVertexPairs.append( aLeftToRightVertexPair );
    }

    //----------------------------------------------------------------------------------

    bool
    Double_Side_Cluster::is_trivial(
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().is_trivial();
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().is_trivial();
        }
        else
        {
            MORIS_ERROR( false, "is_trivial(): can only be LEADER or FOLLOWER" );
            return false;
        }
    }

    //----------------------------------------------------------------------------------

    bool
    Double_Side_Cluster::is_leader_trivial() const
    {
        return this->get_leader_side_cluster().is_trivial();
    }

    //----------------------------------------------------------------------------------

    bool
    Double_Side_Cluster::is_follower_trivial() const
    {
        return this->get_follower_side_cluster().is_trivial();
    }

    //----------------------------------------------------------------------------------

    moris::mtk::Cluster const &
    Double_Side_Cluster::get_leader_side_cluster() const
    {
        return *mLeaderSideCluster;
    }

    //----------------------------------------------------------------------------------

    moris::mtk::Cluster const &
    Double_Side_Cluster::get_follower_side_cluster() const
    {
        return *mFollowerSideCluster;
    }

    //----------------------------------------------------------------------------------

    moris::mtk::Cluster const &
    Double_Side_Cluster::get_cluster(
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster();
        }
        else
        {
            return this->get_follower_side_cluster();
        }
    }

    //----------------------------------------------------------------------------------

    moris::mtk::Vertex const *
    Double_Side_Cluster::get_leader_vertex_pair(
            moris::mtk::Vertex const *aLeaderVertex ) const
    {
        moris_index tLeaderClusterIndex = this->get_leader_side_cluster().get_vertex_cluster_index( aLeaderVertex );

        MORIS_ASSERT( tLeaderClusterIndex < (moris_index)mLeaderToFollowerVertexPairs.size(), "Vertex index out of bounds in pairing." );

        return mLeaderToFollowerVertexPairs( tLeaderClusterIndex );
    }

    //----------------------------------------------------------------------------------

    Vector< moris::mtk::Vertex const * > const &
    Double_Side_Cluster::get_leader_vertex_pairs() const
    {
        return mLeaderToFollowerVertexPairs;
    }

    //----------------------------------------------------------------------------------

    moris_index
    Double_Side_Cluster::get_vertex_cluster_index(
            const Vertex              *aVertex,
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().get_vertex_cluster_index( aVertex );
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().get_vertex_cluster_index( aVertex );
        }
        else
        {
            MORIS_ERROR( false, "get_vertex_cluster_index(): can only be LEADER and FOLLOWER" );
            return 0;
        }
    }

    //----------------------------------------------------------------------------------

    moris::mtk::Cell const &
    Double_Side_Cluster::get_interpolation_cell(
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().get_interpolation_cell();
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().get_interpolation_cell();
        }
        else
        {
            MORIS_ERROR( false, "get_interpolation_cell(): can only be 0 and 1" );

            // This function will be never be used
            return this->get_leader_side_cluster().get_interpolation_cell();
        }
    }

    //----------------------------------------------------------------------------------

    moris::mtk::Cell const &
    Double_Side_Cluster::get_leader_interpolation_cell() const
    {
        return this->get_leader_side_cluster().get_interpolation_cell();
    }

    //----------------------------------------------------------------------------------

    moris::mtk::Cell const &
    Double_Side_Cluster::get_follower_interpolation_cell() const
    {
        return this->get_follower_side_cluster().get_interpolation_cell();
    }

    //----------------------------------------------------------------------------------

    Vector< mtk::Cell const * > const &
    Double_Side_Cluster::get_primary_cells_in_cluster(
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().get_primary_cells_in_cluster();
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().get_primary_cells_in_cluster();
        }
        else
        {
            MORIS_ERROR( false, "get_primary_cells_in_cluster(): can only be LEADER and FOLLOWER" );

            // create a dummy cell that never will be generated
            Vector< mtk::Cell const * > *tDummyCell = new Vector< mtk::Cell const * >( 0 );
            return *tDummyCell;
        }
    }

    //----------------------------------------------------------------------------------

    Vector< mtk::Cell const * > const &
    Double_Side_Cluster::get_leader_integration_cells() const
    {
        return this->get_leader_side_cluster().get_primary_cells_in_cluster();
    }

    //----------------------------------------------------------------------------------

    Vector< mtk::Cell const * > const &
    Double_Side_Cluster::get_follower_integration_cells() const
    {
        return this->get_follower_side_cluster().get_primary_cells_in_cluster();
    }

    //----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Double_Side_Cluster::get_cell_side_ordinals(
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().get_cell_side_ordinals();
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().get_cell_side_ordinals();
        }
        else
        {
            MORIS_ERROR( false, "get_cell_side_ordinals(): can only be LEADER and FOLLOWER" );
            return Matrix< IndexMat >( 0, 0 );
        }
    }

    //----------------------------------------------------------------------------------

    moris_index
    Double_Side_Cluster::get_cell_side_ordinal(
            moris::moris_index         aCellIndexInCluster,
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().get_cell_side_ordinal( aCellIndexInCluster );
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().get_cell_side_ordinal( aCellIndexInCluster );
        }
        else
        {
            MORIS_ERROR( false, "Double_Side_Cluster::get_cell_side_ordinal() - can only be LEADER and FOLLOWER" );
            return 0;
        }
    }

    //----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Double_Side_Cluster::get_leader_integration_cell_side_ordinals() const
    {
        return this->get_leader_side_cluster().get_cell_side_ordinals();
    }

    //----------------------------------------------------------------------------------

    moris_index
    Double_Side_Cluster::get_leader_cell_side_ordinal(
            moris::moris_index aLeaderCellIndexInCluster ) const
    {
        return this->get_leader_side_cluster().get_cell_side_ordinal( aLeaderCellIndexInCluster );
    }

    //----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Double_Side_Cluster::get_follower_integration_cell_side_ordinals() const
    {
        return this->get_follower_side_cluster().get_cell_side_ordinals();
    }

    //----------------------------------------------------------------------------------

    moris_index
    Double_Side_Cluster::get_follower_cell_side_ordinal(
            moris::moris_index aFollowerCellIndexInCluster ) const
    {
        return this->get_follower_side_cluster().get_cell_side_ordinal( aFollowerCellIndexInCluster );
    }

    //----------------------------------------------------------------------------------

    Vector< moris::mtk::Vertex const * >
    Double_Side_Cluster::get_vertices_in_cluster(
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().get_vertices_in_cluster();
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().get_vertices_in_cluster();
        }
        else
        {
            MORIS_ERROR( false, "Double_Side_Cluster::get_vertices_in_cluster() - Can only be LEADER or FOLLOWER" );

            return Vector< moris::mtk::Vertex const * >( 0 );
        }
    }

    //----------------------------------------------------------------------------------

    Vector< moris::mtk::Vertex const * >
    Double_Side_Cluster::get_leader_vertices_in_cluster() const
    {
        return this->get_leader_side_cluster().get_vertices_in_cluster();
    }

    //----------------------------------------------------------------------------------

    Vector< moris::mtk::Vertex const * >
    Double_Side_Cluster::get_follower_vertices_in_cluster() const
    {
        return this->get_follower_side_cluster().get_vertices_in_cluster();
    }

    //----------------------------------------------------------------------------------

    void
    Double_Side_Cluster::compute_mapped_quadrature_points_and_weights_leader( const Matrix< DDRMat > aIntegrationPointsLeader, const Matrix< DDRMat > aIntegrationWeightsLeader )
    {
        // // Get IP element 
        // moris::mtk::Cell const & tInterpolationCell = this->get_leader_interpolation_cell();

        // // Create IP cell interpolation rule
        // mtk::Interpolation_Rule tIPInterpolationRule( tInterpolationCell.get_geometry_type() , mtk::Interpolation_Type::LAGRANGE 
        //         ,tInterpolationCell.get_interpolation_order() , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );
        
        // // Get IG elements
        // Vector< moris::mtk::Cell const * > const& tIGCells =  this->get_leader_integration_cells();

        // // Set the size of the quadrature points and weights matrices
        // mQuadraturePointsLeader.reshape( aIntegrationPointsLeader.n_rows() , tIGCells.size() * aIntegrationPointsLeader.n_cols() );

        // mQuadratureWeightsLeader.reshape( aIntegrationWeightsLeader.n_rows() , tIGCells.size() * aIntegrationWeightsLeader.n_cols() );

        // // Loop over IG Elements
        // for( uint iPrimaryCell = 0; iPrimaryCell < tIGCells.size() ; iPrimaryCell++ )
        // {          
        //     // Get cell info from tIGCell (required by space interpolator)
        //     const mtk::Cell_Info* tCellInfo = tIGCells( iPrimaryCell )->get_cell_info();
            
        //     // Create Interpolation Rule
        //     mtk::Interpolation_Rule tIGInterpolationRule( tIGCells( iPrimaryCell )->get_geometry_type() , mtk::Interpolation_Type::LAGRANGE 
        //         ,tIGCells( iPrimaryCell )->get_interpolation_order() , mtk::Geometry_Type::LINE , mtk::Interpolation_Type::LAGRANGE , mtk::Interpolation_Order::LINEAR );

        //     // Create space interpolator
        //     mtk::Space_Interpolator tIGSpaceInterpolator =  mtk::Space_Interpolator(
        //         tIGInterpolationRule,
        //         tIPInterpolationRule,
        //         tCellInfo->compute_cell_shape( tIGCells( iPrimaryCell ) ),
        //         true ); // since the cell cluster is a sideset this is true

        //     // Get ordinal of side cluster w.r.t. IG cell
        //     moris_index tOrdinal = this->get_leader_cell_side_ordinal( iPrimaryCell );
            
        //     // Get vertices attached to that ordinal
        //     Vector< const mtk::Vertex *> tClusterVertices = tIGCells( iPrimaryCell )->get_geometric_vertices_on_side_ordinal( tOrdinal )

        //     // Declare the local coord matrix
        //     Matrix< DDRMat > tLocalCoords;
        //     tLocalCoords.reshape( tClusterVertices.size(), 2 ) 

        //     // Now get the coordinates of the vertex w.r.t. IP cell
        //     for( uint iVertices = 0 ; iVertices < tClusterVertices.size(); iVertices++ )
        //     {
        //         // get local vertex coordinates
        //         Matrix< DDRMat > tLocalCoord = this->get_leader_vertex_local_coordinate_wrt_interp_cell( tClusterVertices( iVertices ) );

        //         // add them to the total list
        //         tLocalCoords.set_row( iVertices , tLocalCoord );

        //     }
        //     tIGSpaceInterpolator.set_space_param_coeff( tCellLocalCoords );

        //     // Set coeffs for the jacobian
        //     tIGSpaceInterpolator.set_space_coeff( tCellLocalCoords );

        //     // Loop over all quadrature points in local IG element space to IP element local space
        //     for( uint iQuadPoint = 0 ; iQuadPoint < aIntegrationWeights.numel(); iQuadPoint++ )
        //     {
        //         // Set local coord point
        //         tIGSpaceInterpolator.set_space_time( aIntegrationPoints.get_column( iQuadPoint ) );
                
        //         // get mapped quadrature point (from IG element parent element space to IP element parent element space)
        //         Matrix< DDRMat > tMappedQuadraturePoint = tIGSpaceInterpolator.map_integration_point();

        //         // Place in quadrature point vector
        //         mQuadraturePoints.set_column( iPrimaryCell * aIntegrationWeights.numel() + iQuadPoint , ( tMappedQuadraturePoint ) );

        //         // Get determinant
        //         real tDetJ = tIGSpaceInterpolator.space_det_J();

        //         // Get modified weight
        //         real tWStar = aIntegrationWeights( iQuadPoint ) * tDetJ ;

        //         // Set quadrature point as multiplied with det_J
        //         mQuadratureWeights( iPrimaryCell * aIntegrationWeights.numel() + iQuadPoint ) = tWStar ;
                
        //     }

        // }

    }


    //----------------------------------------------------------------------------------

    moris_index
    Double_Side_Cluster::get_follower_vertex_ord_on_facet(
            moris_index               aCellClusterIndex,
            moris::mtk::Vertex const *aFollowerVertex ) const
    {
        return mFollowerSideCluster->get_vertex_ordinal_on_facet( aCellClusterIndex, aFollowerVertex );
    }

    //----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Double_Side_Cluster::get_leader_vertex_indices_in_cluster() const
    {
        return this->get_leader_side_cluster().get_vertex_indices_in_cluster();
    }

    //----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Double_Side_Cluster::get_follower_vertex_indices_in_cluster() const
    {
        return this->get_follower_side_cluster().get_vertex_indices_in_cluster();
    }

    //----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Double_Side_Cluster::get_vertex_indices_in_cluster( mtk::Leader_Follower aLeaderFollower ) const
    {
        switch ( aLeaderFollower )
        {
            // only vertices of the leader element
            case mtk::Leader_Follower::LEADER:
            {
                return this->get_leader_vertex_indices_in_cluster();
            }

            // only vertices of the follower element
            case mtk::Leader_Follower::FOLLOWER:
            {
                return this->get_follower_vertex_indices_in_cluster();
            }

            // use "UNDEFINED" to return a list of all vertices on the side cluster
            case mtk::Leader_Follower::UNDEFINED:
            {
                // number of cells in cluster
                uint tLeaderNumVertices   = mLeaderSideCluster->get_num_vertices_in_cluster();
                uint tFollowerNumVertices = mFollowerSideCluster->get_num_vertices_in_cluster();

                // access the vertices in a side cluster
                Vector< moris::mtk::Vertex const * > const &tLeaderVertices   = mLeaderSideCluster->get_vertices_in_cluster();
                Vector< moris::mtk::Vertex const * > const &tFollowerVertices = mFollowerSideCluster->get_vertices_in_cluster();

                // initialize output
                Matrix< IndexMat > tVertexIndices( 1, tLeaderNumVertices + tFollowerNumVertices );

                // get cell indices and store
                for ( uint iLeaderVertex = 0; iLeaderVertex < tLeaderNumVertices; iLeaderVertex++ )
                {
                    tVertexIndices( iLeaderVertex ) = tLeaderVertices( iLeaderVertex )->get_index();
                }
                for ( uint iFollowerVertex = 0; iFollowerVertex < tFollowerNumVertices; iFollowerVertex++ )
                {
                    tVertexIndices( tLeaderNumVertices + iFollowerVertex ) = tFollowerVertices( iFollowerVertex )->get_index();
                }

                // return the concatenated list
                return tVertexIndices;
            }
            default:
            {
                MORIS_ERROR( false, "Double_Side_Cluster::get_vertex_indices_in_cluster() - Unknown Leader_Follower ENUM." );
                return Matrix< IndexMat >( 1, 1, -1 );
            }

        }    // end switch: Leader_Follower enum

    }    // end function: Double_Side_Cluster::get_vertex_indices_in_cluster()

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().get_vertices_local_coordinates_wrt_interp_cell();
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().get_vertices_local_coordinates_wrt_interp_cell();
        }
        else
        {
            MORIS_ERROR( false, "get_vertices_local_coordinates_wrt_interp_cell(): can only be LEADER and FOLLOWER" );

            return { {} };
        }
    }

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::get_leader_vertices_local_coordinates_wrt_interp_cell() const
    {
        return this->get_leader_side_cluster().get_vertices_local_coordinates_wrt_interp_cell();
    }

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::get_follower_vertices_local_coordinates_wrt_interp_cell() const
    {
        return this->get_follower_side_cluster().get_vertices_local_coordinates_wrt_interp_cell();
    }

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::get_vertex_local_coordinate_wrt_interp_cell(
            moris::mtk::Vertex const  *aVertex,
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().get_vertex_local_coordinate_wrt_interp_cell( aVertex );
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().get_vertex_local_coordinate_wrt_interp_cell( aVertex );
        }
        else
        {
            MORIS_ERROR( false, "get_vertex_local_coordinate_wrt_interp_cell(): can only be LEADER and FOLLOWER" );
            return Matrix< DDRMat >( 0, 0 );
        }
    }

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::get_leader_vertex_local_coordinate_wrt_interp_cell(
            moris::mtk::Vertex const *aVertex ) const
    {
        return this->get_leader_side_cluster().get_vertex_local_coordinate_wrt_interp_cell( aVertex );
    }

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::get_follower_vertex_local_coordinate_wrt_interp_cell(
            moris::mtk::Vertex const *aVertex ) const
    {
        return this->get_follower_side_cluster().get_vertex_local_coordinate_wrt_interp_cell( aVertex );
    }

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::get_cell_local_coords_on_side_wrt_interp_cell(
            moris::moris_index         aClusterLocalIndex,
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().get_cell_local_coords_on_side_wrt_interp_cell( aClusterLocalIndex );
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().get_cell_local_coords_on_side_wrt_interp_cell( aClusterLocalIndex );
        }
        else
        {
            MORIS_ERROR( false, "get_cell_local_coords_on_side_wrt_interp_cell(): can only be LEADER and FOLLOWER" );
            return Matrix< DDRMat >( 0, 0 );
        }
    }

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::get_leader_cell_local_coords_on_side_wrt_interp_cell(
            moris::moris_index aLeaderClusterLocalIndex ) const
    {
        return this->get_leader_side_cluster().get_cell_local_coords_on_side_wrt_interp_cell( aLeaderClusterLocalIndex );
    }

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::get_follower_cell_local_coords_on_side_wrt_interp_cell(
            moris::moris_index aFollowerClusterLocalIndex ) const
    {
        return this->get_follower_side_cluster().get_cell_local_coords_on_side_wrt_interp_cell( aFollowerClusterLocalIndex );
    }

    //----------------------------------------------------------------------------------

    moris_index
    Double_Side_Cluster::get_dim_of_param_coord(
            const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aIsLeader == mtk::Leader_Follower::LEADER )
        {
            return this->get_leader_side_cluster().get_dim_of_param_coord();
        }
        else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
        {
            return this->get_follower_side_cluster().get_dim_of_param_coord();
        }
        else
        {
            MORIS_ERROR( false, "get_dim_of_param_coord(): can only be LEADER and FOLLOWER" );
            return 0;
        }
    }

    //----------------------------------------------------------------------------------

    moris::real
    Double_Side_Cluster::compute_cluster_cell_measure(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cluster const &tCluster = this->get_cluster( aIsLeader );

        return tCluster.compute_cluster_cell_measure( aPrimaryOrVoid, aIsLeader );
    }

    //----------------------------------------------------------------------------------

    moris::real
    Double_Side_Cluster::compute_cluster_group_cell_measure(
            const moris_index          aDiscretizationMeshIndex,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cluster const &tCluster = this->get_cluster( aIsLeader );

        return tCluster.compute_cluster_group_cell_measure( aDiscretizationMeshIndex, aPrimaryOrVoid, aIsLeader );
    }

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::compute_cluster_ig_cell_measures(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cluster const &tCluster = this->get_cluster( aIsLeader );

        return tCluster.compute_cluster_ig_cell_measures( aPrimaryOrVoid, aIsLeader );
    }

    //----------------------------------------------------------------------------------

    moris::real
    Double_Side_Cluster::compute_cluster_cell_measure_derivative(
            const Matrix< DDRMat >    &aPerturbedVertexCoords,
            uint                       aDirection,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cluster const &tCluster = this->get_cluster( aIsLeader );
        return tCluster.compute_cluster_cell_measure_derivative(
                aPerturbedVertexCoords,
                aDirection,
                aPrimaryOrVoid,
                aIsLeader );
    }

    //----------------------------------------------------------------------------------

    moris::real
    Double_Side_Cluster::compute_cluster_group_cell_measure_derivative(
            const moris_index          aDiscretizationMeshIndex,
            const Matrix< DDRMat >    &aPerturbedVertexCoords,
            uint                       aDirection,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cluster const &tCluster = this->get_cluster( aIsLeader );
        return tCluster.compute_cluster_group_cell_measure_derivative(
                aDiscretizationMeshIndex,
                aPerturbedVertexCoords,
                aDirection,
                aPrimaryOrVoid,
                aIsLeader );
    }

    //----------------------------------------------------------------------------------

    moris::real
    Double_Side_Cluster::compute_cluster_cell_side_measure(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cluster const &tCluster = this->get_cluster( aIsLeader );

        return tCluster.compute_cluster_cell_side_measure( aPrimaryOrVoid, aIsLeader );
    }

    //----------------------------------------------------------------------------------

    moris::real
    Double_Side_Cluster::compute_cluster_group_cell_side_measure(
            const moris_index          aDiscretizationMeshIndex,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cluster const &tCluster = this->get_cluster( aIsLeader );

        return tCluster.compute_cluster_group_cell_side_measure( aDiscretizationMeshIndex, aPrimaryOrVoid, aIsLeader );
    }

    //----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Double_Side_Cluster::compute_cluster_ig_cell_side_measures(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cluster const &tCluster = this->get_cluster( aIsLeader );

        return tCluster.compute_cluster_ig_cell_side_measures( aPrimaryOrVoid, aIsLeader );
    }

    //----------------------------------------------------------------

    moris::real
    Double_Side_Cluster::compute_cluster_cell_side_measure_derivative(
            const Matrix< DDRMat >    &aPerturbedVertexCoords,
            uint                       aDirection,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cluster const &tCluster = this->get_cluster( aIsLeader );
        return tCluster.compute_cluster_cell_side_measure_derivative(
                aPerturbedVertexCoords,
                aDirection,
                aPrimaryOrVoid,
                aIsLeader );
    }

    //----------------------------------------------------------------

    moris::real
    Double_Side_Cluster::compute_cluster_group_cell_side_measure_derivative(
            const moris_index          aDiscretizationMeshIndex,
            const Matrix< DDRMat >    &aPerturbedVertexCoords,
            uint                       aDirection,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cluster const &tCluster = this->get_cluster( aIsLeader );
        return tCluster.compute_cluster_group_cell_side_measure_derivative(
                aDiscretizationMeshIndex,
                aPerturbedVertexCoords,
                aDirection,
                aPrimaryOrVoid,
                aIsLeader );
    }

    //----------------------------------------------------------------------------------

    moris_index
    Double_Side_Cluster::get_leader_dim_of_param_coord() const
    {
        return this->get_leader_side_cluster().get_dim_of_param_coord();
    }

    //----------------------------------------------------------------------------------

    moris_index
    Double_Side_Cluster::get_follower_dim_of_param_coord() const
    {
        return this->get_follower_side_cluster().get_dim_of_param_coord();
    }

    //----------------------------------------------------------------------------------

    uint
    Double_Side_Cluster::get_leader_num_vertices_in_cluster() const
    {
        return this->get_leader_side_cluster().get_num_vertices_in_cluster();
    }

    //----------------------------------------------------------------------------------

    uint
    Double_Side_Cluster::get_follower_num_vertices_in_cluster() const
    {
        return this->get_follower_side_cluster().get_num_vertices_in_cluster();
    }

    //----------------------------------------------------------------------------------

    size_t
    Double_Side_Cluster::capacity()
    {
        size_t tCapacity = 0;

        // sum up the member data size
        tCapacity += sizeof( mLeaderSideCluster );
        tCapacity += sizeof( mFollowerSideCluster );
        tCapacity += mLeaderToFollowerVertexPairs.capacity() * ( ( sizeof( void * ) ) + 1 );

        return tCapacity;
    }
}    // namespace moris::mtk
