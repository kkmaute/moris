/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MIG_Periodic_2D.cpp
 *
 */

#include "cl_MIG_Periodic_2D.hpp"

#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_MTK_Side_Cluster.hpp"
#include "moris_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_unique.hpp"
#include <unordered_map>

#include "fn_norm.hpp"
#include "fn_sum.hpp"
#include "fn_inv.hpp"
#include "fn_dot.hpp"
#include "fn_sort.hpp"
#include "fn_rank.hpp"
#include "cl_MTK_Mesh.hpp"
#include "op_equal_equal.hpp"
#include "op_less.hpp"
#include "fn_find.hpp"
#include "fn_join_horiz.hpp"
#include "fn_trans.hpp"

#include "cl_MTK_Cell_ISC.hpp"
#include "cl_MTK_Vertex_ISC_Impl.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_Interpolation.hpp"
#include "cl_MTK_Side_Cluster_ISC_Impl.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_Tracer.hpp"

namespace moris::mig
{
    Periodic_2D::Periodic_2D(
            std::shared_ptr< moris::mtk::Mesh_Manager > aMeshManager,
            moris::moris_index                          aMeshIndex,
            moris::ParameterList                       &aParameterList,
            moris::uint                                 aNumBulkPhases )
            : mMeshManager( aMeshManager )
            , mMeshIndex( aMeshIndex )
            , mNumBulkPhases( aNumBulkPhases )
    {
        // get the periodic mesh set names
        std::string tMeshSideSetNames = aParameterList.get< std::string >( "periodic_side_set_pair" );

        // store the names in mMeshSideSetPairs
        string_to_cell_of_cell( aParameterList.get< std::string >( "periodic_side_set_pair" ), mMeshSideSetPairs );
    }

    //---------------------------------------------------------------------

    Periodic_2D::~Periodic_2D()
    {
    }

    //---------------------------------------------------------------------

    void
    Periodic_2D::perform()
    {
        this->determine_max_size();

        // initialize phase interaction table
        Matrix< IndexMat > tPhaseInteractionTable( mNumBulkPhases, mNumBulkPhases );

        for ( uint Ii = 0; Ii < mNumBulkPhases; Ii++ )
        {
            for ( uint Ij = 0; Ij < mNumBulkPhases; Ij++ )
            {
                tPhaseInteractionTable( Ii, Ij ) = Ij + Ii * mNumBulkPhases;
            }
        }

        // loop over the pairs
        for ( uint tPairCount = 0; tPairCount < mMeshSideSetPairs.size(); tPairCount++ )
        {
            // Keep track of time
            Tracer tTracer( "MTK", "No Type", "Searching For Pairs" );

            // permutation pair indicates the indicial relation between 3D and surface coordinates
            // namely depending on the x-y-z orientation , which indices must be picked
            uint tPermutationOrder = this->permutation_order( tPairCount );

            // obtain the integration mesh
            moris::mtk::Integration_Mesh *tIntegrationMesh = mMeshManager->get_integration_mesh( mMeshIndex );

            // prerequisite offset vector and sets on each periodic surface
            moris::Matrix< DDRMat >    tOffsetVector;
            Vector< std::string > tFirstSideSetNames;
            Vector< std::string > tSecondSideSetNames;

            // find the offset of two surfaces
            this->offset_vector( tOffsetVector, tFirstSideSetNames, tSecondSideSetNames, tPairCount );

            // initialize all the clusters on each side
            Vector< mtk::Cluster const * > tSideClusters1;
            Vector< mtk::Cluster const * > tSideClusters2;

            // cells on the side set are locally indexed and they are assigned a phase to each one
            Vector< moris_index > tLocalCellNumToColor1;

            for ( uint i = 0; i < tFirstSideSetNames.size(); i++ )
            {
                // get the set on the first side
                moris::mtk::Set *tSet = tIntegrationMesh->get_set_by_name( tFirstSideSetNames( i ) );

                // get color(phase) associated with the set
                Matrix< IndexMat > tColor = tSet->get_set_colors();

                // get clusters in the second set
                Vector< moris::mtk::Cluster const * > tSetClusters = tSet->get_clusters_on_set();

                // Assign phases to cluster
                tLocalCellNumToColor1.resize( tLocalCellNumToColor1.size() + tSetClusters.size(), tColor( 0 ) );

                // append the clusters the list of the first side clusters
                tSideClusters1.append( tSetClusters );
            }

            // cells on the side set are locally indexed and they are assigned a phase to each one
            Vector< moris_index > tLocalCellNumToColor2;

            for ( uint i = 0; i < tSecondSideSetNames.size(); i++ )
            {
                // get the set on the first side
                moris::mtk::Set *tSet = tIntegrationMesh->get_set_by_name( tSecondSideSetNames( i ) );

                // get clusters in the second set
                Vector< moris::mtk::Cluster const * > tSetClusters = tSet->get_clusters_on_set();

                // get color(phase) associated with the set
                Matrix< IndexMat > tColor = tSet->get_set_colors();

                // Assign phases to clusters
                tLocalCellNumToColor2.resize( tLocalCellNumToColor2.size() + tSetClusters.size(), tColor( 0 ) );

                // append the clusters the list of the first side clusters
                tSideClusters2.append( tSetClusters );
            }

            // loop over the background cells in order to cut the matching pairs
            for ( const auto &tIter : mBackgroundCellToSideClusterMap1( tPairCount ) )
            {
                // Find all the side clusters that lie in the left and right pair
                Vector< moris_index > const &tSideClustersIndices1 = tIter.second;
                Vector< moris_index > const &tSideClustersIndices2 = mBackgroundCellToSideClusterMap2( tPairCount )[ tIter.first ];

                // Consider different cases to minimize the effort of polygon clipping
                // case 1 where both sides just have one cluster ( they can be trivial or non trivial)
                if ( tSideClustersIndices1.size() == 1 and tSideClustersIndices2.size() == 1 )
                {
                    // parametric coordinates are predefined as the corners and center
                    Vector< Matrix< DDRMat > > tParamCoordsCell2( 1, Matrix< DDRMat >( 1, 2 ) );
                    tParamCoordsCell2( 0 ) = { { -1.0, +1.0 } };

                    // extract the phase of each side cluster
                    moris_index tPhase2 = tLocalCellNumToColor2( tSideClustersIndices2( 0 ) );
                    moris_index tPhase1 = tLocalCellNumToColor1( tSideClustersIndices1( 0 ) );

                    // phase interaction table value
                    moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                    // create a dummy cell in order to specify that all of the cut cells need to be assigned
                    Vector< moris_index > tDummyCel = { 0 };

                    // create the dbl sided cluster based on the cut cell and IP cell information
                    this->create_dbl_sided_cluster(
                            tParamCoordsCell2,
                            tDummyCel,
                            tSideClusters1( tSideClustersIndices1( 0 ) )->get_interpolation_cell(),
                            tSideClusters2( tSideClustersIndices2( 0 ) )->get_interpolation_cell(),
                            tPairCount,
                            tPhaseToPhaseIndex );
                }

                // case 2 where left side has one cluster(trivial/non-trivial) and right has more than one
                else if ( tSideClustersIndices1.size() == 1 and tSideClustersIndices2.size() > 1 )
                {
                    // iterate through the clusters on the right in order to create dbl sided cluster separately
                    for ( uint iCluster = 0; iCluster < tSideClustersIndices2.size(); iCluster++ )
                    {
                        // obtain number of surfaces present in each side cluster
                        uint tPrimaryCellInClusterNum = tSideClusters2( tSideClustersIndices2( iCluster ) )->get_num_primary_cells();

                        // Initialize the cell containing coordinates of the cut surfaces
                        Vector< moris::Matrix< DDRMat > > tParamCoordsCell( tPrimaryCellInClusterNum, moris::Matrix< DDRMat >( 2, 3 ) );

                        // Iterate through each surface to get the coordinates on the side ordinal
                        for ( uint iCell = 0; iCell < tPrimaryCellInClusterNum; iCell++ )
                        {
                            // get the coordinates on the side ordinal and transpose
                            moris::Matrix< DDRMat > t3DcoordMatrix = trans( tSideClusters2( tSideClustersIndices2( iCluster ) )->get_cell_local_coords_on_side_wrt_interp_cell( iCell ) );

                            // convert 3D coordinates to surface coordinates
                            moris::Matrix< DDRMat > tSurfaceCoordMatrix( 1, 2 );
                            tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder );

                            // sort the matrix if the order is not ascending
                            if ( tSurfaceCoordMatrix( 1 ) < tSurfaceCoordMatrix( 0 ) )
                            {
                                tSurfaceCoordMatrix = { { tSurfaceCoordMatrix( 1 ), tSurfaceCoordMatrix( 0 ) } };
                            }

                            // add the facet to the list of coordinates
                            tParamCoordsCell( iCell ) = tSurfaceCoordMatrix;
                        }

                        // create a dummy cell to indicate all the triangle elements will be processed to create dbl sided set
                        Vector< moris_index > tDummyCell( tPrimaryCellInClusterNum );

                        // fill in the dummy cell with dummy values to include all the tri elements
                        for ( uint i = 0; i < tParamCoordsCell.size(); i++ )
                        {
                            tDummyCell( i ) = i;
                        }

                        // get the phase information of the clusters
                        moris_index tPhase2 = tLocalCellNumToColor2( tSideClustersIndices2( iCluster ) );
                        moris_index tPhase1 = tLocalCellNumToColor1( tSideClustersIndices1( 0 ) );

                        // phase interaction table value
                        moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                        // create the dbl sided cluster
                        this->create_dbl_sided_cluster(
                                tParamCoordsCell,
                                tDummyCell,
                                tSideClusters1( tSideClustersIndices1( 0 ) )->get_interpolation_cell(),
                                tSideClusters2( tSideClustersIndices2( iCluster ) )->get_interpolation_cell(),
                                tPairCount,
                                tPhaseToPhaseIndex );
                    }
                }

                // case 3 where right side has one cluster(trivial/non-trivial) and left has more than one
                else if ( tSideClustersIndices1.size() > 1 and tSideClustersIndices2.size() == 1 )
                {
                    // iterate through the clusters on the right in order to create dbl sided cluster separately
                    for ( uint iCluster = 0; iCluster < tSideClustersIndices1.size(); iCluster++ )
                    {
                        // obtain number of surfaces present in each side cluster
                        uint tPrimaryCellInClusterNum = tSideClusters1( tSideClustersIndices1( iCluster ) )->get_num_primary_cells();

                        // Initialize the cell containing coordinates of the cut surfaces
                        Vector< moris::Matrix< DDRMat > > tParamCoordsCell( tPrimaryCellInClusterNum, moris::Matrix< DDRMat >( 2, 3 ) );

                        // Iterate through each surface to get the coordinates on the side ordinal
                        for ( uint iCell = 0; iCell < tPrimaryCellInClusterNum; iCell++ )
                        {
                            // get the coordinates on the side ordinal and transpose
                            moris::Matrix< DDRMat > t3DcoordMatrix = trans( tSideClusters1( tSideClustersIndices1( iCluster ) )->get_cell_local_coords_on_side_wrt_interp_cell( iCell ) );

                            // convert 3D coordinates to surface coordinates
                            moris::Matrix< DDRMat > tSurfaceCoordMatrix( 1, 2 );
                            tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder );

                            // sort the matrix if the order is not ascending
                            if ( tSurfaceCoordMatrix( 1 ) < tSurfaceCoordMatrix( 0 ) )
                            {
                                tSurfaceCoordMatrix = { { tSurfaceCoordMatrix( 1 ), tSurfaceCoordMatrix( 0 ) } };
                            }

                            // add the surface triangle to the list of coordinates
                            tParamCoordsCell( iCell ) = tSurfaceCoordMatrix;
                        }

                        // create a dummy cell to indicate all the triangle elements will be processed to create dbl sided set
                        Vector< moris_index > tDummyCell( tPrimaryCellInClusterNum );

                        // fill in the dummy cell with dummy values to include all the tri elements
                        for ( uint i = 0; i < tParamCoordsCell.size(); i++ )
                        {
                            tDummyCell( i ) = i;
                        }

                        // get the phase information of the each side cluster on the sides
                        moris_index tPhase2 = tLocalCellNumToColor2( tSideClustersIndices2( 0 ) );
                        moris_index tPhase1 = tLocalCellNumToColor1( tSideClustersIndices1( iCluster ) );

                        // phase interaction table value
                        moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                        // iterate through the clusters on the right in order to create dbl sided cluster separately
                        this->create_dbl_sided_cluster(
                                tParamCoordsCell,
                                tDummyCell,
                                tSideClusters1( tSideClustersIndices1( iCluster ) )->get_interpolation_cell(),
                                tSideClusters2( tSideClustersIndices2( 0 ) )->get_interpolation_cell(),
                                tPairCount,
                                tPhaseToPhaseIndex );
                    }
                }

                // case 4 where both sides have more than one cluster within
                else
                {
                    uint tIGCellNum = 0;
                    // loop through the clusters to obtain number of IG surfaces of each cluster
                    for ( uint iCluster = 0; iCluster < tSideClustersIndices1.size(); iCluster++ )
                    {
                        tIGCellNum += tSideClusters1( tSideClustersIndices1( iCluster ) )->get_num_primary_cells();
                    }

                    // Initialize the cut cell surfaces
                    Vector< moris::Matrix< DDRMat > > tParamCoordsCell1( tIGCellNum, moris::Matrix< DDRMat >( 1, 2 ) );

                    // A local map determining that each IG cell/surfaces belongs to which local cluster
                    moris::Matrix< moris::IndexMat > tIGCellToSideClusterMap1( tIGCellNum, 1 );

                    // iterative counter to assign the cut cell surfaces and IG cell to local cluster index map
                    uint iCounter = 0;

                    // loop through each cluster to assign coordinates of the surface
                    for ( uint iCluster = 0; iCluster < tSideClustersIndices1.size(); iCluster++ )
                    {
                        // number of surfaces in each cluster
                        uint tPrimaryCellInClusterNum = tSideClusters1( tSideClustersIndices1( iCluster ) )->get_num_primary_cells();

                        // loop through IG cells in the cluster to get their coordinates
                        for ( uint iCell = 0; iCell < tPrimaryCellInClusterNum; iCell++ )
                        {
                            // get the coordinates on the side ordinal and transpose
                            moris::Matrix< DDRMat > t3DcoordMatrix = trans( tSideClusters1( tSideClustersIndices1( iCluster ) )->get_cell_local_coords_on_side_wrt_interp_cell( iCell ) );

                            // convert 3D coordinates to surface coordinates
                            moris::Matrix< DDRMat > tSurfaceCoordMatrix( 1, 2 );
                            tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder );

                            // sort the matrix if the order is not ascending
                            if ( tSurfaceCoordMatrix( 1 ) < tSurfaceCoordMatrix( 0 ) )
                            {
                                tSurfaceCoordMatrix = { { tSurfaceCoordMatrix( 1 ), tSurfaceCoordMatrix( 0 ) } };
                            }

                            // add the surface triangle to the list of coordinates
                            tParamCoordsCell1( iCounter ) = tSurfaceCoordMatrix;

                            // add the local side cluster index of the IG cell which it belongs to
                            tIGCellToSideClusterMap1( iCounter ) = iCluster;

                            // increase the counter by 1
                            iCounter++;
                        }
                    }

                    // initialize total number of IG cells/surfaces on side 2
                    tIGCellNum = 0;

                    // loop through the clusters to obtain number of IG surfaces of each cluster
                    for ( uint iCluster = 0; iCluster < tSideClustersIndices2.size(); iCluster++ )
                    {
                        tIGCellNum += tSideClusters2( tSideClustersIndices2( iCluster ) )->get_num_primary_cells();
                    }

                    // Initialize the cell of surfaces that will be intersected
                    Vector< moris::Matrix< DDRMat > > tParamCoordsCell2( tIGCellNum, moris::Matrix< DDRMat >( 1, 2 ) );

                    // A local map determining that each IG cell/surfaces belongs to which local cluster
                    moris::Matrix< moris::IndexMat > tIGCellToSideClusterMap2( tIGCellNum, 1 );

                    // iterative counter to assign the cut cell surfaces and IG cell to local cluster index map
                    iCounter = 0;

                    // loop through each cluster to assign coordinates of the surface
                    for ( uint iCluster = 0; iCluster < tSideClustersIndices2.size(); iCluster++ )
                    {
                        // number of surfaces in each cluster
                        uint tPrimaryCellInClusterNum = tSideClusters2( tSideClustersIndices2( iCluster ) )->get_num_primary_cells();

                        // loop through IG cells in the cluster to get their coordinates
                        for ( uint iCell = 0; iCell < tPrimaryCellInClusterNum; iCell++ )
                        {
                            // get the coordinates on the side ordinal and transpose
                            moris::Matrix< DDRMat > t3DcoordMatrix = trans( tSideClusters2( tSideClustersIndices2( iCluster ) )->get_cell_local_coords_on_side_wrt_interp_cell( iCell ) );

                            // convert 3D coordinates to surface coordinates
                            moris::Matrix< DDRMat > tSurfaceCoordMatrix( 1, 2 );
                            tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder );

                            // sort the matrix if the order is not ascending
                            if ( tSurfaceCoordMatrix( 1 ) < tSurfaceCoordMatrix( 0 ) )
                            {
                                tSurfaceCoordMatrix = { { tSurfaceCoordMatrix( 1 ), tSurfaceCoordMatrix( 0 ) } };
                            }

                            // add the surface triangle to the list of coordinates
                            tParamCoordsCell2( iCounter ) = tSurfaceCoordMatrix;

                            // add the local side cluster index of the IG cell which it belongs to
                            tIGCellToSideClusterMap2( iCounter ) = iCluster;

                            // increase the counter by 1
                            iCounter++;
                        }
                    }

                    // initialize the cut polygons and their respective indices that come from their parent element
                    Vector< moris::Matrix< DDRMat > > tCutPolygons;
                    moris::Matrix< moris::IndexMat >       tCutPolygonIdentifier;

                    // Polygon clipping algorithm
                    this->elementwise_bruteforce_search(
                            tParamCoordsCell1, tIGCellToSideClusterMap1, tParamCoordsCell2, tIGCellToSideClusterMap2, tCutPolygons, tCutPolygonIdentifier );

                    // a map from the identifier of each cut cell to all cut cells with the same identifier
                    std::unordered_map< moris_index, Vector< moris_index > > tCutCellIdentifierToCutCell;

                    // populate the map
                    this->group_cut_cells( tCutPolygonIdentifier, tCutCellIdentifierToCutCell );

                    // unique multiplier to distinguish between sub clusters
                    uint tMultiplier = std::max( tParamCoordsCell1.size(), tParamCoordsCell2.size() );

                    // iterate through grouped cut cell in order to create dbl sided cluster with the correct
                    // IP cell and phase information
                    for ( auto &iCutCellGroup : tCutCellIdentifierToCutCell )
                    {
                        // indices of the cut polygon that exist in the particular dbl sided cluster
                        Vector< moris_index > &tIndicesInCutCell = iCutCellGroup.second;

                        // local index of the side cluster on the right
                        uint j = ( iCutCellGroup.first ) % tMultiplier;

                        // local index of the side cluster on the right
                        uint i = std::floor( iCutCellGroup.first / tMultiplier );

                        // phase information of the side clusters
                        moris_index tPhase2 = tLocalCellNumToColor2( tSideClustersIndices2( j ) );
                        moris_index tPhase1 = tLocalCellNumToColor1( tSideClustersIndices1( i ) );

                        // phase interaction table value
                        moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                        // create the dbl sided cluster
                        this->create_dbl_sided_cluster(
                                tCutPolygons,
                                tIndicesInCutCell,
                                tSideClusters1( tSideClustersIndices1( i ) )->get_interpolation_cell(),
                                tSideClusters2( tSideClustersIndices2( j ) )->get_interpolation_cell(),
                                tPairCount,
                                tPhaseToPhaseIndex );
                    }
                    // Find the intersection of two surfaces along with the relative IP cells
                }
            }
        } /*end loop foe the pairs */

        // construct the double sided set and add it to the mesh
        this->construct_add_dbl_sided_set( tPhaseInteractionTable );
    }

    //------------------------------------------------------------------------

    void
    Periodic_2D::create_dbl_sided_cluster(
            Vector< Matrix< DDRMat > > tP,
            Vector< moris_index >     &aIndicesInCutCell,
            moris::mtk::Cell const         &aInterpCell1,
            moris::mtk::Cell const         &aInterpCell2,
            uint                            aPairCount,
            moris_index                     aPhaseToPhase )
    {
        // convert to desired form
        moris::Matrix< DDRMat > tRotation;
        moris::Matrix< DDRMat > tInverseRotation;
        this->rotation_matrix( tRotation, tInverseRotation, aPairCount );

        // make a matrix of intersection points
        Matrix< DDRMat > tUniqueIntersectedPoints;

        // populate the matrix
        for ( size_t i = 0; i < aIndicesInCutCell.size(); i++ )
        {
            tUniqueIntersectedPoints = join_horiz( tUniqueIntersectedPoints, tP( aIndicesInCutCell( i ) ) );
        }

        // create leader and follower side vertices
        this->create_leader_vertices( tUniqueIntersectedPoints, aInterpCell1, aPairCount );
        this->create_follower_vertices( tUniqueIntersectedPoints, aInterpCell2, aPairCount );

        // the outer cell is the side cluster number
        mSideClusterToCells( mNumSideClusters ).resize( tP.size() );
        mSideClusterToCells( mNumSideClusters + 1 ).resize( tP.size() );

        std::iota( mSideClusterToCells( mNumSideClusters ).begin(),
                mSideClusterToCells( mNumSideClusters ).end(),
                mNumCells );

        std::iota( mSideClusterToCells( mNumSideClusters + 1 ).begin(),
                mSideClusterToCells( mNumSideClusters + 1 ).end(),
                mNumCells + tP.size() );

        // loop over the side clusters needs to be created( It is 1 in 2D)
        for ( size_t tClusterNum = 0; tClusterNum < tP.size(); tClusterNum++ )
        {
            // create leader IG cell
            this->create_leader_ig_cell( aInterpCell1, aPairCount );

            // create follower IG cell
            this->create_follower_ig_cell( aInterpCell2, aPairCount );
        }

        // create a matrix for ordinals for follower and leader cells
        moris::Matrix< moris::IndexMat > tLeaderIntegrationCellSideOrdinals( 1, tP.size(), 0 );
        moris::Matrix< moris::IndexMat > tFollowerIntegrationCellSideOrdinals( 1, tP.size(), 0 );

        // parametric coords of leader and follower vertices
        moris::Matrix< DDRMat > tLeaderParamCoords( tUniqueIntersectedPoints.n_cols(), 2 );
        moris::Matrix< DDRMat > tFollowerParamCoords( tUniqueIntersectedPoints.n_cols(), 2 );

        // loop through the
        for ( uint i = 0; i < tUniqueIntersectedPoints.n_cols(); i++ )
        {
            moris::Matrix< DDRMat > tLeaderParamCoord   = { { -1, tUniqueIntersectedPoints( 0, i ) } };
            moris::Matrix< DDRMat > tFollowerParamCoord = { { +1, tUniqueIntersectedPoints( 0, i ) } };
            tLeaderParamCoord                           = tLeaderParamCoord * tRotation;
            tFollowerParamCoord                         = tFollowerParamCoord * tRotation;
            tLeaderParamCoords.get_row( i )             = tLeaderParamCoord.get_row( 0 );
            tFollowerParamCoords.get_row( i )           = tFollowerParamCoord.get_row( 0 );
        }

        // put the coordinate matrix in the right spot
        mVertexParametricCoords( { mNumParamCoords, mNumParamCoords + tUniqueIntersectedPoints.n_cols() - 1 }, { 0, 1 } )                                         = tLeaderParamCoords.matrix_data();
        mVertexParametricCoords( { mNumParamCoords + tUniqueIntersectedPoints.n_cols(), mNumParamCoords + 2 * tUniqueIntersectedPoints.n_cols() - 1 }, { 0, 1 } ) = tFollowerParamCoords.matrix_data();

        // increment number of double sided clusters
        mNumParamCoords += 2 * tUniqueIntersectedPoints.n_cols();

        // save the double sided cluster phase index
        mDoubleSidedClustersIndex( mNumDblSideCluster ) = aPhaseToPhase;

        // increase the count of double sided cluster by 2, 1 for leader and 1 for follower
        mNumSideClusters += 2;

        // increase the count of double sided cluster by 2, 1 for leader and 1 for follower
        mNumCells += 2;

        // increase the count of double sided cluster
        mNumDblSideCluster++;
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Periodic_2D::create_leader_ig_cell( moris::mtk::Cell const &aLeaderInterpCell, uint aPairCount )
    {
        // define the permutation of the vertices
        Vector< moris_index > tPermutation;

        if ( mMeshSideSetPairs( aPairCount )( 0 ) == "4" )
        {
            tPermutation = { { 1, 0, 2 } };
        }
        else
        {
            tPermutation = { { 0, 1, 2 } };
        }

        // we will use all the vertices so reserve size as much of as indices
        mCellToVertexIndices( mNumSideClusters ).resize( mSideClusterToVertexIndices( mNumSideClusters ).size() );

        // use the vertex indices in the cluster to assign for the cell with the order descried above
        std::transform( tPermutation.begin(),
                tPermutation.end(),
                mCellToVertexIndices( mNumSideClusters ).begin(),
                [ this ]( moris_index aIndex ) { return mSideClusterToVertexIndices( mNumSideClusters )( aIndex ); } );

        // store the IP cell( cell cluster ) corresponding to the side cluster
        mSideClusterToIPCell( mNumSideClusters ) = aLeaderInterpCell.get_index();
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Periodic_2D::create_follower_ig_cell( moris::mtk::Cell const &aFollowerInterpCell, uint aPairCount )
    {
        // define the permutation of the vertices
        Vector< moris_index > tPermutation;

        if ( mMeshSideSetPairs( aPairCount )( 0 ) == "4" )
        {
            tPermutation = { 0, 1, 2 };
        }
        else
        {
            tPermutation = { 1, 0, 2 };
        }

        // we will use all the vertices so reserve size as much of as indices
        mCellToVertexIndices( mNumSideClusters + 1 ).resize( mSideClusterToVertexIndices( mNumSideClusters + 1 ).size() );

        // use the vertex indices in the cluster to assign for the cell with the order descried above
        std::transform( tPermutation.begin(),
                tPermutation.end(),
                mCellToVertexIndices( mNumSideClusters + 1 ).begin(),
                [ this ]( moris_index aIndex ) { return mSideClusterToVertexIndices( mNumSideClusters + 1 )( aIndex ); } );

        // store the IP cell( cell cluster ) corresponding to the side cluster
        mSideClusterToIPCell( mNumSideClusters + 1 ) = aFollowerInterpCell.get_index();
    }

    //------------------------------------------------------------------------------------------------------------
    // name the cluster set
    void
    Periodic_2D::construct_add_dbl_sided_set( moris::Matrix< IndexMat > tPhaseInteractionTable )
    {
        // adjust the size of coordinate matrices
        mVertexParametricCoords.resize( mNumParamCoords, 2 );
        mVerticesCoords.resize( 2, mNumVertices );

        // adjust the sizes of the member data
        mSideClusterToVertexIndices.resize( mNumSideClusters );
        mSideClusterToCells.resize( mNumSideClusters );
        mCellToVertexIndices.resize( mNumSideClusters );
        mSideClusterToIPCell.resize( mNumSideClusters );
        mDoubleSidedClustersIndex.resize( mNumDblSideCluster );
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_2D::offset_vector( moris::Matrix< DDRMat > &tOffsetVector, Vector< std::string > &tFirstSideSetNames, Vector< std::string > &tSecondSideSetNames, uint aPairCount )
    {
        // Integration mesh
        moris::mtk::Integration_Mesh *tIntegrationMesh = mMeshManager->get_integration_mesh( mMeshIndex );

        // get all the sidesets
        Vector< moris::mtk::Side_Set * > const &tListOfSideSets = tIntegrationMesh->get_side_sets();

        Vector< std::string > tSetNames( tListOfSideSets.size() );

        std::transform( tListOfSideSets.cbegin(),
                tListOfSideSets.cend(),
                tSetNames.begin(),
                []( moris::mtk::Set *aSet ) -> std::string { return aSet->get_set_name(); } );

        // loop through all the set names and finds sets attached to input side numbers
        std::string tPairFirstName = "SideSet_";
        tPairFirstName.append( mMeshSideSetPairs( aPairCount )( 0 ) );

        // loop through all the set names and finds sets attached to input side numbers
        std::string tPairSecondName = "SideSet_";
        tPairSecondName.append( mMeshSideSetPairs( aPairCount )( 1 ) );

        // loop over the names to get all of the names
        for ( uint i = 0; i < tSetNames.size(); i++ )
        {
            if ( tSetNames( i ).find( tPairFirstName ) != std::string::npos )
            {
                tFirstSideSetNames.push_back( tSetNames( i ) );
            }
            else if ( tSetNames( i ).find( tPairSecondName ) != std::string::npos )
            {
                tSecondSideSetNames.push_back( tSetNames( i ) );
            }
            else
            {
                continue;
            }
        }

        // initialize all the clusters on each side
        Vector< mtk::Cluster const * > tFirstSideClusters;
        Vector< mtk::Cluster const * > tSecondSideClusters;

        for ( uint i = 0; i < tFirstSideSetNames.size(); i++ )
        {
            moris::mtk::Set *tSet = tIntegrationMesh->get_set_by_name( tFirstSideSetNames( i ) );

            // get clusters in the second set
            Vector< moris::mtk::Cluster const * > tSetClusters = tSet->get_clusters_on_set();

            tFirstSideClusters.append( tSetClusters );
        }

        // populate the second side clusters
        for ( uint i = 0; i < tSecondSideSetNames.size(); i++ )
        {
            moris::mtk::Set *tSet = tIntegrationMesh->get_set_by_name( tSecondSideSetNames( i ) );

            // get clusters in the second set
            Vector< moris::mtk::Cluster const * > tSetClusters = tSet->get_clusters_on_set();

            tSecondSideClusters.append( tSetClusters );
        }

        // Second side info
        // get the Integration cell cluster for the side cluster
        Vector< moris::mtk::Cell const * > const &tSecondCells = tSecondSideClusters( 0 )->get_primary_cells_in_cluster();

        // get the side ordinals of the integration cells in the cluster
        moris::Matrix< moris::IndexMat > tSecondCellOrds = tSecondSideClusters( 0 )->get_cell_side_ordinals();

        // each cluster has only 2 vertices
        Vector< moris::mtk::Vertex const * > tSecondVertices      = tSecondCells( 0 )->get_vertices_on_side_ordinal( tSecondCellOrds( 0 ) );
        moris::Matrix< moris::DDRMat >            tSecondVertexCoords1 = tSecondVertices( 0 )->get_coords();
        moris::Matrix< moris::DDRMat >            tSecondVertexCoords2 = tSecondVertices( 1 )->get_coords();

        // calculate the normal of the right segment,  need something better than std::abs to determine the outward
        moris::real                    tMagnitude    = norm( tSecondVertexCoords2 - tSecondVertexCoords1 );
        moris::Matrix< moris::DDRMat > tSecondNormal = { { std::abs( tSecondVertexCoords1( 1 ) - tSecondVertexCoords2( 1 ) ), std::abs( tSecondVertexCoords2( 0 ) - tSecondVertexCoords1( 0 ) ) } };
        tSecondNormal                                = tSecondNormal / tMagnitude;

        // First side info
        // get the Integration cell cluster for the side cluster
        Vector< moris::mtk::Cell const * > const &tFirstCells = tFirstSideClusters( 0 )->get_primary_cells_in_cluster();

        // get the side ordinals of the integration cells in the cluster
        moris::Matrix< moris::IndexMat > tFirstCellOrds = tFirstSideClusters( 0 )->get_cell_side_ordinals();

        // get the vertex information
        Vector< moris::mtk::Vertex const * > tFirstVertices = tFirstCells( 0 )->get_vertices_on_side_ordinal( tFirstCellOrds( 0 ) );

        moris::Matrix< moris::DDRMat > tFirstVertexCoords = tFirstVertices( 0 )->get_coords();

        // distance value between two surfaces
        moris::real tOffsetVectorVal = std::abs( dot( tFirstVertexCoords - tSecondVertexCoords2, tSecondNormal ) );

        // construct the offset vector, needs to be modified
        tOffsetVector = tOffsetVectorVal * tSecondNormal;
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_2D::create_leader_vertices( Matrix< DDRMat > tUniqueIntersectedPoints, moris::mtk::Cell const &aLeaderInterpCell, uint aPairCount )
    {
        moris::Matrix< DDRMat > tRotation;
        moris::Matrix< DDRMat > tInverseRotation;
        this->rotation_matrix( tRotation, tInverseRotation, aPairCount );

        uint tNumSurfaceNodes = tUniqueIntersectedPoints.n_cols();

        // populate the indices
        mSideClusterToVertexIndices( mNumSideClusters ).resize( tNumSurfaceNodes + 1 );

        std::iota( mSideClusterToVertexIndices( mNumSideClusters ).begin(),
                mSideClusterToVertexIndices( mNumSideClusters ).end(),
                mNumVertices );

        // coordinates of nodes attached to element to interpolate for physical coordinates
        moris::Matrix< moris::DDRMat > tCoordinates = aLeaderInterpCell.get_vertex_coords();

        // loop over surface nodes
        for ( uint i = 0; i < tNumSurfaceNodes; i++ )
        {
            // form the parametric coordinates of the nodes
            moris::Matrix< DDRMat > tParamCoordsRelativeToElem = { { -1.0, tUniqueIntersectedPoints( 0, i ) } };

            tParamCoordsRelativeToElem = tParamCoordsRelativeToElem * tRotation;

            // obtain physical coordinates of the node
            Matrix< DDRMat > tNewCoords;
            moris::mtk::Interpolation::bilinear_interpolation( tCoordinates, tParamCoordsRelativeToElem, tNewCoords );

            // save coordinates
            mVerticesCoords( 0, mNumVertices ) = tNewCoords( 0 );
            mVerticesCoords( 1, mNumVertices ) = tNewCoords( 1 );

            // increment vertices count
            mNumVertices++;
        }

        // center node coordinates
        moris::Matrix< moris::DDRMat > tNewNodeCoordinates;
        Matrix< DDRMat >               tParamCoordsRelativeToElem = { { 0.0, 0.0 } };
        moris::mtk::Interpolation::bilinear_interpolation( tCoordinates, tParamCoordsRelativeToElem, tNewNodeCoordinates );

        // transfer coordinates
        mVerticesCoords( 0, mNumVertices ) = tNewNodeCoordinates( 0 );
        mVerticesCoords( 1, mNumVertices ) = tNewNodeCoordinates( 1 );

        // increment vertices count
        mNumVertices++;
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_2D::create_follower_vertices( Matrix< DDRMat > tUniqueIntersectedPoints, moris::mtk::Cell const &aFollowerInterpCell, uint aPairCount )
    {

        moris::Matrix< DDRMat > tRotation;
        moris::Matrix< DDRMat > tInverseRotation;
        this->rotation_matrix( tRotation, tInverseRotation, aPairCount );

        uint tNumSurfaceNodes = tUniqueIntersectedPoints.n_cols();

        // populate the indices
        mSideClusterToVertexIndices( mNumSideClusters + 1 ).resize( tNumSurfaceNodes + 1 );

        std::iota( mSideClusterToVertexIndices( mNumSideClusters + 1 ).begin(),
                mSideClusterToVertexIndices( mNumSideClusters + 1 ).end(),
                mNumVertices );

        // coordinates of nodes attached to element to interpolate for physical coordinates
        moris::Matrix< moris::DDRMat > tCoordinates = aFollowerInterpCell.get_vertex_coords();

        for ( uint i = 0; i < tNumSurfaceNodes; i++ )
        {
            // form the parametric coordinates of the nodes
            moris::Matrix< DDRMat > tParamCoordsRelativeToElem = { { +1.0, tUniqueIntersectedPoints( 0, i ) } };

            // covert to the desired form
            tParamCoordsRelativeToElem = tParamCoordsRelativeToElem * tRotation;

            // obtain physical coordinates of the node
            Matrix< DDRMat > tNewCoords;
            moris::mtk::Interpolation::bilinear_interpolation( tCoordinates, tParamCoordsRelativeToElem, tNewCoords );

            // transfer coordinates
            mVerticesCoords( 0, mNumVertices ) = tNewCoords( 0 );
            mVerticesCoords( 1, mNumVertices ) = tNewCoords( 1 );

            // increment vertices count
            mNumVertices++;
        }

        // center node coordinates
        moris::Matrix< moris::DDRMat > tNewNodeCoordinates;
        Matrix< DDRMat >               tParamCoordsRelativeToElem = { { 0.0, 0.0 } };
        moris::mtk::Interpolation::bilinear_interpolation( tCoordinates, tParamCoordsRelativeToElem, tNewNodeCoordinates );

        // update the id and index
        mVerticesCoords( 0, mNumVertices ) = tNewNodeCoordinates( 0 );
        mVerticesCoords( 1, mNumVertices ) = tNewNodeCoordinates( 1 );

        // increment vertices count
        mNumVertices++;
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_2D::rotation_matrix( moris::Matrix< DDRMat > &aRotation, moris::Matrix< DDRMat > &aInverseRotation, uint aPairCount )
    {
        if ( mMeshSideSetPairs( aPairCount )( 0 ) == "4" )
        {
            aRotation        = { { 1, 0 }, { 0, 1 } };
            aInverseRotation = { { 1, 0 }, { 0, 1 } };
        }
        else if ( mMeshSideSetPairs( aPairCount )( 0 ) == "1" )
        {
            aRotation        = { { 0, 1 }, { 1, 0 } };
            aInverseRotation = { { 0, 1 }, { 1, 0 } };
        }
    }

    // ----------------------------------------------------------------------------

    uint
    Periodic_2D::permutation_order( uint const &aPairCount ) const
    {
        uint tPermutationOrder = 0;

        if ( mMeshSideSetPairs( aPairCount )( 0 ) == "4" )
        {
            tPermutationOrder = 1;
        }
        else if ( mMeshSideSetPairs( aPairCount )( 0 ) == "1" )
        {
            tPermutationOrder = 0;
        }
        return tPermutationOrder;
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_2D::generate_identifier(
            Vector< mtk::Cluster const * >                                  &aSideClusters,
            uint                                                                 &aPairCount,
            std::unordered_map< moris::moris_index, Vector< moris_index > > &aBackgroundCellToSideClusterMap ) const
    {
        // reserve memory for output matrix
        uint tOutputSize = aSideClusters.size();

        // this is an approximate size, it is less than this value
        aBackgroundCellToSideClusterMap.reserve( tOutputSize );

        // moris::mtk::Integration_Mesh*  tIntegrationMesh  = mMeshManager->get_integration_mesh( mMeshIndex );

        uint tIndex = this->permutation_order( aPairCount );

        // loop over the clusters
        for ( uint iI = 0; iI < tOutputSize; iI++ )
        {
            // access the ijk postion
            const luint *tIJK = aSideClusters( iI )->get_interpolation_cell().get_ijk();

            // calculate the identifier
            moris_index tIdentifier = tIJK[ tIndex ];

            // add the cluster index to the relevant bg cell
            aBackgroundCellToSideClusterMap[ tIdentifier ].push_back( iI );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_2D::elementwise_bruteforce_search(
            Vector< Matrix< moris::DDRMat > > const &aParamCoordsCell1,
            moris::Matrix< IndexMat > const              &aIGCellToSideClusterMap1,
            Vector< Matrix< moris::DDRMat > > const &aParamCoordsCell2,
            moris::Matrix< IndexMat > const              &aIGCellToSideClusterMap2,
            Vector< Matrix< moris::DDRMat > >       &aIntersectedAreas,
            moris::Matrix< IndexMat >                    &aIntersectedAreasIdentifier ) const
    {
        // multiplier used to assign an id to each cut surfaces based on the parent cells
        uint tMultiplier = std::max( aParamCoordsCell1.size(), aParamCoordsCell2.size() );

        // max size of the cut  polygons
        uint tMaxSize = aParamCoordsCell1.size() * aParamCoordsCell2.size();

        // reserve the space for max
        aIntersectedAreas.reserve( tMaxSize );

        // set the max size in order to avoid resizing at each step
        aIntersectedAreasIdentifier.set_size( tMaxSize, 1 );

        // counter to assign for cut polygon and their identifier
        uint iCounter = 0;

        // Loop over second mesh
        for ( uint iI = 0; iI < aParamCoordsCell1.size(); iI++ )
        {
            // Loop over the first mesh
            for ( uint iJ = 0; iJ < aParamCoordsCell2.size(); iJ++ )
            {
                // initialize matrix as it needs to be refilled
                moris::Matrix< moris::DDRMat > tP;

                // find the intersection of 2 element triangulation
                this->Intersect( aParamCoordsCell1( iI ), aParamCoordsCell2( iJ ), tP );

                // If it is a polygon add to the output
                if ( tP.n_cols() > 1 )
                {
                    // add the cut polygon to the cell and assign identifier to it
                    aIntersectedAreas.push_back( tP );
                    aIntersectedAreasIdentifier( iCounter ) = aIGCellToSideClusterMap1( iI ) * tMultiplier + aIGCellToSideClusterMap2( iJ );

                    // increase the counter by 1
                    iCounter++;
                }
            }
        }

        // resize the cell and matrix to the actual size
        aIntersectedAreas.shrink_to_fit();
        aIntersectedAreasIdentifier.resize( iCounter, 1 );
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_2D::Intersect(
            moris::Matrix< moris::DDRMat > const &aFirstLineCoords,
            moris::Matrix< moris::DDRMat > const &aSecondLineCoords,
            moris::Matrix< moris::DDRMat >       &aIntersectedPoints ) const
    {
        // check if two lines overlap
        if ( ( ( aFirstLineCoords( 1 ) - aSecondLineCoords( 0 ) ) > 0 ) and ( ( aSecondLineCoords( 1 ) - aFirstLineCoords( 0 ) ) > 0 ) )
        {
            // overlap 2 line segments
            aIntersectedPoints = { { std::max( aFirstLineCoords( 0 ), aSecondLineCoords( 0 ) ), std::min( aFirstLineCoords( 1 ), aSecondLineCoords( 1 ) ) } };
        }
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_2D::group_cut_cells(
            moris::Matrix< IndexMat > const                               &aCutTrianglesIdentifier,
            std::unordered_map< moris_index, Vector< moris_index > > &aCutCellIdentifierToCutCellIndex ) const
    {
        // loop over all the cut cells
        for ( uint iI = 0; iI < aCutTrianglesIdentifier.n_rows(); iI++ )
        {
            aCutCellIdentifierToCutCellIndex[ aCutTrianglesIdentifier( iI ) ].push_back( iI );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_2D::determine_max_size()
    {
        // resize the maps based on the spatial dimension
        mBackgroundCellToSideClusterMap1.resize( 2 );
        mBackgroundCellToSideClusterMap2.resize( 2 );

        // initialize the maximum size of the double side clusters
        uint tNumVertices            = 0;
        uint tNumDoubleSidedClusters = 0;

        // loop over the pairs to determine the sizes
        for ( uint tPairCount = 0; tPairCount < mMeshSideSetPairs.size(); tPairCount++ )
        {
            // obtain the integration mesh
            moris::mtk::Integration_Mesh *tIntegrationMesh = mMeshManager->get_integration_mesh( mMeshIndex );

            // prerequisite offset vector and sets on each periodic surface
            moris::Matrix< DDRMat >    tOffsetVector;
            Vector< std::string > tFirstSideSetNames;
            Vector< std::string > tSecondSideSetNames;

            // find the offset of two surfaces
            this->offset_vector( tOffsetVector, tFirstSideSetNames, tSecondSideSetNames, tPairCount );

            // initialize all the clusters on each side
            Vector< mtk::Cluster const * > tSideClusters1;
            Vector< mtk::Cluster const * > tSideClusters2;

            // cells on the side set are locally indexed and they are assigned a phase to each one
            Vector< moris_index > tLocalCellNumToColor1;

            for ( uint i = 0; i < tFirstSideSetNames.size(); i++ )
            {
                // get the set on the first side
                moris::mtk::Set *tSet = tIntegrationMesh->get_set_by_name( tFirstSideSetNames( i ) );

                // get color(phase) associated with the set
                Matrix< IndexMat > tColor = tSet->get_set_colors();

                // get clusters in the second set
                Vector< moris::mtk::Cluster const * > tSetClusters = tSet->get_clusters_on_set();

                // append the clusters the list of the first side clusters
                tSideClusters1.append( tSetClusters );
            }

            // cells on the side set are locally indexed and they are assigned a phase to each one
            Vector< moris_index > tLocalCellNumToColor2;

            for ( uint i = 0; i < tSecondSideSetNames.size(); i++ )
            {
                // get the set on the first side
                moris::mtk::Set *tSet = tIntegrationMesh->get_set_by_name( tSecondSideSetNames( i ) );

                // get clusters in the second set
                Vector< moris::mtk::Cluster const * > tSetClusters = tSet->get_clusters_on_set();

                // append the clusters the list of the first side clusters
                tSideClusters2.append( tSetClusters );
            }

            // populate the map and identifier
            this->generate_identifier( tSideClusters1, tPairCount, mBackgroundCellToSideClusterMap1( tPairCount ) );
            this->generate_identifier( tSideClusters2, tPairCount, mBackgroundCellToSideClusterMap2( tPairCount ) );

            for ( const auto &tKeyValue : mBackgroundCellToSideClusterMap1( tPairCount ) )
            {
                auto iIJK            = tKeyValue.first;
                auto iClusterIndices = tKeyValue.second;

                tNumDoubleSidedClusters += iClusterIndices.size() * mBackgroundCellToSideClusterMap2( tPairCount )[ iIJK ].size();
            }

            uint tNumVertices1 = 0;
            for ( const auto &iCluster : tSideClusters1 )
                tNumVertices1 += iCluster->get_num_vertices_in_cluster() + 1;

            uint tNumVertices2 = 0;
            for ( const auto &iCluster : tSideClusters2 )
                tNumVertices2 += iCluster->get_num_vertices_in_cluster() + 1;

            tNumVertices += ( 2 * std::max( tNumVertices2, tNumVertices1 ) );
        }

        mSideClusterToVertexIndices.resize( 2 * tNumDoubleSidedClusters );

        mSideClusterToCells.resize( 2 * tNumDoubleSidedClusters );

        mSideClusterToIPCell.resize( 2 * tNumDoubleSidedClusters );

        mCellToVertexIndices.resize( 2 * tNumDoubleSidedClusters );

        mDoubleSidedClustersIndex.resize( tNumDoubleSidedClusters );

        mVerticesCoords.set_size( 2, tNumVertices );

        mVertexParametricCoords.set_size( tNumVertices, 2 );
    }
}    // namespace moris::mig
