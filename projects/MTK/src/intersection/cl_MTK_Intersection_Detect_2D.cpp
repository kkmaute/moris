/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Intersection_Detect_2D.cpp
 *
 */

#include "cl_MTK_Intersection_Detect_2D.hpp"
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

namespace moris
{
    namespace mtk
    {

        Intersection_Detect_2D::Intersection_Detect_2D(
                std::shared_ptr< moris::mtk::Mesh_Manager > aMeshManager,
                moris::moris_index                          aMeshIndex,
                moris::Parameter_List                      &aParameterList,
                moris::uint                                 aNumBulkPhases )
                : mMeshManager( aMeshManager )
                , mMeshIndex( aMeshIndex )
                , mNumBulkPhases( aNumBulkPhases )
        {
            // get the periodic mesh set names
            std::string tMeshSideSetNames = aParameterList.get< std::string >( "periodic_side_set_pair" );

            // store the names in mMeshSideSetPairs
            string_to_cell_of_cell( aParameterList.get< std::string >( "periodic_side_set_pair" ), mMeshSideSetPairs );

            // initialize the mesh intersection data to store the id and index of the added entities
            mIntersectedMeshData.set_up_external_entity_data( mMeshManager->get_integration_mesh( mMeshIndex ) );
        }

        //---------------------------------------------------------------------

        Intersection_Detect_2D::~Intersection_Detect_2D()
        {
            // delete all the pointer data
            for ( auto p : mDoubleSidedClusters )
            {
                delete p;
            }
            for ( auto p : mLeaderSidedClusters )
            {
                delete p;
            }
            for ( auto p : mFollowerSidedClusters )
            {
                delete p;
            }
            for ( auto p : mLeaderSideCells )
            {
                delete p;
            }
            for ( auto p : mFollowerSideCells )
            {
                delete p;
            }
            for ( auto p : mLeaderVertices )
            {
                delete p;
            }
            for ( auto p : mFollowerVertices )
            {
                delete p;
            }
        }

        //---------------------------------------------------------------------

        void
        Intersection_Detect_2D::perform()
        {
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

                // Initialize the background cell to side cluster map
                std::unordered_map< moris::moris_index, Vector< moris_index > > tBackgroundCellToSideClusterMap1;
                std::unordered_map< moris::moris_index, Vector< moris_index > > tBackgroundCellToSideClusterMap2;

                // populate the map and identifier
                this->generate_identifier( tSideClusters1, tPairCount, tBackgroundCellToSideClusterMap1 );
                this->generate_identifier( tSideClusters2, tPairCount, tBackgroundCellToSideClusterMap2 );

                // loop over the background cells in order to cut the matching pairs
                for ( const auto &tIter : tBackgroundCellToSideClusterMap1 )
                {
                    // Find all the side clusters that lie in the left and right pair
                    Vector< moris_index > const &tSideClustersIndices1 = tIter.second;
                    Vector< moris_index > const &tSideClustersIndices2 = tBackgroundCellToSideClusterMap2[ tIter.first ];

                    // Consider different cases to minmize the effort of polygon clipping
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
                        this->create_dbl_sided_cluster( tParamCoordsCell2, tDummyCel, tSideClusters1( tSideClustersIndices1( 0 ) )->get_interpolation_cell(), tSideClusters2( tSideClustersIndices2( 0 ) )->get_interpolation_cell(), tPairCount, tPhaseToPhaseIndex );
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
                            this->create_dbl_sided_cluster( tParamCoordsCell, tDummyCell, tSideClusters1( tSideClustersIndices1( 0 ) )->get_interpolation_cell(), tSideClusters2( tSideClustersIndices2( iCluster ) )->get_interpolation_cell(), tPairCount, tPhaseToPhaseIndex );
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
                            this->create_dbl_sided_cluster( tParamCoordsCell, tDummyCell, tSideClusters1( tSideClustersIndices1( iCluster ) )->get_interpolation_cell(), tSideClusters2( tSideClustersIndices2( 0 ) )->get_interpolation_cell(), tPairCount, tPhaseToPhaseIndex );
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

                        // iterate through grouped cut cell in order to create dbl sideded cluster with the correct
                        // IP cell and phase information
                        for ( auto &iCutCellGroup : tCutCellIdentifierToCutCell )
                        {
                            // indices of the cut polygon that exist in the particular dbl sided cluster
                            Vector< moris_index > &tIndicesinCutCell = iCutCellGroup.second;

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
                            this->create_dbl_sided_cluster( tCutPolygons, tIndicesinCutCell, tSideClusters1( tSideClustersIndices1( i ) )->get_interpolation_cell(), tSideClusters2( tSideClustersIndices2( j ) )->get_interpolation_cell(), tPairCount, tPhaseToPhaseIndex );
                        }
                        // Find the intersection of two surfaces along with the relative IP cells
                    }
                }
            } /*end loop foe the pairs */

            // construct the double sided set and add it to the mesh
            this->constrcuct_add_dbl_sided_set( tPhaseInteractionTable );
        }

        //------------------------------------------------------------------------

        void
        Intersection_Detect_2D::create_dbl_sided_cluster( Vector< Matrix< DDRMat > > tP,
                Vector< moris_index >                                               &aIndicesinCutCell,
                moris::mtk::Cell const                                                   &aInterpCell1,
                moris::mtk::Cell const                                                   &aInterpCell2,
                uint                                                                      aPairCount,
                moris_index                                                               aPhaseToPhase )
        {

            // convert to desired form
            moris::Matrix< DDRMat > tRotation;
            moris::Matrix< DDRMat > tInverseRotation;
            this->rotation_matrix( tRotation, tInverseRotation, aPairCount );

            // make a matrix of intersection points
            Matrix< DDRMat > tUniqueIntersectedPoints;

            // populate the matrix
            for ( size_t i = 0; i < aIndicesinCutCell.size(); i++ )
            {
                tUniqueIntersectedPoints = join_horiz( tUniqueIntersectedPoints, tP( aIndicesinCutCell( i ) ) );
            }

            // create leader and follower side vertices
            Vector< moris::mtk::Vertex * > tLeaderVerticesCell   = this->create_leader_vertices( tUniqueIntersectedPoints, aInterpCell1, aPairCount );
            Vector< moris::mtk::Vertex * > tFollowerVerticesCell = this->create_follower_vertices( tUniqueIntersectedPoints, aInterpCell2, aPairCount );

            // initialize leader,follower side cluster data
            Vector< moris::mtk::Cell const * > tLeaderIntegCells;
            Vector< moris::mtk::Cell const * > tFollowerIntegCells;

            // loop over the side clusters needs to be created( It is 1 in 2D)
            for ( size_t tClusterNum = 0; tClusterNum < tP.size(); tClusterNum++ )
            {
                // create leader IG cell
                moris::mtk::Cell const *tLeaderIgCell = this->create_leader_ig_cell( tLeaderVerticesCell, aInterpCell1, aPairCount );
                tLeaderIntegCells.push_back( tLeaderIgCell );

                // create follower IG cell
                moris::mtk::Cell const *tFollowerIgCell = this->create_follower_ig_cell( tFollowerVerticesCell, aInterpCell2, aPairCount );
                tFollowerIntegCells.push_back( tFollowerIgCell );
            }

            // create a matrix for ordinals for follower and leader cells
            moris::Matrix< moris::IndexMat > tLeaderIntegrationCellSideOrdinals( 1, tP.size(), 0 );
            moris::Matrix< moris::IndexMat > tFollowerIntegrationCellSideOrdinals( 1, tP.size(), 0 );

            // parametric coords of leader and follower vertices
            moris::Matrix< DDRMat > tLeaderParamCoords( tUniqueIntersectedPoints.n_cols(), 2 );
            moris::Matrix< DDRMat > tFollowerParamCoords( tUniqueIntersectedPoints.n_cols(), 2 );

            // loop through
            for ( uint i = 0; i < tUniqueIntersectedPoints.n_cols(); i++ )
            {
                moris::Matrix< DDRMat > tLeaderParamCoord   = { { -1, tUniqueIntersectedPoints( 0, i ) } };
                moris::Matrix< DDRMat > tFollowerParamCoord = { { +1, tUniqueIntersectedPoints( 0, i ) } };
                tLeaderParamCoord                           = tLeaderParamCoord * tRotation;
                tFollowerParamCoord                         = tFollowerParamCoord * tRotation;
                tLeaderParamCoords.get_row( i )             = tLeaderParamCoord.get_row( 0 );
                tFollowerParamCoords.get_row( i )           = tFollowerParamCoord.get_row( 0 );
            }

            // create constant version of the vertices on cluster
            Vector< moris::mtk::Vertex const * > tLeaderVerticesConst;
            Vector< moris::mtk::Vertex const * > tFollowerVerticesConst;

            for ( uint i = 0; i < tLeaderVerticesCell.size() - 1; i++ )
            {
                tLeaderVerticesConst.push_back( tLeaderVerticesCell( i ) );
                tFollowerVerticesConst.push_back( tFollowerVerticesCell( i ) );
            }

            moris::mtk::Side_Cluster_ISC *tLeaderSideCluster = new moris::mtk::Side_Cluster_ISC( false,
                    &aInterpCell1,
                    tLeaderIntegCells,
                    tLeaderIntegrationCellSideOrdinals,
                    tLeaderVerticesConst,
                    tLeaderParamCoords );

            moris::mtk::Side_Cluster_ISC *tFollowerSideCluster = new moris::mtk::Side_Cluster_ISC( false,
                    &aInterpCell2,
                    tFollowerIntegCells,
                    tFollowerIntegrationCellSideOrdinals,
                    tFollowerVerticesConst,
                    tFollowerParamCoords );

            // create double side cluster
            mtk::Double_Side_Cluster *tDblSideCluster = new mtk::Double_Side_Cluster(
                    tLeaderSideCluster,
                    tFollowerSideCluster,
                    tFollowerSideCluster->mVerticesInCluster );

            // store double side cluster
            mDoubleSidedClusters.push_back( tDblSideCluster );
            mDoubleSidedClustersIndex.push_back( aPhaseToPhase );

            // store side clusters
            mLeaderSidedClusters.push_back( tLeaderSideCluster );
            mFollowerSidedClusters.push_back( tFollowerSideCluster );

            // store the vertices to prevent memory leak
            mLeaderVertices.append( tLeaderVerticesConst );
            mFollowerVertices.append( tFollowerVerticesConst );

            // Append the integration cells
            mLeaderSideCells.append( tLeaderIntegCells );
            mFollowerSideCells.append( tFollowerIntegCells );
        }

        //------------------------------------------------------------------------------------------------------------

        moris::mtk::Cell const *
        Intersection_Detect_2D::create_leader_ig_cell( Vector< moris::mtk::Vertex * > aLeaderVertices, moris::mtk::Cell const &aLeaderInterpCell, uint aPairCount )
        {
            Matrix< IndexMat > tPermutation;

            if ( mMeshSideSetPairs( aPairCount )( 0 ) == "4" )
            {
                tPermutation = { { 1, 0, 2 } };
            }
            else
            {
                tPermutation = { { 0, 1, 2 } };
            }

            // arrange the permutation based on the set
            Vector< moris::mtk::Vertex * > tLeaderVertices( aLeaderVertices.size() );
            for ( uint i = 0; i < tPermutation.n_cols(); i++ )
            {
                tLeaderVertices( i ) = aLeaderVertices( tPermutation( i ) );
            }

            // allocate id and index from intersection data
            moris_index tCellIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::ELEMENT );
            moris_id    tCellId    = mIntersectedMeshData.allocate_entity_ids_external_entity_data( 1, EntityRank::ELEMENT );

            // cell info for TRI 3
            mtk::Cell_Info_Factory                   tCellInfoFactory;
            std::shared_ptr< moris::mtk::Cell_Info > tLinearCellInfo = tCellInfoFactory.create_cell_info_sp( Geometry_Type::TRI, mtk::Interpolation_Order::LINEAR );

            // create the new integration cell
            moris::mtk::Cell *tIgCell = new Cell_ISC(
                    tCellId,
                    tCellIndex,
                    aLeaderInterpCell.get_owner(),
                    tLinearCellInfo,
                    tLeaderVertices );

            mIntersectedMeshData.update_first_available_index_external_data( tCellIndex + 1, EntityRank::ELEMENT );

            return tIgCell;
        }

        //------------------------------------------------------------------------------------------------------------

        moris::mtk::Cell const *
        Intersection_Detect_2D::create_follower_ig_cell( Vector< moris::mtk::Vertex * > aFollowerVertices, moris::mtk::Cell const &aFollowerInterpCell, uint aPairCount )
        {

            Matrix< IndexMat > tPermutation;

            if ( mMeshSideSetPairs( aPairCount )( 0 ) == "4" )
            {
                tPermutation = { { 0, 1, 2 } };
            }
            else
            {
                tPermutation = { { 1, 0, 2 } };
            }

            // arrange the permutation based on the set
            Vector< moris::mtk::Vertex * > tFollowerVertices( aFollowerVertices.size() );
            for ( uint i = 0; i < tPermutation.n_cols(); i++ )
            {
                tFollowerVertices( i ) = aFollowerVertices( tPermutation( i ) );
            }

            // allocate id and index from the data
            moris_index tCellIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::ELEMENT );
            moris_id    tCellId    = mIntersectedMeshData.allocate_entity_ids_external_entity_data( 1, EntityRank::ELEMENT );

            // cell info  and geometry for TRI 3
            mtk::Cell_Info_Factory                   tCellInfoFactory;
            std::shared_ptr< moris::mtk::Cell_Info > tLinearCellInfo = tCellInfoFactory.create_cell_info_sp( Geometry_Type::TRI, mtk::Interpolation_Order::LINEAR );

            // create the IG cell
            moris::mtk::Cell *tIgCell = new Cell_ISC(
                    tCellId,
                    tCellIndex,
                    aFollowerInterpCell.get_owner(),
                    tLinearCellInfo,
                    tFollowerVertices );

            mIntersectedMeshData.update_first_available_index_external_data( tCellIndex + 1, EntityRank::ELEMENT );

            return tIgCell;
        }

        //------------------------------------------------------------------------------------------------------------
        // name the cluster set
        void
        Intersection_Detect_2D::constrcuct_add_dbl_sided_set( moris::Matrix< IndexMat > tPhaseInteractionTable )
        {
            Tracer tTracer( "MTK", "Double Sided Set", "Construction of Double Sided Sets" );

            // initialize all the double sided sets
            Vector< Vector< mtk::Cluster const * > > tDoubleSideSets( tPhaseInteractionTable.numel() );

            // get the integration mesh
            moris::mtk::Integration_Mesh *tIntegMesh = mMeshManager->get_integration_mesh( mMeshIndex );

            uint tSpatialDim = tIntegMesh->get_spatial_dim();

            // assign double sided clusters to the right double sided set
            for ( uint i = 0; i < mDoubleSidedClusters.size(); i++ )
            {
                tDoubleSideSets( mDoubleSidedClustersIndex( i ) ).push_back( mDoubleSidedClusters( i ) );
            }

            // add double sided set to the mesh and name them
            for ( size_t aBulkPhase1 = 0; aBulkPhase1 < tPhaseInteractionTable.n_cols(); aBulkPhase1++ )
            {
                for ( size_t aBulkPhase2 = 0; aBulkPhase2 < tPhaseInteractionTable.n_cols(); aBulkPhase2++ )
                {
                    uint tPhaseToPhaseIndex = tPhaseInteractionTable( aBulkPhase1, aBulkPhase2 );

                    // check if double sided set is empty
                    if ( tDoubleSideSets( tPhaseToPhaseIndex ).size() > 0 )
                    {
                        // name the double sided set
                        std::string tDoubleSideSetName = "P" + std::to_string( aBulkPhase1 ) + std::to_string( aBulkPhase2 );

                        // colors of the double sided set,  needs to be modified
                        moris::Matrix< moris::IndexMat > tColors = { { 0 } };

                        // Construct the double side set
                        auto *tDblSideSet = new moris::mtk::Double_Side_Set( tDoubleSideSetName, tDoubleSideSets( tPhaseToPhaseIndex ), tColors, tSpatialDim );

                        // add double sided periodic boundary condition to the integration mesh
                        tIntegMesh->add_double_side_set( tDblSideSet );
                    }
                }
            }

            // FIXME: communicate double sided side sets constructed here
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect_2D::offset_vector( moris::Matrix< DDRMat > &tOffsetVector, Vector< std::string > &tFirstSideSetNames, Vector< std::string > &tSecondSideSetNames, uint aPairCount )
        {
            // Integration mesh
            moris::mtk::Integration_Mesh *tIntegrationMesh = mMeshManager->get_integration_mesh( mMeshIndex );

            // get all the sidesets
            Vector< std::string > tSetNames = tIntegrationMesh->get_set_names( EntityRank::EDGE );

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

        Vector< moris::mtk::Vertex * >
        Intersection_Detect_2D::create_leader_vertices( Matrix< DDRMat > tUniqueIntersectedPoints, moris::mtk::Cell const &aLeaderInterpCell, uint aPairCount )
        {
            moris::Matrix< DDRMat > tRotation;
            moris::Matrix< DDRMat > tInverseRotation;
            this->rotation_matrix( tRotation, tInverseRotation, aPairCount );

            uint tNumSurfaceNodes = tUniqueIntersectedPoints.n_cols();

            // get the interpolation mesh
            moris::mtk::Interpolation_Mesh *tInterpMesh = mMeshManager->get_interpolation_mesh( mMeshIndex );

            // initialize the vertices of the cell ( +1 for the head of the tet4 )
            Vector< moris::mtk::Vertex * > tVertcies( tNumSurfaceNodes + 1 );

            // allocate ids for all the vertices
            moris_id tVertexId = mIntersectedMeshData.allocate_entity_ids_external_entity_data( tNumSurfaceNodes + 1, EntityRank::NODE );

            // get the indices of the interpolation cell nodes
            moris::Matrix< moris::IndexMat > tElementNodes =
                    tInterpMesh->get_entity_connected_to_entity_loc_inds(
                            aLeaderInterpCell.get_index(),
                            mtk::EntityRank::ELEMENT,
                            mtk::EntityRank::NODE );

            // coordinates of nodes attached to element to interpolate for physical coordinates
            moris::Matrix< moris::DDRMat > tCoordinates( tElementNodes.numel(), 2 );

            for ( uint i = 0; i < tElementNodes.numel(); i++ )
            {
                moris::Matrix< moris::DDRMat > tmp = tInterpMesh->get_node_coordinate( tElementNodes( i ) );
                tCoordinates.get_row( i )          = tmp.get_row( 0 );
            }

            // loop over surface nodes
            for ( uint i = 0; i < tNumSurfaceNodes; i++ )
            {
                // form the parametric coordinates of the nodes
                moris::Matrix< DDRMat > tParamCoordsRelativeToElem = { { -1.0, tUniqueIntersectedPoints( 0, i ) } };

                tParamCoordsRelativeToElem = tParamCoordsRelativeToElem * tRotation;

                // obtain physical coordinates of the node
                Matrix< DDRMat > tNewCoords;
                moris::mtk::Interpolation::bilinear_interpolation( tCoordinates, tParamCoordsRelativeToElem, tNewCoords );

                moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE );

                // construct the new vertex
                tVertcies( i ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewCoords );

                // update the id and index
                mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex + 1, EntityRank::NODE );
                tVertexId++;
            }

            // center node coordinates
            moris::Matrix< moris::DDRMat > tNewNodeCoordinates;
            Matrix< DDRMat >               tParamCoordsRelativeToElem = { { 0.0, 0.0 } };
            moris::mtk::Interpolation::bilinear_interpolation( tCoordinates, tParamCoordsRelativeToElem, tNewNodeCoordinates );

            moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE );
            mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex + 1, EntityRank::NODE );

            // add the top vertex of TRI#
            tVertcies( tNumSurfaceNodes ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewNodeCoordinates );

            return tVertcies;
        }

        // ----------------------------------------------------------------------------

        Vector< moris::mtk::Vertex * >
        Intersection_Detect_2D::create_follower_vertices( Matrix< DDRMat > tUniqueIntersectedPoints, moris::mtk::Cell const &aFollowerInterpCell, uint aPairCount )
        {

            moris::Matrix< DDRMat > tRotation;
            moris::Matrix< DDRMat > tInverseRotation;
            this->rotation_matrix( tRotation, tInverseRotation, aPairCount );

            uint tNumSurfaceNodes = tUniqueIntersectedPoints.n_cols();

            // get the interpolation mesh
            moris::mtk::Interpolation_Mesh *tInterpMesh = mMeshManager->get_interpolation_mesh( mMeshIndex );

            // initialize the vertices of the cell ( +1 for the head of the tet4 )
            Vector< moris::mtk::Vertex * > tVertcies( tNumSurfaceNodes + 1 );

            // allocate ids for all the vertices
            moris_id tVertexId = mIntersectedMeshData.allocate_entity_ids_external_entity_data( tNumSurfaceNodes + 1, EntityRank::NODE );

            // get the indices of the interpolation cell nodes
            moris::Matrix< moris::IndexMat > tElementNodes =
                    tInterpMesh->get_entity_connected_to_entity_loc_inds(
                            aFollowerInterpCell.get_index(),
                            mtk::EntityRank::ELEMENT,
                            mtk::EntityRank::NODE );

            // coordinates of nodes attached to element to interpolate for physical coordinates
            moris::Matrix< moris::DDRMat > tCoordinates( tElementNodes.numel(), 2 );

            for ( uint i = 0; i < tElementNodes.numel(); i++ )
            {
                moris::Matrix< moris::DDRMat > tmp = tInterpMesh->get_node_coordinate( tElementNodes( i ) );
                tCoordinates.get_row( i )          = tmp.get_row( 0 );
            }

            for ( uint i = 0; i < tNumSurfaceNodes; i++ )
            {
                // form the parametric coordinates of the nodes
                moris::Matrix< DDRMat > tParamCoordsRelativeToElem = { { +1.0, tUniqueIntersectedPoints( 0, i ) } };

                // covert to the desired form
                tParamCoordsRelativeToElem = tParamCoordsRelativeToElem * tRotation;

                // obtain physical coordinates of the node
                Matrix< DDRMat > tNewCoords;
                moris::mtk::Interpolation::bilinear_interpolation( tCoordinates, tParamCoordsRelativeToElem, tNewCoords );

                moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE );

                // construct the new vertex
                tVertcies( i ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewCoords );

                // update the id and index
                mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex + 1, EntityRank::NODE );
                tVertexId++;
            }

            // center node coordinates
            moris::Matrix< moris::DDRMat > tNewNodeCoordinates;
            Matrix< DDRMat >               tParamCoordsRelativeToElem = { { 0.0, 0.0 } };
            moris::mtk::Interpolation::bilinear_interpolation( tCoordinates, tParamCoordsRelativeToElem, tNewNodeCoordinates );

            // update the id and index
            moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE );
            mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex + 1, EntityRank::NODE );

            // add the top vertex of TRI3
            tVertcies( tNumSurfaceNodes ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewNodeCoordinates );

            return tVertcies;
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect_2D::rotation_matrix( moris::Matrix< DDRMat > &aRotation, moris::Matrix< DDRMat > &aInverseRotation, uint aPairCount )
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
        Intersection_Detect_2D::permutation_order( uint const &aPairCount ) const
        {
            uint tPermuationOrder = 0;

            if ( mMeshSideSetPairs( aPairCount )( 0 ) == "4" )
            {
                tPermuationOrder = 1;
            }
            else if ( mMeshSideSetPairs( aPairCount )( 0 ) == "1" )
            {
                tPermuationOrder = 0;
            }
            return tPermuationOrder;
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect_2D::generate_identifier( Vector< mtk::Cluster const * > &aSideClusters,
                uint                                                                     &aPairCount,
                std::unordered_map< moris::moris_index, Vector< moris_index > >     &aBackgroundCellToSideClusterMap ) const
        {
            // reserve memory for output matrix
            uint tOutputSize = aSideClusters.size();

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
        Intersection_Detect_2D::elementwise_bruteforce_search(
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
        Intersection_Detect_2D::Intersect(
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
        Intersection_Detect_2D::group_cut_cells( moris::Matrix< IndexMat > const &aCutTrianglesIdentifier,
                std::unordered_map< moris_index, Vector< moris_index > >    &aCutCellIdentifierToCutCellIndex ) const
        {
            // loop over all the cut cells
            for ( uint iI = 0; iI < aCutTrianglesIdentifier.n_rows(); iI++ )
            {
                aCutCellIdentifierToCutCellIndex[ aCutTrianglesIdentifier( iI ) ].push_back( iI );
            }
        }

        // ----------------------------------------------------------------------------

    } /* end namespace mtk */

} /* end namespace moris */
