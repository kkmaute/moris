/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MIG_Periodic_3D.cpp
 *
 */

#include "cl_MIG_Periodic_3D.hpp"

#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_MTK_Side_Cluster.hpp"
#include "moris_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_unique.hpp"
#include <unordered_map>
#include <utility>

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

#include "cl_Stopwatch.hpp"//CHR/src

namespace moris::mig
{
    Periodic_3D::Periodic_3D(
            std::shared_ptr< moris::mtk::Mesh_Manager > aMeshManager,
            moris::moris_index                          aMeshIndex,
            moris::Parameter_List                      &aParameterList,
            moris::uint                                 aNumBulkPhases )
            : mMeshManager( std::move( aMeshManager ) )
            , mMeshIndex( aMeshIndex )
            , mNumBulkPhases( aNumBulkPhases )
    {
        // get the periodic mesh set names
        std::string tMeshSideSetNames = aParameterList.get< std::string >( "periodic_side_set_pair" );

        // store the names in mMeshSideSetPairs
        string_to_vector_of_vectors( aParameterList.get< std::string >( "periodic_side_set_pair" ), mMeshSideSetPairs );
    }

    //---------------------------------------------------------------------

    Periodic_3D::~Periodic_3D()
    {
    }

    //---------------------------------------------------------------------

    void
    Periodic_3D::perform()
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
            std::pair< moris_id, moris_id > tPermutationOrder = this->permutation_pair( tPairCount );

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
                Vector< moris_index > const &tSideClustersIndices2 = mBackgroundCellToSideClusterMap2( tPairCount )[tIter.first];

                // Consider different cases to minmize the effort of polygon clipping
                // case 1 where both sides just have one cluster ( they can be trivial or non trivial)
                if ( tSideClustersIndices1.size() == 1 and tSideClustersIndices2.size() == 1 )
                {
                    // parametric coordinates are predefined as the corners and center
                    Vector< Matrix< DDRMat > > tParamCoordsCell2( 4, Matrix< DDRMat >( 2, 3 ) );
                    tParamCoordsCell2( 0 ) = { { -1.0, +1.0, 0.0 }, { -1.0, -1.0, 0.0 } };
                    tParamCoordsCell2( 1 ) = { { -1.0, 0.0, -1.0 }, { -1.0, 0.0, +1.0 } };
                    tParamCoordsCell2( 2 ) = { { -1.0, 0.0, 1.0 }, { +1.0, 0.0, 1.0 } };
                    tParamCoordsCell2( 3 ) = { { 1.0, 0.0, +1.0 }, { 1.0, 0.0, -1.0 } };

                    // extract the phase of each side cluster
                    moris_index tPhase2 = tLocalCellNumToColor2( tSideClustersIndices2( 0 ) );
                    moris_index tPhase1 = tLocalCellNumToColor1( tSideClustersIndices1( 0 ) );

                    // phase interaction table value
                    moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                    // create a dummy cell in order to specify that all of the cut cells need to be assigned
                    Vector< moris_index > tDummyCel = { 0, 1, 2, 3 };

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
                            moris::Matrix< DDRMat > tSurfaceCoordMatrix( 2, 3 );

                            tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder.first );
                            tSurfaceCoordMatrix.get_row( 1 ) = t3DcoordMatrix.get_row( tPermutationOrder.second );

                            // sort them counter clock wise in order to be treated the same as other coordinates
                            this->SortAndRemove( tSurfaceCoordMatrix );

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
                            moris::Matrix< DDRMat > tSurfaceCoordMatrix( 2, 3 );
                            tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder.first );
                            tSurfaceCoordMatrix.get_row( 1 ) = t3DcoordMatrix.get_row( tPermutationOrder.second );

                            // sort them counter clock wise in order to be treated the same as other coordinates
                            this->SortAndRemove( tSurfaceCoordMatrix );

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
                        this->create_dbl_sided_cluster( tParamCoordsCell,
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
                            moris::Matrix< DDRMat > tSurfaceCoordMatrix( 2, 3 );
                            tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder.first );
                            tSurfaceCoordMatrix.get_row( 1 ) = t3DcoordMatrix.get_row( tPermutationOrder.second );

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
                            moris::Matrix< DDRMat > tSurfaceCoordMatrix( 2, 3 );
                            tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder.first );
                            tSurfaceCoordMatrix.get_row( 1 ) = t3DcoordMatrix.get_row( tPermutationOrder.second );

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
                        this->create_dbl_sided_cluster( tCutPolygons,
                            tIndicesinCutCell,
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
    Periodic_3D::create_dbl_sided_cluster( Vector< Matrix< DDRMat > > tP,
        Vector< moris_index >                                        &aIndicesinCutCell,
        moris::mtk::Cell const                                            &aInterpCell1,
        moris::mtk::Cell const                                            &aInterpCell2,
        uint                                                               aPairCount,
        moris_index                                                        aPhaseToPhase )
    {
        // requires that we give a unique coordinates and associated index to each of them
        Matrix< DDRMat > tUniqueIntersectedPoints;

        // join all the coordinates together
        for ( size_t i = 0; i < aIndicesinCutCell.size(); i++ )
        {
            tUniqueIntersectedPoints = join_horiz( tUniqueIntersectedPoints, tP( aIndicesinCutCell( i ) ) );
        }

        // remove duplicates of the point
        uint i       = 0;
        uint j       = 1;
        uint tNumCol = tUniqueIntersectedPoints.n_cols();

        while ( j < tNumCol )
        {
            if ( norm( tUniqueIntersectedPoints.get_column( i ) - tUniqueIntersectedPoints.get_column( j ) ) > 0.00000000001 )
            {
                i++;

                tUniqueIntersectedPoints.get_column( i ) = tUniqueIntersectedPoints.get_column( j );

                j++;
            }
            else
            {
                j++;
            }
        }

        tUniqueIntersectedPoints.resize( 2, i + 1 );

        // intersection cell from the algorithm comes as coordinates, convert them to unique coordinate indices
        Vector< Matrix< IndexMat > > tPVertexIndex( aIndicesinCutCell.size(), Matrix< IndexMat >( 1, 3 ) );

        // loop over the cell of coordinates to assign indicies
        for ( size_t i = 0; i < aIndicesinCutCell.size(); i++ )
        {
            if ( tP( aIndicesinCutCell( i ) ).n_cols() != 3 )
            {
                tPVertexIndex( i ) = Matrix< IndexMat >( 1, tP( aIndicesinCutCell( i ) ).n_cols() );
            }

            for ( size_t j = 0; j < tP( aIndicesinCutCell( i ) ).n_cols(); j++ )
            {
                for ( size_t k = 0; k < tUniqueIntersectedPoints.n_cols(); k++ )
                {
                    if ( norm( tP( aIndicesinCutCell( i ) ).get_column( j ) - tUniqueIntersectedPoints.get_column( k ) ) < 0.00000000001 )
                    {
                        tPVertexIndex( i )( j ) = k;
                        break;
                    }
                }
            }
        }

        moris::Matrix< moris::DDRMat > tIPCellCoordinates = aInterpCell1.get_vertex_coords();

        // number of nodes at the surface
        uint tNumSurfaceNodes = tUniqueIntersectedPoints.n_cols();

        // paramteric coordinates leader side 3D
        Matrix< DDRMat > tParametricCoordinates3D( tNumSurfaceNodes, 3 );

        // permutation of the 2D coordinates to 3D
        std::pair< moris_id, moris_id > tCoordOrder = this->permutation_pair( aPairCount );

        // 3rd coordinate that
        Matrix< DDRMat > tDummyCoords( tNumSurfaceNodes, 1, -1 );

        // fill in the coordinate matrix
        tParametricCoordinates3D.get_column( 3 - tCoordOrder.first - tCoordOrder.second ) = tDummyCoords.get_column( 0 );
        tParametricCoordinates3D.get_column( tCoordOrder.first )                          = trans( tUniqueIntersectedPoints.get_row( 0 ) );
        tParametricCoordinates3D.get_column( tCoordOrder.second )                         = trans( tUniqueIntersectedPoints.get_row( 1 ) );

        // initialize the 3D physical coordinates
        Matrix< DDRMat > tPhysicalCoordinates3D( tNumSurfaceNodes, 3 );

        // interpolate the physical coordinates
        moris::mtk::Interpolation::trilinear_interpolation_multivalue( tIPCellCoordinates, tParametricCoordinates3D, tPhysicalCoordinates3D );

        // head vertex of the tetrahedron which lies in the center of the cube
        Matrix< DDRMat > tCoordinatesTop;
        Matrix< DDRMat > tHeadnodeCoords( 1, 3, 0 );
        moris::mtk::Interpolation::trilinear_interpolation( tIPCellCoordinates, tHeadnodeCoords, tCoordinatesTop );

        // geerate the leader vertex data
        this->create_leader_vertices( tPhysicalCoordinates3D, aInterpCell1, aPairCount, tCoordinatesTop );

        // extract the corner coordinates of the IP cell
        tIPCellCoordinates = aInterpCell2.get_vertex_coords();

        // Duplicate the parametric coordinates for the salve side
        Matrix< DDRMat > tFollowerParametricCoordinates3D = tParametricCoordinates3D.copy();

        // Use the same matrices and change the 3rd dimension coordinates
        tDummyCoords.fill( 1.0 );
        tFollowerParametricCoordinates3D.get_column( 3 - tCoordOrder.first - tCoordOrder.second ) = tDummyCoords.get_column( 0 );

        // get physical coordinates of the follower side
        moris::mtk::Interpolation::trilinear_interpolation_multivalue( tIPCellCoordinates, tFollowerParametricCoordinates3D, tPhysicalCoordinates3D );

        // physcal coordinates of the
        moris::mtk::Interpolation::trilinear_interpolation( tIPCellCoordinates, tHeadnodeCoords, tCoordinatesTop );

        // generate the follower vertex data
        this->create_follower_vertices( tPhysicalCoordinates3D, aInterpCell2, aPairCount, tCoordinatesTop );

        tic tTimerCell;

        //////loop to determine number of the cells created
        uint tNumCellsOnOneSide = 0;
        for ( size_t tClusterNum = 0; tClusterNum < aIndicesinCutCell.size(); tClusterNum++ )
        {
            // if polygon break it further to tri by fan triangulation
            if ( tPVertexIndex( tClusterNum ).n_cols() > 3 )
            {
                // find number of newly aded surfaces if it is a polygon
                uint tNumEdges = tPVertexIndex( tClusterNum ).n_cols();

                tNumCellsOnOneSide += ( tNumEdges - 3 );
            }
            else
            {
                tNumCellsOnOneSide++;
            }
        }

        // the outer cell is the side cluster number
        mSideClusterToCells( mNumSideClusters ).resize( tNumCellsOnOneSide );
        mSideClusterToCells( mNumSideClusters + 1 ).resize( tNumCellsOnOneSide );

        std::iota( mSideClusterToCells( mNumSideClusters ).begin(),
            mSideClusterToCells( mNumSideClusters ).end(),
            mNumCells );

        std::iota( mSideClusterToCells( mNumSideClusters + 1 ).begin(),
            mSideClusterToCells( mNumSideClusters + 1 ).end(),
            mNumCells + tNumCellsOnOneSide );

        mCellToVertexIndices( mNumSideClusters ).reserve( tNumCellsOnOneSide * 4 );
        mCellToVertexIndices( mNumSideClusters + 1 ).reserve( tNumCellsOnOneSide * 4 );

        // initialize leader,follower side IG cells
        // Vector< moris::mtk::Cell const * > tLeaderIntegCells;
        // Vector< moris::mtk::Cell const * > tFollowerIntegCells;

        // iterate through each intersection surface and make Ig cells
        for ( size_t tClusterNum = 0; tClusterNum < aIndicesinCutCell.size(); tClusterNum++ )
        {
            // if polygon break it further to tri by fan triangulation
            if ( tPVertexIndex( tClusterNum ).n_cols() > 3 )
            {
                // find number of newly aded surfaces if it is a polygon
                uint tNumEdges = tPVertexIndex( tClusterNum ).n_cols();

                uint tNumTriangles = tNumEdges - 2;

                // loop over the subtriangles created by fan triangulation
                for ( uint tTRI = 0; tTRI < tNumTriangles; tTRI++ )
                {
                    // vertex Index of the cretae sub triangles
                    Matrix< IndexMat > tTmpVertexIndex = { { tPVertexIndex( tClusterNum )( 0 ), tPVertexIndex( tClusterNum )( tTRI + 1 ), tPVertexIndex( tClusterNum )( tTRI + 2 ) } };

                    // leader IG cell
                    this->create_leader_ig_cell( tTmpVertexIndex, aInterpCell1, aPairCount );

                    // follower IG cell
                    this->create_follower_ig_cell( tTmpVertexIndex, aInterpCell2, aPairCount );
                }
            }
            else
            {
                // leader IG cells
                this->create_leader_ig_cell( tPVertexIndex( tClusterNum ), aInterpCell1, aPairCount );

                // follower IG cells
                this->create_follower_ig_cell( tPVertexIndex( tClusterNum ), aInterpCell2, aPairCount );
            }
        }

        mVertexParametricCoords( { mNumParamCoords, mNumParamCoords + tUniqueIntersectedPoints.n_cols() - 1 }, { 0, 2 } )                                         = tParametricCoordinates3D.matrix_data();
        mVertexParametricCoords( { mNumParamCoords + tUniqueIntersectedPoints.n_cols(), mNumParamCoords + 2 * tUniqueIntersectedPoints.n_cols() - 1 }, { 0, 2 } ) = tFollowerParametricCoordinates3D.matrix_data();

        mSideClusterToIPCell( mNumSideClusters )     = aInterpCell1.get_index();
        mSideClusterToIPCell( mNumSideClusters + 1 ) = aInterpCell2.get_index();

        mNumParamCoords += 2 * tUniqueIntersectedPoints.n_cols();

        // save the double sided cluster phase index
        mDoubleSidedClustersIndex( mNumDblSideCluster ) = aPhaseToPhase;

        // increase the count of double sided cluster by 2, 1 for leader and 1 for follower
        mNumSideClusters += 2;

        // increase the count of double sided cluster by 2, 1 for leader and 1 for follower
        mNumCells += 2 * tNumCellsOnOneSide;

        // increase the count of double sided cluster
        mNumDblSideCluster++;
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Periodic_3D::create_leader_ig_cell( Matrix< IndexMat > const &tVertexIndex, moris::mtk::Cell const &aLeaderInterpCell, uint aPairCount )
    {
        // Sort the nodes such that it is consistent with Exodus convention
        Matrix< IndexMat > tVertexIndexCounterClockWise;

        if ( mMeshSideSetPairs( aPairCount )( 0 ) == "5" )
        {
            tVertexIndexCounterClockWise = tVertexIndex;
        }
        else
        {
            tVertexIndexCounterClockWise = { { tVertexIndex( 0 ), tVertexIndex( 2 ), tVertexIndex( 1 ) } };
        }

        // cell of Ig vertices
        Vector< moris_index > tIgCellVertices( 4 );

        // surface vertices
        for ( uint i = 0; i < tVertexIndex.n_cols(); i++ )
        {
            tIgCellVertices( i ) = mSideClusterToVertexIndices( mNumSideClusters )( tVertexIndexCounterClockWise( i ) );
        }

        // head node of tet4
        tIgCellVertices( 3 ) = mSideClusterToVertexIndices( mNumSideClusters )( mSideClusterToVertexIndices( mNumSideClusters ).size() - 1 );

        mCellToVertexIndices( mNumSideClusters ).append( tIgCellVertices );
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Periodic_3D::create_follower_ig_cell( Matrix< IndexMat > const &tVertexIndex, moris::mtk::Cell const &aFollowerInterpCell, uint aPairCount )
    {
        // Sort the nodes such that it is consistent with Exodus convention
        Matrix< IndexMat > tVertexIndexCounterClockWise;

        if ( mMeshSideSetPairs( aPairCount )( 0 ) == "5" )
        {
            tVertexIndexCounterClockWise = { { tVertexIndex( 0 ), tVertexIndex( 2 ), tVertexIndex( 1 ) } };
        }
        else
        {
            tVertexIndexCounterClockWise = tVertexIndex;
        }

        // cell of Ig vertices
        Vector< moris_index > tIgCellVertices( 4 );

        // surface vertices
        for ( uint i = 0; i < tVertexIndex.n_cols(); i++ )
        {
            tIgCellVertices( i ) = mSideClusterToVertexIndices( mNumSideClusters + 1 )( tVertexIndexCounterClockWise( i ) );
        }

        // head node of tet4
        tIgCellVertices( 3 ) = mSideClusterToVertexIndices( mNumSideClusters + 1 )( mSideClusterToVertexIndices( mNumSideClusters ).size() - 1 );

        mCellToVertexIndices( mNumSideClusters + 1 ).append( tIgCellVertices );
    }

    //------------------------------------------------------------------------------------------------------------
    // name the cluster set
    void
    Periodic_3D::construct_add_dbl_sided_set( const moris::Matrix< IndexMat > &tPhaseInteractionTable )
    {
        // resize the coordinate matrix
        mVertexParametricCoords.resize( mNumParamCoords, 3 );
        mVerticesCoords.resize( 3, mNumVertices );

        // resize the double sided cluster data and side cluster data
        mSideClusterToVertexIndices.resize( mNumSideClusters );
        mSideClusterToCells.resize( mNumSideClusters );
        mCellToVertexIndices.resize( mNumSideClusters );
        mSideClusterToIPCell.resize( mNumSideClusters );
        mDoubleSidedClustersIndex.resize( mNumDblSideCluster );
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_3D::offset_vector( moris::Matrix< DDRMat > &tOffsetVector, Vector< std::string > &tFirstSideSetNames, Vector< std::string > &tSecondSideSetNames, uint aPairCount )
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

        // populate the first side clusters
        for ( uint i = 0; i < tFirstSideSetNames.size(); i++ )
        {
            moris::mtk::Set *tSet = tIntegrationMesh->get_set_by_name( tFirstSideSetNames( i ) );

            // get clusters in the second set
            Vector< moris::mtk::Cluster const * > tSetClusters = tSet->get_clusters_on_set();

            tFirstSideClusters.append( tSetClusters );

            if ( tFirstSideClusters.size() > 0 )
            {
                break;
            }
        }

        // populate the second side clusters
        for ( uint i = 0; i < tSecondSideSetNames.size(); i++ )
        {
            moris::mtk::Set *tSet = tIntegrationMesh->get_set_by_name( tSecondSideSetNames( i ) );

            // get clusters in the second set
            Vector< moris::mtk::Cluster const * > tSetClusters = tSet->get_clusters_on_set();

            tSecondSideClusters.append( tSetClusters );

            if ( tSecondSideClusters.size() > 0 )
            {
                break;
            }
        }

        Vector< moris::mtk::Cell const * > const &tSecondCells = tSecondSideClusters( 0 )->get_primary_cells_in_cluster();

        // get the side ordinals of the integration cells in the cluster
        moris::Matrix< moris::IndexMat > tSecondCellOrds = tSecondSideClusters( 0 )->get_cell_side_ordinals();

        // vertices on the second one
        Vector< moris::mtk::Vertex const * > tSecondVertices = tSecondCells( 0 )->get_vertices_on_side_ordinal( tSecondCellOrds( 0 ) );

        moris::Matrix< moris::DDRMat > tSecondVertexCoords = tSecondVertices( 0 )->get_coords();

        moris::Matrix< moris::DDRMat > tSecondSideNormal = tSecondCells( 0 )->compute_outward_side_normal( tSecondCellOrds( 0 ) );

        // get the Integration cell cluster for the side cluster
        Vector< moris::mtk::Cell const * > const &tFirstCells = tFirstSideClusters( 0 )->get_primary_cells_in_cluster();

        // get the side ordinals of the integration cells in the cluster
        moris::Matrix< moris::IndexMat > tFirstCellOrds = tFirstSideClusters( 0 )->get_cell_side_ordinals();

        // get the vertex information
        Vector< moris::mtk::Vertex const * > tFirstVertices = tFirstCells( 0 )->get_vertices_on_side_ordinal( tFirstCellOrds( 0 ) );

        moris::Matrix< moris::DDRMat > tFirstVertexCoords = tFirstVertices( 0 )->get_coords();

        // form the output
        tOffsetVector = dot( tSecondVertexCoords - tFirstVertexCoords, tSecondSideNormal ) * trans( tSecondSideNormal );
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_3D::create_leader_vertices( Matrix< DDRMat > &tPhysicalCoordinates3D, moris::mtk::Cell const &aLeaderInterpCell, uint aPairCount, Matrix< DDRMat > &tNewNodeCoordinates )
    {

        uint tNumSurfaceNodes = tPhysicalCoordinates3D.n_rows();

        // populate the indices
        mSideClusterToVertexIndices( mNumSideClusters ).resize( tNumSurfaceNodes + 1 );

        // generate consecutive vertices for te vertices
        std::iota( mSideClusterToVertexIndices( mNumSideClusters ).begin(),
            mSideClusterToVertexIndices( mNumSideClusters ).end(),
            mNumVertices );

        // loop over surface nodes
        for ( uint i = 0; i < tNumSurfaceNodes; i++ )
        {
            // obtain physical coordinates of the node
            Matrix< DDRMat > tNewCoords = tPhysicalCoordinates3D.get_row( i );

            // save coordinates
            mVerticesCoords( 0, mNumVertices ) = tNewCoords( 0 );
            mVerticesCoords( 1, mNumVertices ) = tNewCoords( 1 );
            mVerticesCoords( 2, mNumVertices ) = tNewCoords( 2 );

            // increment vertices count
            mNumVertices++;
        }

        // save coordinates
        mVerticesCoords( 0, mNumVertices ) = tNewNodeCoordinates( 0 );
        mVerticesCoords( 1, mNumVertices ) = tNewNodeCoordinates( 1 );
        mVerticesCoords( 2, mNumVertices ) = tNewNodeCoordinates( 2 );

        // increment vertices count
        mNumVertices++;
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_3D::create_follower_vertices( Matrix< DDRMat > &tPhysicalCoordinates3D, moris::mtk::Cell const &aFollowerInterpCell, uint aPairCount, Matrix< DDRMat > &tNewNodeCoordinates )
    {

        uint tNumSurfaceNodes = tPhysicalCoordinates3D.n_rows();

        // populate the indices
        mSideClusterToVertexIndices( mNumSideClusters + 1 ).resize( tNumSurfaceNodes + 1 );

        // generate consecutive vertices for te vertices
        std::iota( mSideClusterToVertexIndices( mNumSideClusters + 1 ).begin(),
            mSideClusterToVertexIndices( mNumSideClusters + 1 ).end(),
            mNumVertices );

        for ( uint i = 0; i < tNumSurfaceNodes; i++ )
        {
            // obtain physical coordinates of the node
            Matrix< DDRMat > tNewCoords = tPhysicalCoordinates3D.get_row( i );

            // save coordinates
            mVerticesCoords( 0, mNumVertices ) = tNewCoords( 0 );
            mVerticesCoords( 1, mNumVertices ) = tNewCoords( 1 );
            mVerticesCoords( 2, mNumVertices ) = tNewCoords( 2 );

            // increment vertices count
            mNumVertices++;
        }

        // center node coordinate
        mVerticesCoords( 0, mNumVertices ) = tNewNodeCoordinates( 0 );
        mVerticesCoords( 1, mNumVertices ) = tNewNodeCoordinates( 1 );
        mVerticesCoords( 2, mNumVertices ) = tNewNodeCoordinates( 2 );

        // increment vertices count
        mNumVertices++;
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_3D::rotation_matrix( moris::Matrix< DDRMat > &aRotation, moris::Matrix< DDRMat > &aInverseRotation, uint aPairCount )
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

    std::pair< moris_id, moris_id >
    Periodic_3D::permutation_pair( uint const &aPairCount ) const
    {
        std::pair< moris_id, moris_id > tPermuationOrder;

        if ( mMeshSideSetPairs( aPairCount )( 0 ) == "4" )
        {
            tPermuationOrder = std::make_pair( 2, 1 );
        }
        else if ( mMeshSideSetPairs( aPairCount )( 0 ) == "1" )
        {
            tPermuationOrder = std::make_pair( 0, 2 );
        }
        else
        {
            tPermuationOrder = std::make_pair( 0, 1 );
        }
        return tPermuationOrder;
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_3D::generate_identifier( Vector< mtk::Cluster const * >    &aSideClusters,
        uint                                                                 &aPairCount,
        std::unordered_map< moris::moris_index, Vector< moris_index > > &aBackgroundCellToSideClusterMap ) const
    {
        // reserve memory for output matrix
        uint tOutputSize = aSideClusters.size();

        moris::mtk::Integration_Mesh *tIntegrationMesh = mMeshManager->get_integration_mesh( mMeshIndex );

        // multiplier for the identification number
        uint tMultiplier = tIntegrationMesh->get_num_entities( mtk::EntityRank::NODE );

        uint tFirst;
        uint tSecond;
        // depending on the side pick relevant ijk
        if ( mMeshSideSetPairs( aPairCount )( 0 ) == "4" )
        {
            // pick j and k
            tFirst  = 1;
            tSecond = 2;
        }
        else if ( mMeshSideSetPairs( aPairCount )( 0 ) == "1" )
        {
            // pick i and k
            tFirst  = 0;
            tSecond = 2;
        }
        else
        {
            // pick i and j
            tFirst  = 0;
            tSecond = 1;
        }

        // loop over the clusters
        for ( uint iI = 0; iI < tOutputSize; iI++ )
        {
            // access the ijk postion
            const luint *tIJK = aSideClusters( iI )->get_interpolation_cell().get_ijk();

            // calculate the identifier
            moris_index tIdentifier = tIJK[tFirst] * tMultiplier + tIJK[tSecond];

            // add the cluster index to the relevant bg cell
            aBackgroundCellToSideClusterMap[tIdentifier].push_back( iI );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_3D::elementwise_bruteforce_search(
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

        // Initialize the output
        moris::Matrix< moris::DDUMat > tnc;

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
                this->Intersect( aParamCoordsCell1( iI ), aParamCoordsCell2( iJ ), tP, tnc );

                // If it is a polygon add to the output
                if ( tP.n_cols() > 2 )
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
    Periodic_3D::group_cut_cells( moris::Matrix< IndexMat > const     &aCutTrianglesIdentifier,
        std::unordered_map< moris_index, Vector< moris_index > > &aCutCellIdentifierToCutCellIndex ) const
    {
        // loop over all the cut cells
        for ( uint iI = 0; iI < aCutTrianglesIdentifier.n_rows(); iI++ )
        {
            aCutCellIdentifierToCutCellIndex[aCutTrianglesIdentifier( iI )].push_back( iI );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Periodic_3D::determine_max_size()
    {
        // resize the maps based on the spatial dimension
        mBackgroundCellToSideClusterMap1.resize( 3 );
        mBackgroundCellToSideClusterMap2.resize( 3 );

        // initialize the maximum size of the double side clusters
        uint tNumVertices            = 0;
        uint tNumDoubleSidedClusters = 0;

        // loop over the pairs to detmine the sizes
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

            for ( const auto & tKeyValue : mBackgroundCellToSideClusterMap1( tPairCount ) )
            {
                auto iIJK = tKeyValue.first;
                auto iClusterIndices = tKeyValue.second;

                tNumDoubleSidedClusters += iClusterIndices.size() * mBackgroundCellToSideClusterMap2( tPairCount )[iIJK].size();
            }

            //number of vertices on the left side
            uint tNumVertices1 = 0;
            for ( const auto &iCluster : tSideClusters1 )
            {
                tNumVertices1 += iCluster->get_num_vertices_in_cluster() + 1;
            }

            //number of vertices on the left side
            uint tNumVertices2 = 0;
            for ( const auto &iCluster : tSideClusters2 )
            {
                tNumVertices2 += iCluster->get_num_vertices_in_cluster() + 1;
            }

            //accumulate the vertices size
            tNumVertices += ( 2 * std::max( tNumVertices2, tNumVertices1 ) );
        }

        //resize side cluster and dbl side cluster data
        mSideClusterToVertexIndices.resize( 2 * tNumDoubleSidedClusters );
        mSideClusterToCells.resize( 2 * tNumDoubleSidedClusters );
        mSideClusterToIPCell.resize( 2 * tNumDoubleSidedClusters );
        mCellToVertexIndices.resize( 2 * tNumDoubleSidedClusters );
        mDoubleSidedClustersIndex.resize( tNumDoubleSidedClusters );

        //conservatively multiply number of vertices by 5
        tNumVertices *= 5;

        //set size for coordinate matrix
        mVerticesCoords.set_size( 3, tNumVertices );
        mVertexParametricCoords.set_size( tNumVertices, 3 );
    }

    void
    Periodic_3D::EdgeIntersect(
        moris::Matrix< DDRMat > const &aFirstTRICoords,
        moris::Matrix< DDRMat > const &aSecondTRICoords,
        moris::Matrix< DDRMat >       &aIntersectedPoints,
        moris::Matrix< DDUMat >       &aIntersectVec ) const
    {
        uint tNumIntersections = 0;

        aIntersectVec.set_size( 1, 3, 0 );

        for ( uint i = 0; i < 3; i++ )
        {
            for ( uint j = 0; j < 3; j++ )
            {
                // form the matrix on RHS
                moris::Matrix< DDRMat > b = aSecondTRICoords.get_column( j ) - aFirstTRICoords.get_column( i );

                // matrix on the LHS
                moris::Matrix< DDRMat > A( 2, 2 );
                A.set_column( 0, aFirstTRICoords.get_column( ( i + 1 ) % 3 ) - aFirstTRICoords.get_column( i ) );
                A.set_column( 1, -aSecondTRICoords.get_column( ( j + 1 ) % 3 ) + aSecondTRICoords.get_column( j ) );

                // solve the system
                if ( std::abs( A( 0, 0 ) * A( 1, 1 ) - A( 0, 1 ) * A( 1, 0 ) ) > 0.0000001 )
                {
                    moris::Matrix< DDRMat > R = inv( A ) * b;

                    moris::real eps = 0.01;

                    // Intersection Condition
                    if ( R( 0 ) >= 0 and ( R( 0 ) - 1 ) <= eps and R( 1 ) >= 0 and ( R( 1 ) - 1 ) <= eps )
                    {
                        // grow the matrix to insert the new intersection
                        aIntersectedPoints.resize( 2, tNumIntersections + 1 );

                        aIntersectedPoints.get_column( tNumIntersections ) = aFirstTRICoords.get_column( i ) + R( 0 ) * A.get_column( 0 );

                        aIntersectVec( i ) = 1;

                        // increase number of intersection by 1
                        tNumIntersections++;
                    }
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------

    void
    Periodic_3D::PointsXInY(
        moris::Matrix< moris::DDRMat > const &aFirstTRICoords,
        moris::Matrix< moris::DDRMat > const &aSecondTRICoords,
        moris::Matrix< moris::DDRMat >       &aIntersectedPoints ) const
    {
        uint tNumIntersections = 0;

        // Interior points
        moris::Matrix< moris::DDRMat > v0 = aSecondTRICoords.get_column( 1 ) - aSecondTRICoords.get_column( 0 );
        moris::Matrix< moris::DDRMat > v1 = aSecondTRICoords.get_column( 2 ) - aSecondTRICoords.get_column( 0 );

        // Baricenteric Coordinates
        real d00 = dot( v0, v0 );
        real d01 = dot( v0, v1 );
        real d11 = dot( v1, v1 );

        real id = 1.0 / ( d00 * d11 - d01 * d01 );

        for ( uint i = 0; i < 3; i++ )
        {
            moris::Matrix< moris::DDRMat > v2 = aFirstTRICoords.get_column( i ) - aSecondTRICoords.get_column( 0 );

            real d02 = dot( v0, v2 );
            real d12 = dot( v1, v2 );

            real u = ( d11 * d02 - d01 * d12 ) * id;
            real v = ( d00 * d12 - d01 * d02 ) * id;

            real eps = 0.001;

            if ( u >= -eps and v >= -eps and ( u + v <= 1 + eps ) )
            {
                aIntersectedPoints.resize( 2, tNumIntersections + 1 );

                aIntersectedPoints.get_column( tNumIntersections ) = aFirstTRICoords.get_column( i );

                tNumIntersections++;
            }
        }
    }

    //------------------------------------------------------------------------

    void
    Periodic_3D::SortAndRemove( moris::Matrix< moris::DDRMat > &aIntersectedPoints ) const
    {
        real eps = 0.0001;

        uint tNumColumn = aIntersectedPoints.n_cols();

        if ( tNumColumn > 0 )
        {
            // C
            moris::Matrix< moris::DDRMat > C( 2, 1 );
            C.get_row( 0 ) = sum( aIntersectedPoints.get_row( 0 ) ) / tNumColumn;
            C.get_row( 1 ) = sum( aIntersectedPoints.get_row( 1 ) ) / tNumColumn;

            moris::Matrix< moris::DDRMat > ao( 1, tNumColumn );

            // order polygon corners counter
            for ( uint i = 0; i < tNumColumn; i++ )
            {
                moris::Matrix< moris::DDRMat > d = aIntersectedPoints.get_column( i ) - C;

                ao( i ) = std::atan2( d( 1 ), d( 0 ) );
            }

            moris::Matrix< moris::DDRMat > aoSorted;

            moris::Matrix< moris::DDUMat > IdMatrix( 1, tNumColumn );

            moris::sort( ao, aoSorted, "ascend", 1 );

            // find the index matrix
            for ( uint i = 0; i < tNumColumn; i++ )
            {
                for ( uint j = 0; j < tNumColumn; j++ )
                {
                    if ( ao( j ) == aoSorted( i ) )
                    {
                        IdMatrix( i ) = j;
                        break;
                    }
                }
            }

            // sort the points based on IdMatrix
            moris::Matrix< moris::DDRMat > tSortedIntersectedPoints( 2, tNumColumn );

            for ( uint i = 0; i < tNumColumn; i++ )
            {
                tSortedIntersectedPoints.get_column( i ) = aIntersectedPoints.get_column( IdMatrix( i ) );
            }

            // remove duplicates of the point
            uint i = 0;
            uint j = 1;

            while ( j < tNumColumn )
            {
                if ( norm( tSortedIntersectedPoints.get_column( i ) - tSortedIntersectedPoints.get_column( j ) ) > eps )
                {
                    i++;

                    tSortedIntersectedPoints.get_column( i ) = tSortedIntersectedPoints.get_column( j );

                    j++;
                }
                else
                {
                    j++;
                }
            }

            tSortedIntersectedPoints.resize( 2, i + 1 );

            aIntersectedPoints = tSortedIntersectedPoints;
        }
    }

    //------------------------------------------------------------------------

    void
    Periodic_3D::Intersect(
        moris::Matrix< moris::DDRMat > const &aFirstTRICoords,
        moris::Matrix< moris::DDRMat > const &aSecondTRICoords,
        moris::Matrix< moris::DDRMat >       &aIntersectedPoints,
        moris::Matrix< moris::DDUMat >       &aIntersectVec ) const
    {
        this->EdgeIntersect( aFirstTRICoords, aSecondTRICoords, aIntersectedPoints, aIntersectVec );

        Matrix< DDRMat > P1;

        this->PointsXInY( aFirstTRICoords, aSecondTRICoords, P1 );

        if ( P1.n_cols() > 1 )
        {
            aIntersectVec.set_size( 1, 3, 1 );
        }

        aIntersectedPoints = join_horiz( aIntersectedPoints, P1 );

        this->PointsXInY( aSecondTRICoords, aFirstTRICoords, P1 );

        aIntersectedPoints = join_horiz( aIntersectedPoints, P1 );

        this->SortAndRemove( aIntersectedPoints );
    }

}// namespace moris::mig

