/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Periodic_Boundary_Condition_Helper.cpp
 *
 */

#include "cl_MTK_Periodic_Boundary_Condition_Helper.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Set.hpp"

#include "fn_norm.hpp"
#include "fn_sum.hpp"
#include "cl_MTK_Mesh.hpp"
namespace moris
{
    namespace mtk
    {

        Periodic_Boundary_Condition_Helper::Periodic_Boundary_Condition_Helper(
                std::shared_ptr< Mesh_Manager > aMeshManager,
                moris_index                     aMeshIndex,
                moris::ParameterList&           aParameterList )
                : mMeshManager( aMeshManager )
                , mMeshIndex( aMeshIndex )
        {
            // get the periodic mesh set names
            std::string tMeshSideSetNames = aParameterList.get< std::string >( "periodic_side_set_pair" );

            // store the names in mMeshSideSetPairs
            string_to_cell_of_cell( aParameterList.get< std::string >( "periodic_side_set_pair" ), mMeshSideSetPairs );
        }

        void
        Periodic_Boundary_Condition_Helper::setup_periodic_boundary_conditions()
        {
            // Access the integration mesh
            moris::mtk::Integration_Mesh* tIntegrationMesh = mMeshManager->get_integration_mesh( mMeshIndex );

            // output the mesh
            //            moris::mtk::Writer_Exodus tWriter(tIntegrationMesh);
            //            tWriter.write_mesh("", "1DBeam.exo", "", "temp.exo");
            //            tWriter.close_file();

            // Double sided side set clusters
            moris::Cell< moris::Cell< Cluster const * > > tDoubleSideSetClusters( mMeshSideSetPairs.size() );

            // Double sided side set
            moris::Cell< mtk::Double_Side_Set* > tDblSideSet( mMeshSideSetPairs.size() );

            // iterate through periodic pairs to
            for ( uint tPairCounter = 0; tPairCounter < mMeshSideSetPairs.size(); tPairCounter++ )
            {
                // get the first side set from the mesh
                moris::mtk::Set* tSet1 = tIntegrationMesh->get_set_by_name( mMeshSideSetPairs( tPairCounter )( 0 ) );

                // Get the spatial dim of the mesh
                uint tSpatialDim = tIntegrationMesh->get_spatial_dim();
                // get clusters in the the first set
                moris::Cell< moris::mtk::Cluster const * > tSetClusters = tSet1->get_clusters_on_set();

                // Cell to store 1st side IDs of each integration cell - each member of the cell holds one cell
                moris::Cell< moris::Matrix< moris::IdMat > > tFirstIds( tSetClusters.size() );

                // get the second side set from the mesh
                moris::mtk::Set* tSet2 = tIntegrationMesh->get_set_by_name( mMeshSideSetPairs( tPairCounter )( 1 ) );

                // get clusters in the second set
                moris::Cell< moris::mtk::Cluster const * > tSetClusters2 = tSet2->get_clusters_on_set();

                // Cell to store 2nd side IDs of each integration cell - each member of the cell holds one cell
                moris::Cell< moris::Matrix< moris::IdMat > > tSecondIds( tSetClusters.size() );

                // Ids of the cell clusters corresponding to each other
                moris::Matrix< moris::DDUMat > tPairedIndices;
                tPairedIndices.set_size( tSetClusters.size(), 1, tSetClusters.size() + 1 );

                uint tNumVertexinCluster = tSetClusters( 0 )->get_num_vertices_in_cluster();
                // Permutation of the nodes matrix
                moris::Matrix< moris::IdMat > tPermutation;
                tPermutation.set_size( tNumVertexinCluster, 1, 0 );

                // flag to find the permutation
                bool tFoundPermutation = false;
                switch ( tSpatialDim )
                {
                    case 3:
                    {
                        // X,Y and Z coordinates of the the first cluster to store and average
                        moris::Matrix< moris::DDRMat > tXVec;
                        moris::Matrix< moris::DDRMat > tYVec;
                        moris::Matrix< moris::DDRMat > tZVec;

                        // set size
                        tXVec.set_size( tNumVertexinCluster, 1, 0 );
                        tYVec.set_size( tNumVertexinCluster, 1, 0 );
                        tZVec.set_size( tNumVertexinCluster, 1, 0 );

                        // X,Y and Z coordinates of the the first cluster to store and average
                        moris::Matrix< moris::DDRMat > tXVec2;
                        moris::Matrix< moris::DDRMat > tYVec2;
                        moris::Matrix< moris::DDRMat > tZVec2;

                        // set size
                        tXVec2.set_size( tNumVertexinCluster, 1, 0 );
                        tYVec2.set_size( tNumVertexinCluster, 1, 0 );
                        tZVec2.set_size( tNumVertexinCluster, 1, 0 );

                        // checks to see if there it fits the criteria that surfaces have to be parallel to x-y-z and each other
                        for ( uint tClusterIndex = 0; tClusterIndex < tSetClusters.size(); tClusterIndex++ )
                        {
                            // First Cluster
                            moris::Matrix< moris::DDRMat > tVertexCoordsInCluster = tSetClusters( tClusterIndex )->get_vertex_coords_in_cluster();
                            tXVec += tVertexCoordsInCluster.get_column( 0 );
                            tYVec += tVertexCoordsInCluster.get_column( 1 );
                            tZVec += tVertexCoordsInCluster.get_column( 2 );

                            // second cluster
                            moris::Matrix< moris::DDRMat > tVertexCoordsInCluster2 = tSetClusters2( tClusterIndex )->get_vertex_coords_in_cluster();
                            tXVec2 += tVertexCoordsInCluster2.get_column( 0 );
                            tYVec2 += tVertexCoordsInCluster2.get_column( 1 );
                            tZVec2 += tVertexCoordsInCluster2.get_column( 2 );
                        }

                        // X, Y and Z average of 1st cluster set(plane)
                        moris::real tXCoordAvg = sum( tXVec ) / ( tNumVertexinCluster * tSetClusters.size() );
                        moris::real tYCoordAvg = sum( tYVec ) / ( tNumVertexinCluster * tSetClusters.size() );
                        moris::real tZCoordAvg = sum( tZVec ) / ( tNumVertexinCluster * tSetClusters.size() );

                        // X, Y and Z average of 2nd cluster set(plane)
                        moris::real tXCoordAvg2 = sum( tXVec2 ) / ( tNumVertexinCluster * tSetClusters.size() );
                        moris::real tYCoordAvg2 = sum( tYVec2 ) / ( tNumVertexinCluster * tSetClusters.size() );
                        moris::real tZCoordAvg2 = sum( tZVec2 ) / ( tNumVertexinCluster * tSetClusters.size() );

                        // value of the offset and its direction to check
                        moris::real tOffset;
                        std::string tOffsetDir;
                        moris::uint tOffDirEnum;

                        // iterate through possible cases for a periodic cube
                        if ( ( std::abs( tXCoordAvg - tXCoordAvg2 ) < 0.001 ) && ( std::abs( tYCoordAvg - tYCoordAvg2 ) < 0.001 ) )
                        {
                            tOffset     = std::abs( tZCoordAvg - tZCoordAvg2 );
                            tOffsetDir  = "z";
                            tOffDirEnum = 2;
                        }
                        else if ( ( std::abs( tYCoordAvg - tYCoordAvg2 ) < 0.001 ) && ( std::abs( tZCoordAvg - tZCoordAvg2 ) < 0.001 ) )
                        {
                            tOffset     = std::abs( tXCoordAvg - tXCoordAvg2 );
                            tOffsetDir  = "x";
                            tOffDirEnum = 0;
                        }
                        else if ( ( std::abs( tXCoordAvg - tXCoordAvg2 ) < 0.001 ) && ( std::abs( tZCoordAvg - tZCoordAvg2 ) < 0.001 ) )
                        {
                            tOffset     = std::abs( tYCoordAvg - tYCoordAvg2 );
                            tOffsetDir  = "y";
                            tOffDirEnum = 1;
                        }
                        else
                        {
                            MORIS_ERROR( false, "This version of the boundary condition is not implemented, the surfaces should be in x-y-z direction and parallel to each other" );
                        }

                        // Iterate through 1st set of clusters
                        for ( uint tFirstClusterIndex = 0; tFirstClusterIndex < tSetClusters.size(); tFirstClusterIndex++ )
                        {
                            moris::Matrix< moris::IdMat > tVertexIdsInCluster = tSetClusters( tFirstClusterIndex )->get_vertex_ids_in_cluster();

                            // found a pair flag to save some iteration cost
                            for ( uint tSecondClusterIndex = 0; tSecondClusterIndex < tSetClusters2.size(); tSecondClusterIndex++ )
                            {
                                bool tFoundPairFlag = false;
                                for ( uint tFoundIndex = 0; tFoundIndex < tSetClusters2.size(); tFoundIndex++ )
                                {
                                    uint tIndex     = tPairedIndices( tFoundIndex, 0 );
                                    uint tIndexDiff = tIndex - tSecondClusterIndex;
                                    if ( tIndexDiff == 0 )
                                    {
                                        tFoundPairFlag = true;
                                        break;
                                    }
                                }
                                if ( tFoundPairFlag )
                                {
                                    continue;
                                }

                                // vertex Coords of the 1st cluster set member
                                moris::Matrix< moris::DDRMat > tVertexCoordsInCluster = tSetClusters( tFirstClusterIndex )->get_vertex_coords_in_cluster();

                                // vertex IDs of the 2nd cluster set member
                                moris::Matrix< moris::IdMat > tVertexIdsInCluster2 = tSetClusters2( tSecondClusterIndex )->get_vertex_ids_in_cluster();

                                // vertex Coords of the 2nd cluster set member
                                moris::Matrix< moris::DDRMat > tVertexCoordsInCluster2 = tSetClusters2( tSecondClusterIndex )->get_vertex_coords_in_cluster();

                                // possible cases of offset in different directions, subtract the offset to see two clusters have the same Coords
                                switch ( tOffDirEnum )
                                {
                                    case 0:
                                    {
                                        if ( tVertexCoordsInCluster2( 0, 0 ) > tVertexCoordsInCluster( 0, 0 ) )
                                        {
                                            Matrix< DDRMat > tOffsetVec;
                                            tOffsetVec.set_size( tNumVertexinCluster, 1, tOffset );
                                            tVertexCoordsInCluster2.get_column( 0 ) = tVertexCoordsInCluster2.get_column( 0 ) - tOffsetVec;
                                        }
                                        else
                                        {
                                            Matrix< DDRMat > tOffsetVec;
                                            tOffsetVec.set_size( tNumVertexinCluster, 1, tOffset );
                                            tVertexCoordsInCluster.get_column( 0 ) = tVertexCoordsInCluster.get_column( 0 ) - tOffsetVec;
                                        }

                                        break;
                                    }

                                    case 1:
                                    {
                                        if ( tVertexCoordsInCluster2( 0, 1 ) > tVertexCoordsInCluster( 0, 1 ) )
                                        {
                                            Matrix< DDRMat > tOffsetVec;
                                            tOffsetVec.set_size( tNumVertexinCluster, 1, tOffset );
                                            tVertexCoordsInCluster2.get_column( 1 ) = tVertexCoordsInCluster2.get_column( 1 ) - tOffsetVec;
                                        }
                                        else
                                        {
                                            Matrix< DDRMat > tOffsetVec;
                                            tOffsetVec.set_size( tNumVertexinCluster, 1, tOffset );
                                            tVertexCoordsInCluster.get_column( 1 ) = tVertexCoordsInCluster.get_column( 1 ) - tOffsetVec;
                                        }

                                        break;
                                    }

                                    case 2:
                                    {
                                        if ( tVertexCoordsInCluster2( 0, 2 ) > tVertexCoordsInCluster( 0, 2 ) )
                                        {
                                            Matrix< DDRMat > tOffsetVec;
                                            tOffsetVec.set_size( tNumVertexinCluster, 1, tOffset );
                                            tVertexCoordsInCluster2.get_column( 2 ) = tVertexCoordsInCluster2.get_column( 2 ) - tOffsetVec;
                                        }
                                        else
                                        {
                                            Matrix< DDRMat > tOffsetVec;
                                            tOffsetVec.set_size( tNumVertexinCluster, 1, tOffset );
                                            tVertexCoordsInCluster.get_column( 2 ) = tVertexCoordsInCluster.get_column( 2 ) - tOffsetVec;
                                        }

                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false, "Not Implemented yet" );
                                    }
                                }
                                // the matrix of Coords for the cluster should have similar norms
                                if ( std::abs( norm( tVertexCoordsInCluster2 ) - norm( tVertexCoordsInCluster ) ) < 0.0000000001 )
                                {
                                    // go through te nodes once to find the right permutation
                                    if ( tFoundPermutation == false )
                                    {
                                        for ( uint tFirstNode = 0; tFirstNode < tNumVertexinCluster; tFirstNode++ )
                                        {
                                            for ( uint tSecondNode = 0; tSecondNode < tNumVertexinCluster; tSecondNode++ )
                                            {
                                                Matrix< moris::DDRMat > tDiff = tVertexCoordsInCluster.get_row( tFirstNode ) - tVertexCoordsInCluster2.get_row( tSecondNode );
                                                if ( norm( tDiff ) < 0.001 )
                                                {
                                                    tPermutation.get_row( tFirstNode ) = tSecondNode;

                                                    break;
                                                }
                                            }

                                            tFoundPermutation = true;
                                        }
                                    }

                                    // add the Id's of the second cluster to be matched
                                    tSecondIds( tFirstClusterIndex ) = tVertexIdsInCluster2;

                                    // ID of the 2nd cluster matching the fisrt one
                                    tPairedIndices.get_row( tFirstClusterIndex ) = tSecondClusterIndex;

                                    break;
                                }
                            }

                            // add the IDs of the first cluster to be paired later
                            tFirstIds( tFirstClusterIndex ) = tVertexIdsInCluster;
                        }

                        // resize the cell of the set of cluster- each cluster set has tSetClusters2.size() clusters in it
                        tDoubleSideSetClusters( tPairCounter ).resize( tSetClusters2.size() );

                        // iterate through the cluster to pair the cluster and nodes
                        for ( uint tClusterPairCounter = 0; tClusterPairCounter < tSetClusters.size(); tClusterPairCounter++ )
                        {
                            // ID's of the second cluster
                            moris::Matrix< moris::IdMat > tSecondPairIds = tSecondIds( tClusterPairCounter );

                            // vertex pairing for individual clusters
                            moris::Cell< moris::mtk::Vertex const * > tVertexPairing( tNumVertexinCluster );
                            tVertexPairing( 0 ) = &tIntegrationMesh->get_mtk_vertex( tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id( tSecondPairIds( 0, tPermutation( 0, 0 ) ), EntityRank::NODE ) );
                            tVertexPairing( 1 ) = &tIntegrationMesh->get_mtk_vertex( tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id( tSecondPairIds( 0, tPermutation( 1, 0 ) ), EntityRank::NODE ) );
                            tVertexPairing( 2 ) = &tIntegrationMesh->get_mtk_vertex( tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id( tSecondPairIds( 0, tPermutation( 2, 0 ) ), EntityRank::NODE ) );
                            tVertexPairing( 3 ) = &tIntegrationMesh->get_mtk_vertex( tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id( tSecondPairIds( 0, tPermutation( 3, 0 ) ), EntityRank::NODE ) );

                            // construct the double side cluster
                            Double_Side_Cluster* tDoubleSideCluster = new Double_Side_Cluster( tSetClusters( tClusterPairCounter ), tSetClusters2( tPairedIndices( tClusterPairCounter, 0 ) ), tVertexPairing );

                            // add the cluster to the mesh, the only purpose is to be able to delete later in the deconstructor
                            tIntegrationMesh->add_double_sided_cluster( tDoubleSideCluster );

                            // save the constructed double sided cluster
                            tDoubleSideSetClusters( tPairCounter )( tClusterPairCounter ) = tDoubleSideCluster;
                        }

                        // name the cluster set
                        std::string               tDoubleSideSetName = "Periodic" + std::to_string( tPairCounter );
                        moris::Matrix< IndexMat > tColors            = { { 0 } };

                        // Construct the double side set
                        tDblSideSet( tPairCounter ) = new moris::mtk::Double_Side_Set( tDoubleSideSetName, tDoubleSideSetClusters( tPairCounter ), tColors, tSpatialDim );

                        // add double sided periodic boundary condition to the integration mesh
                        tIntegrationMesh->add_double_side_set( tDblSideSet( tPairCounter ) );

                        break;
                    }

                    case 2:
                    {
                        // X,Y and Z coordinates of the the first cluster to store and average
                        moris::Matrix< moris::DDRMat > tXVec;
                        moris::Matrix< moris::DDRMat > tYVec;

                        // set size
                        tXVec.set_size( tNumVertexinCluster, 1, 0 );
                        tYVec.set_size( tNumVertexinCluster, 1, 0 );

                        // X,Y and Z coordinates of the the first cluster to store and average
                        moris::Matrix< moris::DDRMat > tXVec2;
                        moris::Matrix< moris::DDRMat > tYVec2;

                        // set size
                        tXVec2.set_size( tNumVertexinCluster, 1, 0 );
                        tYVec2.set_size( tNumVertexinCluster, 1, 0 );

                        // checks to see if there it fits the criteria that surfaces have to be parallel to x-y-z and each other
                        for ( uint tClusterIndex = 0; tClusterIndex < tSetClusters.size(); tClusterIndex++ )
                        {
                            // First Cluster
                            moris::Matrix< moris::DDRMat > tVertexCoordsInCluster = tSetClusters( tClusterIndex )->get_vertex_coords_in_cluster();
                            tXVec += tVertexCoordsInCluster.get_column( 0 );
                            tYVec += tVertexCoordsInCluster.get_column( 1 );

                            // second cluster
                            moris::Matrix< moris::DDRMat > tVertexCoordsInCluster2 = tSetClusters2( tClusterIndex )->get_vertex_coords_in_cluster();
                            tXVec2 += tVertexCoordsInCluster2.get_column( 0 );
                            tYVec2 += tVertexCoordsInCluster2.get_column( 1 );
                        }

                        // X, Y and Z average of 1st cluster set(plane)
                        moris::real tXCoordAvg = sum( tXVec ) / ( tNumVertexinCluster * tSetClusters.size() );
                        moris::real tYCoordAvg = sum( tYVec ) / ( tNumVertexinCluster * tSetClusters.size() );

                        // X, Y and Z average of 2nd cluster set(plane)
                        moris::real tXCoordAvg2 = sum( tXVec2 ) / ( tNumVertexinCluster * tSetClusters.size() );
                        moris::real tYCoordAvg2 = sum( tYVec2 ) / ( tNumVertexinCluster * tSetClusters.size() );

                        // value of the offset and its direction to check
                        moris::real tOffset;
                        std::string tOffsetDir;
                        moris::uint tOffDirEnum;

                        // iterate through possible cases for a periodic cube
                        if ( ( std::abs( tYCoordAvg - tYCoordAvg2 ) < 0.001 ) )
                        {
                            tOffset     = std::abs( tXCoordAvg - tXCoordAvg2 );
                            tOffsetDir  = "x";
                            tOffDirEnum = 0;
                        }
                        else if ( ( std::abs( tXCoordAvg - tXCoordAvg2 ) < 0.001 ) )
                        {
                            tOffset     = std::abs( tYCoordAvg - tYCoordAvg2 );
                            tOffsetDir  = "y";
                            tOffDirEnum = 1;
                        }
                        else
                        {
                            MORIS_ERROR( false, "This version of the boundary condition is not implemented, the surfaces should be in x-y-z direction and parallel to each other" );
                        }

                        // Iterate through 1st set of clusters
                        for ( uint tFirstClusterIndex = 0; tFirstClusterIndex < tSetClusters.size(); tFirstClusterIndex++ )
                        {
                            moris::Matrix< moris::IdMat > tVertexIdsInCluster = tSetClusters( tFirstClusterIndex )->get_vertex_ids_in_cluster();

                            // found a pair flag to save some iteration cost
                            for ( uint tSecondClusterIndex = 0; tSecondClusterIndex < tSetClusters2.size(); tSecondClusterIndex++ )
                            {
                                bool tFoundPairFlag = false;
                                for ( uint tFoundIndex = 0; tFoundIndex < tSetClusters2.size(); tFoundIndex++ )
                                {
                                    uint tIndex     = tPairedIndices( tFoundIndex, 0 );
                                    uint tIndexDiff = tIndex - tSecondClusterIndex;
                                    if ( tIndexDiff == 0 )
                                    {
                                        tFoundPairFlag = true;
                                        break;
                                    }
                                }
                                if ( tFoundPairFlag )
                                {
                                    continue;
                                }

                                // vertex Coords of the 1st cluster set member
                                moris::Matrix< moris::DDRMat > tVertexCoordsInCluster = tSetClusters( tFirstClusterIndex )->get_vertex_coords_in_cluster();

                                // vertex IDs of the 2nd cluster set member
                                moris::Matrix< moris::IdMat > tVertexIdsInCluster2 = tSetClusters2( tSecondClusterIndex )->get_vertex_ids_in_cluster();

                                // vertex Coords of the 2nd cluster set member
                                moris::Matrix< moris::DDRMat > tVertexCoordsInCluster2 = tSetClusters2( tSecondClusterIndex )->get_vertex_coords_in_cluster();

                                // possible cases of offset in different directions, subtract the offset to see two clusters have the same Coords
                                switch ( tOffDirEnum )
                                {
                                    case 0:
                                    {
                                        if ( tVertexCoordsInCluster2( 0, 0 ) > tVertexCoordsInCluster( 0, 0 ) )
                                        {
                                            Matrix< DDRMat > tOffsetVec;
                                            tOffsetVec.set_size( tNumVertexinCluster, 1, tOffset );
                                            tVertexCoordsInCluster2.get_column( 0 ) = tVertexCoordsInCluster2.get_column( 0 ) - tOffsetVec;
                                        }
                                        else
                                        {
                                            Matrix< DDRMat > tOffsetVec;
                                            tOffsetVec.set_size( tNumVertexinCluster, 1, tOffset );
                                            tVertexCoordsInCluster.get_column( 0 ) = tVertexCoordsInCluster.get_column( 0 ) - tOffsetVec;
                                        }

                                        break;
                                    }

                                    case 1:
                                    {
                                        if ( tVertexCoordsInCluster2( 0, 1 ) > tVertexCoordsInCluster( 0, 1 ) )
                                        {
                                            Matrix< DDRMat > tOffsetVec;
                                            tOffsetVec.set_size( tNumVertexinCluster, 1, tOffset );
                                            tVertexCoordsInCluster2.get_column( 1 ) = tVertexCoordsInCluster2.get_column( 1 ) - tOffsetVec;
                                        }
                                        else
                                        {
                                            Matrix< DDRMat > tOffsetVec;
                                            tOffsetVec.set_size( tNumVertexinCluster, 1, tOffset );
                                            tVertexCoordsInCluster.get_column( 1 ) = tVertexCoordsInCluster.get_column( 1 ) - tOffsetVec;
                                        }

                                        break;
                                    }

                                    default:
                                    {
                                        MORIS_ERROR( false, "Not Implemented yet" );
                                    }
                                }
                                // the matrix of Coords for the cluster should have similar norms
                                if ( std::abs( norm( tVertexCoordsInCluster2 ) - norm( tVertexCoordsInCluster ) ) < 0.0000000001 )
                                {
                                    // go through te nodes once to find the right permutation
                                    if ( tFoundPermutation == false )
                                    {
                                        for ( uint tFirstNode = 0; tFirstNode < tNumVertexinCluster; tFirstNode++ )
                                        {
                                            for ( uint tSecondNode = 0; tSecondNode < tNumVertexinCluster; tSecondNode++ )
                                            {
                                                Matrix< moris::DDRMat > tDiff = tVertexCoordsInCluster.get_row( tFirstNode ) - tVertexCoordsInCluster2.get_row( tSecondNode );
                                                if ( norm( tDiff ) < 0.001 )
                                                {
                                                    tPermutation.get_row( tFirstNode ) = tSecondNode;

                                                    break;
                                                }
                                            }

                                            tFoundPermutation = true;
                                        }
                                    }

                                    // add the Id's of the second cluster to be matched
                                    tSecondIds( tFirstClusterIndex ) = tVertexIdsInCluster2;

                                    // ID of the 2nd cluster matching the fisrt one
                                    tPairedIndices.get_row( tFirstClusterIndex ) = tSecondClusterIndex;

                                    break;
                                }
                            }

                            // add the IDs of the first cluster to be paired later
                            tFirstIds( tFirstClusterIndex ) = tVertexIdsInCluster;
                        }

                        // resize the cell of the set of cluster- each cluster set has tSetClusters2.size() clusters in it
                        tDoubleSideSetClusters( tPairCounter ).resize( tSetClusters2.size() );

                        // iterate through the cluster to pair the cluster and nodes
                        for ( uint tClusterPairCounter = 0; tClusterPairCounter < tSetClusters.size(); tClusterPairCounter++ )
                        {
                            // ID's of the second cluster
                            moris::Matrix< moris::IdMat > tSecondPairIds = tSecondIds( tClusterPairCounter );

                            // vertex pairing for individual clusters
                            moris::Cell< moris::mtk::Vertex const * > tVertexPairing( tNumVertexinCluster );
                            tVertexPairing( 0 ) = &tIntegrationMesh->get_mtk_vertex( tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id( tSecondPairIds( 0, tPermutation( 0, 0 ) ), EntityRank::NODE ) );
                            tVertexPairing( 1 ) = &tIntegrationMesh->get_mtk_vertex( tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id( tSecondPairIds( 0, tPermutation( 1, 0 ) ), EntityRank::NODE ) );

                            // construct the double side cluster
                            Double_Side_Cluster* tDoubleSideCluster = new Double_Side_Cluster( tSetClusters( tClusterPairCounter ), tSetClusters2( tPairedIndices( tClusterPairCounter, 0 ) ), tVertexPairing );

                            // add the cluster to the mesh, the only purpose is to be able to delete later in the deconstructor
                            tIntegrationMesh->add_double_sided_cluster( tDoubleSideCluster );

                            // save the constructed double sided cluster
                            tDoubleSideSetClusters( tPairCounter )( tClusterPairCounter ) = tDoubleSideCluster;
                        }

                        // name the cluster set
                        std::string               tDoubleSideSetName = "Periodic" + std::to_string( tPairCounter );
                        moris::Matrix< IndexMat > tColors            = { { 0 } };

                        // Construct the double side set
                        tDblSideSet( tPairCounter ) = new moris::mtk::Double_Side_Set( tDoubleSideSetName, tDoubleSideSetClusters( tPairCounter ), tColors, tSpatialDim );

                        // add double sided periodic boundary condition to the integration mesh
                        tIntegrationMesh->add_double_side_set( tDblSideSet( tPairCounter ) );

                        break;
                    }
                }
            }

            // FIXME: communicate double sided side sets constructed here
        }
    }    // namespace mtk
}    // namespace moris
