/*
 * cl_MTK_Intersection_Detect_2D.cpp
 *
 *  Created on: Jun 7, 2021
 *      Author: momo
 */
#include "cl_MTK_Intersection_Detect_2D.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Integration_Mesh.hpp"


#include "cl_MTK_Side_Cluster.hpp"
#include "typedefs.hpp"
#include "cl_Mesh_Enums.hpp"
#include "fn_unique.hpp"
#include<unordered_map>

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
                std::shared_ptr<moris::mtk::Mesh_Manager> aMeshManager,
                moris::moris_index                        aMeshIndex,
                moris::ParameterList &                    aParameterList,
                moris::uint                               aNumBulkPhases):mMeshManager(aMeshManager),
                        mMeshIndex(aMeshIndex),
                        mNumBulkPhases(aNumBulkPhases)
                        {
            //get the periodic mesh set names
            std::string tMeshSideSetNames = aParameterList.get<std::string>( "periodic_side_set_pair" );

            //store the names in mMeshSideSetPairs
            string_to_cell_of_cell(aParameterList.get<std::string>( "periodic_side_set_pair" ), mMeshSideSetPairs);

            //initialize the mesh intersection data to store the id and index of the added entities
            mIntersectedMeshData.set_up_external_entity_data(mMeshManager->get_integration_mesh( mMeshIndex ));
                        }

        //---------------------------------------------------------------------

        Intersection_Detect_2D::~Intersection_Detect_2D()
        {
            //delete all the pointer data
            for(auto p : mDoubleSidedClusters)
            {
                delete p;
            }
            for(auto p : mMasterSidedClusters)
            {
                delete p;
            }
            for(auto p : mSlaveSidedClusters)
            {
                delete p;
            }
            for(auto p : mMasterSideCells)
            {
                delete p;
            }
            for(auto p : mSlaveSideCells)
            {
                delete p;
            }
            for(auto p : mMasterVertices)
            {
                delete p;
            }
            for(auto p : mSlaveVertices)
            {
                delete p;
            }
        }

        //---------------------------------------------------------------------

        void
        Intersection_Detect_2D::perform()
        {
            // initialize phase interaction table
            Matrix<IndexMat> tPhaseInteractionTable( mNumBulkPhases, mNumBulkPhases );

            for( uint Ii = 0; Ii < mNumBulkPhases ; Ii++ )
            {
                for(uint Ij = 0 ; Ij < mNumBulkPhases ; Ij++)
                {
                    tPhaseInteractionTable( Ii, Ij )= Ij + Ii *  mNumBulkPhases;
                }
            }

            //loop over the pairs
            for( uint tPairCount = 0 ; tPairCount < mMeshSideSetPairs.size() ; tPairCount++ )
            {
                //Keep track of time
                Tracer tTracer("MTK", "No Type", "Searching For Pairs");

                //find the rotation and inverse rotation
                //rotation goes from 1d parametric to 2d parametric
                //inverse other way around
                moris::Matrix<DDRMat >     tRotation;
                moris::Matrix<DDRMat >     tInverseRotation;
                this->rotation_matrix( tRotation, tInverseRotation, tPairCount );

                //obtain the integration mesh
                moris::mtk::Integration_Mesh* tIntegrationMesh = mMeshManager->get_integration_mesh( mMeshIndex );

                //prerequisite offset vector and sets on each periodic surface
                moris::Matrix<DDRMat >     tOffsetVector;
                moris::Cell< std::string > tFirstSideSetNames;
                moris::Cell< std::string > tSecondSideSetNames;

                //find the offset of two surfaces
                this->offset_vector( tOffsetVector, tFirstSideSetNames, tSecondSideSetNames, tPairCount );

                //find the offset matrix
                moris::Matrix < DDRMat > tVecTransform( 2, 1, 1);
                moris::Matrix< moris::DDRMat >  tOffsetMatrix = tVecTransform * tOffsetVector ;

                //initialize all the clusters on each side
                moris::Cell<mtk::Cluster const *> tSideClusters1;
                moris::Cell<mtk::Cluster const *> tSideClusters2;

                //cells on the side set are locally indexed and they are assigned a phase to each one
                moris::Cell< moris_index > tLocalCellNumToColor1 ;

                for (uint i = 0; i < tFirstSideSetNames.size() ; i++)
                {
                    //get the set on the first side
                    moris::mtk::Set* tSet = tIntegrationMesh->get_set_by_name( tFirstSideSetNames(i) );

                    //get color(phase) associated with the set
                    Matrix<IndexMat> tColor = tSet->get_set_colors();

                    // get clusters in the second set
                    moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet->get_clusters_on_set();

                    //Assign phases to cluster
                    tLocalCellNumToColor1.resize( tLocalCellNumToColor1.size()+ tSetClusters.size(), tColor( 0 ) );

                    //append the clusters the list of the first side clusters
                    tSideClusters1.append( tSetClusters );
                }

                //cells on the side set are locally indexed and they are assigned a phase to each one
                moris::Cell< moris_index > tLocalCellNumToColor2 ;

                for (uint i = 0; i < tSecondSideSetNames.size() ; i++)
                {
                    //get the set on the first side
                    moris::mtk::Set* tSet = tIntegrationMesh->get_set_by_name( tSecondSideSetNames( i ) );

                    // get clusters in the second set
                    moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet->get_clusters_on_set();

                    //get color(phase) associated with the set
                    Matrix<IndexMat> tColor = tSet->get_set_colors();

                    //Assign phases to clusters
                    tLocalCellNumToColor2.resize( tLocalCellNumToColor2.size()+ tSetClusters.size(), tColor( 0 ) );

                    //append the clusters the list of the first side clusters
                    tSideClusters2.append(tSetClusters);
                }

                // Ids of the cell clusters corresponding to each other
                moris::Matrix< moris::DDBMat > tPairedClusterIndicesFlag;
                tPairedClusterIndicesFlag.set_size( tSideClusters2.size(), mNumBulkPhases, 0 );

                //loop over the clusters on the first side
                for(uint tClusterNum1 = 0; tClusterNum1 < tSideClusters1.size() ;  tClusterNum1++ )
                {
                    //loop over the clusters on the second side
                    for(uint tClusterNum2 = 0; tClusterNum2 < tSideClusters2.size() ;  tClusterNum2++)
                    {
                        //check if the second side cluster has been visited and paired
                        if( tPairedClusterIndicesFlag( tClusterNum2, tLocalCellNumToColor1 ( tClusterNum1 ) ) == 1 )
                        {
                            continue;
                        }
                        //interpolation cell 1st side
                        moris::mtk::Cell const & tInterpCell1 = tSideClusters1( tClusterNum1 )->get_interpolation_cell();

                        //interpolation cell 2nd side
                        moris::mtk::Cell const & tInterpCell2 = tSideClusters2( tClusterNum2 )->get_interpolation_cell();

                        //coordinates to compare
                        Matrix<DDRMat> tInterpCoords1 = tInterpCell1.get_cell_geometric_coords_on_side_ordinal( std::stoi( mMeshSideSetPairs( tPairCount )( 0 ) ) - 1 );

                        Matrix<DDRMat> tInterpCoords2 = tInterpCell2.get_cell_geometric_coords_on_side_ordinal( std::stoi( mMeshSideSetPairs( tPairCount )( 1 ) ) - 1 );

                        //subtract the offset to compare the coords
                        tInterpCoords2 = tInterpCoords2 - tOffsetMatrix;

                        // matching interpolation cell criteria
                        bool tMatched = false;

                        //criteria for matching cells
                        if( std::abs( norm( tInterpCoords2 ) - norm( tInterpCoords1 ) ) < 0.0000000001 )
                        {
                            uint s = 0 ;

                            //check they actually have the same points
                            for ( uint i = 0 ; i < 2 ; i++ )
                            {
                                for( uint j = 0 ; j < 2 ; j++)
                                {
                                    if ( norm(tInterpCoords2.get_row(i)  - tInterpCoords1.get_row(j) ) < 0.00000000001 )
                                    {
                                        s++;
                                        break;
                                    }
                                }
                            }

                            //if they match
                            if( s == 2 )
                            {
                                tMatched = true;

                            }

                            //if the interpolation cells match
                            if( tMatched )
                            {

                                //mark that second side cluster is paired so
                                tPairedClusterIndicesFlag( tClusterNum2, tLocalCellNumToColor1 ( tClusterNum1 )) = 1;

                                //the case where left and right clusters are trivial
                                if ( tSideClusters1( tClusterNum1 )->is_trivial() and tSideClusters2( tClusterNum2 )->is_trivial() )
                                {
                                    // Ids of the cell clusters corresponding to each other
                                    moris::Matrix< moris::DDBMat > tCompMatrix;
                                    tCompMatrix.set_size( 1, mNumBulkPhases, 1 );

                                    tPairedClusterIndicesFlag.get_row( tClusterNum2 ) = tCompMatrix.get_row( 0 );

                                    //obtain the phases of each side
                                    moris_index tPhase2  = tLocalCellNumToColor2 ( tClusterNum2 );
                                    moris_index tPhase1  = tLocalCellNumToColor1 ( tClusterNum1 );

                                    //phase interaction table value
                                    moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                                    //parametric 1d coordinates of the intersection points
                                    moris::Cell< Matrix<DDRMat> >tParamCoordsCell2 (1);
                                    tParamCoordsCell2(0) = { {-1.0, 1.0 } };

                                    //creates the pair of the single sided clusters
                                    this->make_new_pairs( tParamCoordsCell2, tInterpCell1, tInterpCell2, tPairCount, tPhaseToPhaseIndex);
                                }

                                //case where master side is trivial
                                else if ( tSideClusters1( tClusterNum1 )->is_trivial() == true and tSideClusters2( tClusterNum2 )->is_trivial() == false )
                                {

                                    // Ids of the cell clusters corresponding to each other
                                    moris::Matrix< moris::DDBMat > tCompMatrix;
                                    tCompMatrix.set_size( 1, mNumBulkPhases, 1 );

                                    tPairedClusterIndicesFlag.get_row( tClusterNum2 ) = tCompMatrix.get_row( 0 );

                                    //get the Integration cell cluster for the side cluster
                                    moris::Cell<moris::mtk::Cell const *> const & tCells2 = tSideClusters2( tClusterNum2 )->get_primary_cells_in_cluster();

                                    // get the side ordinals of the integration cells in the cluster
                                    moris::Matrix<moris::IndexMat> tCellOrds2 = tSideClusters2( tClusterNum2 )->get_cell_side_ordinals();

                                    //matrix to store the local coords
                                    moris::Cell < moris::Matrix <DDRMat> > tParamCoordsCell2;

                                    for( uint i = 0 ; i < tCells2.size() ; i++ )
                                    {
                                        //find vertices on the side
                                        moris::Cell<moris::mtk::Vertex const *> tVertices = tCells2( i )->get_vertices_on_side_ordinal( tCellOrds2( i ) );

                                        //parametric coordinates of each line
                                        moris::Matrix <DDRMat> tParamCoordsMatrix2( 1, 2 );

                                        for ( uint j = 0 ; j < tVertices.size() ; j++ )
                                        {
                                            //temporary matrix of 2d parametric coordinates
                                            Matrix<DDRMat> tTmpCoordMatrix = tSideClusters2( tClusterNum2 )->get_vertex_local_coordinate_wrt_interp_cell( tVertices( j ) ).get_row( 0 );

                                            //convert to desired form
                                            tTmpCoordMatrix = tTmpCoordMatrix * tInverseRotation ;

                                            //extract the ones needed
                                            tParamCoordsMatrix2( 0, j ) = tTmpCoordMatrix( 1 );
                                        }

                                        //sort the matrix if the order is not ascending
                                        if ( tParamCoordsMatrix2 (1) < tParamCoordsMatrix2 (0) )
                                        {
                                            tParamCoordsMatrix2 = { {tParamCoordsMatrix2 (1),tParamCoordsMatrix2 (0) } };
                                        }

                                        //store each triangle in the the cell
                                        tParamCoordsCell2.push_back( tParamCoordsMatrix2 );
                                    }

                                    //obtain the phases of each side
                                    moris_index tPhase2  = tLocalCellNumToColor2 ( tClusterNum2 );
                                    moris_index tPhase1  = tLocalCellNumToColor1 ( tClusterNum1 );

                                    //phase interaction table value
                                    moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                                    //creates the pair of the single sided clusters
                                    this->make_new_pairs( tParamCoordsCell2, tInterpCell1, tInterpCell2, tPairCount, tPhaseToPhaseIndex );
                                }

                                //case where slave side is trivial
                                else if ( tSideClusters1( tClusterNum1 )->is_trivial() == false and tSideClusters2( tClusterNum2 )->is_trivial() == true )
                                {
                                    //get the Integration cell cluster for the side cluster
                                    moris::Cell<moris::mtk::Cell const *> const & tCells1 = tSideClusters1( tClusterNum1 )->get_primary_cells_in_cluster();

                                    // get the side ordinals of the integration cells in the cluster
                                    moris::Matrix<moris::IndexMat> tCellOrds1 = tSideClusters1( tClusterNum1 )->get_cell_side_ordinals();

                                    //parametric coordinates of each line
                                    moris::Cell< moris::Matrix <DDRMat> > tParamCoordsCell1;

                                    for( uint i = 0 ; i < tCells1.size() ; i++ )
                                    {
                                        //find vertices on the side
                                        moris::Cell<moris::mtk::Vertex const *> tVertices = tCells1(i)->get_vertices_on_side_ordinal( tCellOrds1(i) );

                                        moris::Matrix <DDRMat> tParamCoordsMatrix1( 1, 2 );

                                        for ( uint j = 0 ; j < tVertices.size() ; j++ )
                                        {
                                            //temporary matrix of 2d parametric coordinates
                                            Matrix<DDRMat> tTmpCoordMatrix = tSideClusters1( tClusterNum1 )->get_vertex_local_coordinate_wrt_interp_cell( tVertices( j ) ).get_row( 0 );

                                            //convert to desired form
                                            tTmpCoordMatrix = tTmpCoordMatrix * tInverseRotation ;

                                            //extract the ones needed
                                            tParamCoordsMatrix1( 0, j ) = tTmpCoordMatrix(1);
                                        }

                                        //sort the matrix if the order is not ascending
                                        if ( tParamCoordsMatrix1 (1) < tParamCoordsMatrix1 (0) )
                                        {
                                            tParamCoordsMatrix1 = { {tParamCoordsMatrix1 (1),tParamCoordsMatrix1 (0) } };
                                        }

                                        //store each triangle in the the cell
                                        tParamCoordsCell1.push_back(tParamCoordsMatrix1);
                                    }

                                    //obtain the phases of each side
                                    moris_index tPhase2  = tLocalCellNumToColor2 ( tClusterNum2 );
                                    moris_index tPhase1  = tLocalCellNumToColor1 ( tClusterNum1 );

                                    //phase interaction table value
                                    moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                                    //creates the pair of the single sided clusters
                                    this->make_new_pairs( tParamCoordsCell1, tInterpCell1, tInterpCell2 ,tPairCount,tPhaseToPhaseIndex);

                                }

                                else
                                {
                                    //get the Integration cell cluster for the side cluster
                                    moris::Cell<moris::mtk::Cell const *> const & tCells1 = tSideClusters1( tClusterNum1 )->get_primary_cells_in_cluster();

                                    // get the side ordinals of the integration cells in the cluster
                                    moris::Matrix<moris::IndexMat> tCellOrds1 = tSideClusters1( tClusterNum1 )->get_cell_side_ordinals();

                                    //matrix to store the local coords
                                    moris::Cell< moris::Matrix <DDRMat> > tParamCoordsCell1;

                                    for( uint i = 0 ; i < tCells1.size() ; i++ )
                                    {
                                        //vertices on the side
                                        moris::Cell<moris::mtk::Vertex const *> tVertices = tCells1(i)->get_geometric_vertices_on_side_ordinal( tCellOrds1(i) );

                                        //parametric coords matrix of each line
                                        moris::Matrix <DDRMat> tParamCoordsMatrix1( 1, 2 );

                                        for ( uint j = 0 ; j < tVertices.size() ; j++ )
                                        {
                                            //temporary matrix of 2d parametric coordinates
                                            Matrix<DDRMat> tTmpCoordMatrix = tSideClusters1( tClusterNum1 )->get_vertex_local_coordinate_wrt_interp_cell( tVertices( j ) ).get_row( 0 );

                                            //convert to desired form
                                            tTmpCoordMatrix = tTmpCoordMatrix * tInverseRotation ;

                                            //extract the ones needed
                                            tParamCoordsMatrix1( 0, j ) = tTmpCoordMatrix(1);
                                        }

                                        //store each triangle in the the cell
                                        tParamCoordsCell1.push_back(tParamCoordsMatrix1);
                                    }

                                    //get the Integration cell cluster for the side cluster
                                    moris::Cell<moris::mtk::Cell const *> const & tCells2 = tSideClusters2( tClusterNum2 )->get_primary_cells_in_cluster();

                                    // get the side ordinals of the integration cells in the cluster
                                    moris::Matrix<moris::IndexMat> tCellOrds2 = tSideClusters2( tClusterNum2 )->get_cell_side_ordinals();

                                    //matrix to store the local coords
                                    moris::Cell < moris::Matrix <DDRMat> > tParamCoordsCell2;

                                    for( uint i = 0 ; i < tCells2.size() ; i++ )
                                    {
                                        //vertices on the side
                                        moris::Cell<moris::mtk::Vertex const *> tVertices = tCells2(i)->get_geometric_vertices_on_side_ordinal( tCellOrds2(i) );

                                        //parametric coords matrix of each line
                                        moris::Matrix <DDRMat> tParamCoordsMatrix2(1,2);

                                        for ( uint j = 0 ; j < tVertices.size() ; j++ )
                                        {

                                            Matrix<DDRMat> tTmpCoordMatrix = tSideClusters2( tClusterNum2 )->get_vertex_local_coordinate_wrt_interp_cell( tVertices( j ) ).get_row( 0 );

                                            //convert to desired form
                                            tTmpCoordMatrix = tTmpCoordMatrix * tInverseRotation ;

                                            //extract the ones needed
                                            tParamCoordsMatrix2( 0, j ) = tTmpCoordMatrix(1);
                                        }

                                        // //store each triangle in the the cell
                                        tParamCoordsCell2.push_back( tParamCoordsMatrix2 );
                                    }

                                    //sort the parametric coordinates
                                    moris::Matrix< moris::DDRMat > tParamCoords2Sorted;
                                    moris::Matrix< moris::DDRMat > tParamCoords1Sorted;

                                    // sort the parametric coords of each side
                                    moris::sort( tParamCoordsCell2(0), tParamCoords2Sorted ,"ascend", 1 );
                                    moris::sort( tParamCoordsCell1(0), tParamCoords1Sorted ,"ascend", 1 );

                                    //define the overlap interval
                                    moris::Matrix< moris::DDRMat > tOverLap(1,2);

                                    //check to see if 2 line segments overlap
                                    if ( ( (tParamCoords1Sorted(1) - tParamCoords2Sorted(0) ) >= 0 ) and ( (tParamCoords2Sorted(1) - tParamCoords1Sorted(0) ) >= 0 ) )
                                    {
                                        //overlap 2 line segments
                                        tOverLap = { { std::max(tParamCoords1Sorted(0) , tParamCoords2Sorted(0) ), std::min(tParamCoords1Sorted(1) , tParamCoords2Sorted(1))} };

                                        if( std::abs( tOverLap(0) - tOverLap(1) ) < MORIS_REAL_EPS )
                                        {
                                            continue;
                                        }
                                        //obtain the phases of each side
                                        moris_index tPhase2  = tLocalCellNumToColor2 ( tClusterNum2 );
                                        moris_index tPhase1  = tLocalCellNumToColor1 ( tClusterNum1 );

                                        //phase interaction table value
                                        moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                                        //make the output in the desired form
                                        moris::Cell< moris::Matrix< moris::DDRMat > > tOverLapCell(1);
                                        tOverLapCell(0) = tOverLap;

                                        //creates the pair of the single sided clusters
                                        this->make_new_pairs(  tOverLapCell, tInterpCell1, tInterpCell2 ,tPairCount, tPhaseToPhaseIndex);
                                    }
                                }
                            } /*end loop for matching IP cells */
                        }
                    } /*end loop for clusters on the 2nd side */
                }/*end loop for clusters on the 1st side */
            }/*end loop foe the pairs */

            //construct the double sided set and add it to the mesh
            this->constrcuct_add_dbl_sided_set(tPhaseInteractionTable );

        }

        //------------------------------------------------------------------------

        void
        Intersection_Detect_2D::make_new_pairs( moris::Cell< Matrix<DDRMat> > tP,
                moris::mtk::Cell const & aInterpCell1,
                moris::mtk::Cell const & aInterpCell2,
                uint aPairCount,
                moris_index aPhaseToPhase)
        {

            //convert to desired form
            moris::Matrix<DDRMat >     tRotation;
            moris::Matrix<DDRMat >     tInverseRotation;
            this->rotation_matrix( tRotation, tInverseRotation, aPairCount );

            // make a matrix of intersection points
            Matrix< DDRMat > tUniqueIntersectedPoints;

            //populate the matrix
            for (size_t i = 0 ; i < tP.size(); i++ )
            {
                tUniqueIntersectedPoints = join_horiz( tUniqueIntersectedPoints, tP( i ) );
            }

            //create master and slave side vertices
            moris::Cell<moris::mtk::Vertex* > tMasterVerticesCell = this->create_master_vertices(tUniqueIntersectedPoints,aInterpCell1, aPairCount);
            moris::Cell<moris::mtk::Vertex* > tSlaveVerticesCell = this->create_slave_vertices(tUniqueIntersectedPoints,aInterpCell2, aPairCount);

            //initialize master,slave side cluster data
            moris::Cell<moris::mtk::Cell const *> tMasterIntegCells;
            moris::Cell<moris::mtk::Cell const *> tSlaveIntegCells;

            //loop over the side clusters needs to be created( It is 1 in 2D)
            for( size_t tClusterNum = 0 ; tClusterNum < tP.size() ; tClusterNum++ )
            {
                //create master IG cell
                moris::mtk::Cell const * tMasterIgCell = this->create_master_ig_cell( tMasterVerticesCell, aInterpCell1, aPairCount) ;
                tMasterIntegCells.push_back(tMasterIgCell);

                //create slave IG cell
                moris::mtk::Cell const * tSlaveIgCell = this->create_slave_ig_cell( tSlaveVerticesCell, aInterpCell2, aPairCount ) ;
                tSlaveIntegCells.push_back(tSlaveIgCell);
            }

            //create a matrix for ordinals for slave and master cells
            moris::Matrix<moris::IndexMat>          tMasterIntegrationCellSideOrdinals(1,tP.size(), 0);
            moris::Matrix<moris::IndexMat>          tSlaveIntegrationCellSideOrdinals(1,tP.size(), 0);

            //parametric coords of master and slave vertices
            moris::Matrix< DDRMat > tMasterParamCoords(tUniqueIntersectedPoints.n_cols(), 2 );
            moris::Matrix< DDRMat > tSlaveParamCoords(tUniqueIntersectedPoints.n_cols(), 2 ) ;

            //loop through
            for(uint  i = 0 ; i < tUniqueIntersectedPoints.n_cols() ; i++ )
            {
                moris::Matrix< DDRMat > tMasterParamCoord = { { -1,  tUniqueIntersectedPoints( 0, i )} };
                moris::Matrix< DDRMat > tSlaveParamCoord = { { +1,  tUniqueIntersectedPoints( 0, i ) } };
                tMasterParamCoord = tMasterParamCoord * tRotation;
                tSlaveParamCoord  = tSlaveParamCoord * tRotation;
                tMasterParamCoords.get_row( i ) = tMasterParamCoord.get_row( 0 );
                tSlaveParamCoords.get_row( i )  = tSlaveParamCoord.get_row( 0 );
            }



            //create constant version of the vertices on cluster
            moris::Cell<moris::mtk::Vertex const *> tMasterVerticesConst;
            moris::Cell<moris::mtk::Vertex const *> tSlaveVerticesConst;

            for(uint i = 0 ; i < tMasterVerticesCell.size() - 1 ; i++)
            {
                tMasterVerticesConst.push_back(tMasterVerticesCell(i));
                tSlaveVerticesConst.push_back(tSlaveVerticesCell(i));
            }

            moris::mtk::Side_Cluster_ISC* tMasterSideCluster = new moris::mtk::Side_Cluster_ISC(false,
                    & aInterpCell1,
                    tMasterIntegCells,
                    tMasterIntegrationCellSideOrdinals,
                    tMasterVerticesConst,
                    tMasterParamCoords);

            moris::mtk::Side_Cluster_ISC* tSlaveSideCluster  = new moris::mtk::Side_Cluster_ISC(false,
                    & aInterpCell2,
                    tSlaveIntegCells,
                    tSlaveIntegrationCellSideOrdinals,
                    tSlaveVerticesConst,
                    tSlaveParamCoords);

            // create double side cluster
            mtk::Double_Side_Cluster* tDblSideCluster = new mtk::Double_Side_Cluster(
                    tMasterSideCluster,
                    tSlaveSideCluster,
                    tSlaveSideCluster->mVerticesInCluster);

            //store double side cluster
            mDoubleSidedClusters.push_back(tDblSideCluster);
            mDoubleSidedClustersIndex.push_back(aPhaseToPhase);

            //store side clusters
            mMasterSidedClusters.push_back( tMasterSideCluster );
            mSlaveSidedClusters.push_back( tSlaveSideCluster );

            //store the vertices to prevent memory leak
            mMasterVertices.append(tMasterVerticesConst);
            mSlaveVertices.append( tSlaveVerticesConst);

            // Append the integration cells
            mMasterSideCells.append(tMasterIntegCells);
            mSlaveSideCells.append(tSlaveIntegCells);
        }


        //------------------------------------------------------------------------------------------------------------

        moris::mtk::Cell const *
        Intersection_Detect_2D::create_master_ig_cell( moris::Cell<moris::mtk::Vertex *> aMasterVertices,  moris::mtk::Cell const & aMasterInterpCell, uint aPairCount)
        {
            Matrix<IndexMat> tPermutation;

            if( mMeshSideSetPairs(aPairCount)(0) == "4")
            {
                tPermutation = { { 1, 0, 2 } };
            }
            else
            {
                tPermutation = { { 0, 1, 2 } };
            }

            //arrange the permutation based on the set
            moris::Cell<moris::mtk::Vertex *> tMasterVertices( aMasterVertices.size() );
            for(uint  i = 0 ; i < tPermutation.n_cols(); i++)
            {
                tMasterVertices( i ) = aMasterVertices( tPermutation( i ) );
            }

            //allocate id and index from intersection data
            moris_index tCellIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::ELEMENT) ;
            moris_id tCellId= mIntersectedMeshData.allocate_entity_ids_external_entity_data(1, EntityRank::ELEMENT);

            // cell info for TRI 3
            mtk::Cell_Info_Factory tCellInfoFactory;
            std::shared_ptr<moris::mtk::Cell_Info> tLinearCellInfo = tCellInfoFactory.create_cell_info_sp(Geometry_Type::TRI,mtk::Interpolation_Order::LINEAR);

            //create the new integration cell
            moris::mtk::Cell* tIgCell = new Cell_ISC(
                    tCellId,
                    tCellIndex,
                    aMasterInterpCell.get_owner(),
                    tLinearCellInfo,
                    tMasterVertices);

            mIntersectedMeshData.update_first_available_index_external_data( tCellIndex+1, EntityRank::ELEMENT);

            return tIgCell;
        }

        //------------------------------------------------------------------------------------------------------------

        moris::mtk::Cell const *
        Intersection_Detect_2D::create_slave_ig_cell( moris::Cell<moris::mtk::Vertex *> aSlaveVertices ,  moris::mtk::Cell const & aSlaveInterpCell, uint aPairCount)
        {

            Matrix<IndexMat> tPermutation;

            if( mMeshSideSetPairs(aPairCount)(0) == "4")
            {
                tPermutation = { { 0, 1, 2 } };
            }
            else
            {
                tPermutation = { { 1, 0, 2 } };
            }

            //arrange the permutation based on the set
            moris::Cell<moris::mtk::Vertex *> tSlaveVertices( aSlaveVertices.size() );
            for(uint  i = 0 ; i < tPermutation.n_cols(); i++)
            {
                tSlaveVertices( i ) = aSlaveVertices( tPermutation( i ) );
            }

            //allocate id and index from the data
            moris_index tCellIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::ELEMENT) ;
            moris_id tCellId = mIntersectedMeshData.allocate_entity_ids_external_entity_data(1, EntityRank::ELEMENT);

            // cell info  and geometry for TRI 3
            mtk::Cell_Info_Factory tCellInfoFactory;
            std::shared_ptr<moris::mtk::Cell_Info> tLinearCellInfo = tCellInfoFactory.create_cell_info_sp(Geometry_Type::TRI,mtk::Interpolation_Order::LINEAR);

            //create the IG cell
            moris::mtk::Cell* tIgCell = new Cell_ISC(
                    tCellId,
                    tCellIndex,
                    aSlaveInterpCell.get_owner(),
                    tLinearCellInfo,
                    tSlaveVertices);

            mIntersectedMeshData.update_first_available_index_external_data( tCellIndex+1, EntityRank::ELEMENT);

            return tIgCell;
        }

        //------------------------------------------------------------------------------------------------------------
        //name the cluster set
        void Intersection_Detect_2D::constrcuct_add_dbl_sided_set ( moris::Matrix < IndexMat > tPhaseInteractionTable )
        {
            Tracer tTracer("MTK", "Double Sided Set", "Construction of Double Sided Sets");

            //initialize all the double sided sets
            moris::Cell<moris::Cell<mtk::Cluster const*>> tDoubleSideSets( tPhaseInteractionTable.numel() );

            //get the integration mesh
            moris::mtk::Integration_Mesh*  tIntegMesh  = mMeshManager->get_integration_mesh( mMeshIndex );

            uint tSpatialDim = tIntegMesh->get_spatial_dim();

            //assign double sided clusters to the right double sided set
            for( uint i = 0 ; i < mDoubleSidedClusters.size(); i++ )
            {
                tDoubleSideSets(mDoubleSidedClustersIndex(i)).push_back(mDoubleSidedClusters(i));
            }

            //add double sided set to the mesh and name them
            for( size_t aBulkPhase1 = 0 ; aBulkPhase1 <  tPhaseInteractionTable.n_cols() ; aBulkPhase1++)
            {
                for( size_t aBulkPhase2 = 0 ; aBulkPhase2 <  tPhaseInteractionTable.n_cols() ; aBulkPhase2++ )
                {
                    uint tPhaseToPhaseIndex = tPhaseInteractionTable( aBulkPhase1, aBulkPhase2 );

                    //check if double sided set is empty
                    if( tDoubleSideSets(tPhaseToPhaseIndex).size() > 0)
                    {
                        //name the double sided set
                        std::string                  tDoubleSideSetName     = "P"+std::to_string(aBulkPhase1)+std::to_string(aBulkPhase2);

                        //colors of the double sided set,  needs to be modified
                        moris::Matrix<moris::IndexMat>      tColors         = {{0}};

                        // Construct the double side set
                        moris::mtk::Set* tDblSideSet = new moris::mtk::Double_Side_Set(tDoubleSideSetName,  tDoubleSideSets(tPhaseToPhaseIndex), tColors, tSpatialDim);

                        // add double sided periodic boundary condition to the integration mesh
                        tIntegMesh->add_double_side_set( tDblSideSet );
                    }
                }
            }
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect_2D::offset_vector(moris::Matrix<DDRMat > &tOffsetVector,  moris::Cell< std::string > & tFirstSideSetNames,moris::Cell< std::string > &tSecondSideSetNames,uint aPairCount )
        {
            // Integration mesh
            moris::mtk::Integration_Mesh*  tIntegrationMesh  = mMeshManager->get_integration_mesh( mMeshIndex );

            //get all the sidesets
            moris::Cell<std::string> tSetNames = tIntegrationMesh->get_set_names(EntityRank::EDGE);

            //loop through all the set names and finds sets attached to input side numbers
            std::string tPairFirstName = "SideSet_";
            tPairFirstName.append(mMeshSideSetPairs(aPairCount)(0));

            //loop through all the set names and finds sets attached to input side numbers
            std::string tPairSecondName = "SideSet_";
            tPairSecondName.append(mMeshSideSetPairs(aPairCount)(1));

            //loop over the names to get all of the names
            for (uint i = 0; i < tSetNames.size() ; i++)
            {
                if( tSetNames(i).find(tPairFirstName) != std::string::npos )
                {
                    tFirstSideSetNames.push_back(tSetNames(i));
                }
                else if ( tSetNames(i).find(tPairSecondName) != std::string::npos )
                {
                    tSecondSideSetNames.push_back(tSetNames(i));
                }
                else
                {
                    continue;
                }

            }

            //initialize all the clusters on each side
            moris::Cell<mtk::Cluster const *> tFirstSideClusters;
            moris::Cell<mtk::Cluster const *> tSecondSideClusters;

            for(uint i = 0 ; i < tFirstSideSetNames.size() ; i++)
            {
                moris::mtk::Set* tSet = tIntegrationMesh->get_set_by_name( tFirstSideSetNames( i ) );

                // get clusters in the second set
                moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet->get_clusters_on_set();

                tFirstSideClusters.append(tSetClusters);
            }

            //populate the second side clusters
            for(uint i = 0 ; i < tSecondSideSetNames.size(); i++ )
            {
                moris::mtk::Set* tSet = tIntegrationMesh->get_set_by_name( tSecondSideSetNames( i ) );

                // get clusters in the second set
                moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet->get_clusters_on_set();

                tSecondSideClusters.append(tSetClusters);
            }


            //Second side info
            //get the Integration cell cluster for the side cluster
            moris::Cell<moris::mtk::Cell const *> const & tSecondCells = tSecondSideClusters( 0 )->get_primary_cells_in_cluster();

            // get the side ordinals of the integration cells in the cluster
            moris::Matrix<moris::IndexMat> tSecondCellOrds = tSecondSideClusters( 0 )->get_cell_side_ordinals();

            //each cluster has only 2 vertices
            moris::Cell<moris::mtk::Vertex const *> tSecondVertices = tSecondCells( 0 )->get_vertices_on_side_ordinal( tSecondCellOrds( 0 ) );
            moris::Matrix< moris::DDRMat > tSecondVertexCoords1 = tSecondVertices(0)->get_coords();
            moris::Matrix< moris::DDRMat > tSecondVertexCoords2 = tSecondVertices(1)->get_coords();


            //calculate the normal of the right segment,  need something better than std::abs to determine the outward
            moris::real tMagnitude = norm( tSecondVertexCoords2 - tSecondVertexCoords1  );
            moris::Matrix<moris::DDRMat> tSecondNormal = {{std::abs( tSecondVertexCoords1( 1 ) - tSecondVertexCoords2( 1 ) ) , std::abs(tSecondVertexCoords2( 0 ) - tSecondVertexCoords1( 0 ) ) } };
            tSecondNormal =  tSecondNormal / tMagnitude ;

            //First side info
            //get the Integration cell cluster for the side cluster
            moris::Cell<moris::mtk::Cell const *> const & tFirstCells = tFirstSideClusters( 0 )->get_primary_cells_in_cluster();

            // get the side ordinals of the integration cells in the cluster
            moris::Matrix<moris::IndexMat> tFirstCellOrds = tFirstSideClusters( 0 )->get_cell_side_ordinals();

            //get the vertex information
            moris::Cell<moris::mtk::Vertex const *> tFirstVertices = tFirstCells( 0 )->get_vertices_on_side_ordinal( tFirstCellOrds( 0 ) );

            moris::Matrix< moris::DDRMat > tFirstVertexCoords = tFirstVertices(0)->get_coords();

            // distance value between two surfaces
            moris::real tOffsetVectorVal =  std::abs( dot( tFirstVertexCoords-tSecondVertexCoords2 , tSecondNormal) );

            // construct the offset vector, needs to be modified
            tOffsetVector = tOffsetVectorVal*tSecondNormal;
        }

        // ----------------------------------------------------------------------------

        moris::Cell< moris::mtk::Vertex *>
        Intersection_Detect_2D::create_master_vertices(Matrix < DDRMat> tUniqueIntersectedPoints,moris::mtk::Cell const & aMasterInterpCell, uint aPairCount )
        {
            moris::Matrix<DDRMat >     tRotation;
            moris::Matrix<DDRMat >     tInverseRotation;
            this->rotation_matrix( tRotation, tInverseRotation, aPairCount );

            uint tNumSurfaceNodes = tUniqueIntersectedPoints.n_cols();

            //get the interpolation mesh
            moris::mtk::Interpolation_Mesh* tInterpMesh = mMeshManager->get_interpolation_mesh( mMeshIndex );

            //initialize the vertices of the cell ( +1 for the head of the tet4 )
            moris::Cell <moris::mtk::Vertex * > tVertcies( tNumSurfaceNodes+1 ) ;

            // allocate ids for all the vertices
            moris_id tVertexId = mIntersectedMeshData.allocate_entity_ids_external_entity_data( tNumSurfaceNodes+1,  EntityRank::NODE);

            //get the indices of the interpolation cell nodes
            moris::Matrix<moris::IndexMat> tElementNodes =
                    tInterpMesh->get_entity_connected_to_entity_loc_inds(
                            aMasterInterpCell.get_index(),
                            moris::EntityRank::ELEMENT,
                            moris::EntityRank::NODE);

            // coordinates of nodes attached to element to interpolate for physical coordinates
            moris::Matrix<moris::DDRMat> tCoordinates( tElementNodes.numel() , 2);

            for (uint  i = 0 ; i < tElementNodes.numel() ; i++ )
            {
                moris::Matrix<moris::DDRMat> tmp =  tInterpMesh->get_node_coordinate( tElementNodes( i ) );
                tCoordinates.get_row(i) = tmp.get_row(0);
            }

            //loop over surface nodes
            for( uint i = 0 ; i < tNumSurfaceNodes ; i++)
            {
                //form the parametric coordinates of the nodes
                moris::Matrix<DDRMat> tParamCoordsRelativeToElem = { { -1.0, tUniqueIntersectedPoints( 0, i ) } };

                tParamCoordsRelativeToElem = tParamCoordsRelativeToElem * tRotation ;

                //obtain physical coordinates of the node
                Matrix <DDRMat > tNewCoords;
                moris::mtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem, tNewCoords);

                moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;

                //construct the new vertex
                tVertcies ( i ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex,  tNewCoords ) ;

                //update the id and index
                mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );
                tVertexId++;
            }

            // center node coordinates
            moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
            Matrix <DDRMat > tParamCoordsRelativeToElem = {{0.0, 0.0}};
            moris::mtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem, tNewNodeCoordinates);


            moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;
            mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );

            //add the top vertex of TRI#
            tVertcies( tNumSurfaceNodes ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewNodeCoordinates) ;

            return tVertcies;
        }

        // ----------------------------------------------------------------------------

        moris::Cell< moris::mtk::Vertex *>
        Intersection_Detect_2D::create_slave_vertices(Matrix < DDRMat> tUniqueIntersectedPoints,moris::mtk::Cell const & aSlaveInterpCell, uint aPairCount )
        {

            moris::Matrix<DDRMat >     tRotation;
            moris::Matrix<DDRMat >     tInverseRotation;
            this->rotation_matrix( tRotation, tInverseRotation, aPairCount );


            uint tNumSurfaceNodes = tUniqueIntersectedPoints.n_cols();

            //get the interpolation mesh
            moris::mtk::Interpolation_Mesh* tInterpMesh = mMeshManager->get_interpolation_mesh( mMeshIndex );

            //initialize the vertices of the cell ( +1 for the head of the tet4 )
            moris::Cell <moris::mtk::Vertex * > tVertcies( tNumSurfaceNodes+1 ) ;

            // allocate ids for all the vertices
            moris_id tVertexId = mIntersectedMeshData.allocate_entity_ids_external_entity_data( tNumSurfaceNodes+1,  EntityRank::NODE);

            //get the indices of the interpolation cell nodes
            moris::Matrix<moris::IndexMat> tElementNodes =
                    tInterpMesh->get_entity_connected_to_entity_loc_inds(
                            aSlaveInterpCell.get_index(),
                            moris::EntityRank::ELEMENT,
                            moris::EntityRank::NODE);

            // coordinates of nodes attached to element to interpolate for physical coordinates
            moris::Matrix<moris::DDRMat> tCoordinates( tElementNodes.numel() , 2);

            for (uint  i = 0 ; i < tElementNodes.numel() ; i++ )
            {
                moris::Matrix<moris::DDRMat> tmp =  tInterpMesh->get_node_coordinate( tElementNodes( i ) );
                tCoordinates.get_row(i) = tmp.get_row(0);
            }

            for( uint i = 0 ; i < tNumSurfaceNodes ; i++)
            {
                //form the parametric coordinates of the nodes
                moris::Matrix<DDRMat> tParamCoordsRelativeToElem = { { +1.0, tUniqueIntersectedPoints( 0, i ) } };

                //covert to the desired form
                tParamCoordsRelativeToElem = tParamCoordsRelativeToElem * tRotation;

                //obtain physical coordinates of the node
                Matrix <DDRMat > tNewCoords;
                moris::mtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem, tNewCoords);

                moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE);

                //construct the new vertex
                tVertcies ( i ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewCoords ) ;

                //update the id and index
                mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );
                tVertexId++;
            }

            // center node coordinates
            moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
            Matrix <DDRMat > tParamCoordsRelativeToElem = {{0.0, 0.0}};
            moris::mtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem, tNewNodeCoordinates);

            //update the id and index
            moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;
            mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );

            //add the top vertex of TRI3
            tVertcies( tNumSurfaceNodes ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewNodeCoordinates) ;

            return tVertcies;
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect_2D::rotation_matrix(moris::Matrix< DDRMat > &aRotation, moris::Matrix< DDRMat > &aInverseRotation, uint aPairCount )
        {
            if( mMeshSideSetPairs(aPairCount)(0) == "4")
            {
                aRotation        = { { 1, 0}, { 0, 1} };
                aInverseRotation = { { 1, 0}, { 0, 1}  };
            }
            else if ( mMeshSideSetPairs(aPairCount)(0) == "1" )
            {
                aRotation        = { { 0,1 }, { 1, 0 } };
                aInverseRotation = { { 0,1 }, { 1, 0 } };
            }
        }
    } /* end namespace mtk */

} /* end namespace moris */
