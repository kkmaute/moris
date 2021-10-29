/*
 * cl_MTK_Intersection_Detect.cpp
 *
 *  Created on: Jun 7, 2021
 *      Author: momo
 */


#include "cl_MTK_Intersection_Detect.hpp"

#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

#include "cl_Interpolation.hpp"

#include<unordered_map>
//LINALG/src
#include "fn_unique.hpp"
#include "fn_norm.hpp"
#include "fn_sum.hpp"
#include "fn_inv.hpp"
#include "fn_dot.hpp"
#include "fn_sort.hpp"
#include "op_equal_equal.hpp"
#include "op_less.hpp"
#include "fn_find.hpp"
#include "fn_join_horiz.hpp"
#include "fn_trans.hpp"
// New Vertex, Cell and Side Cluster
#include "cl_MTK_Cell_ISC.hpp"
#include "cl_MTK_Vertex_ISC_Impl.hpp"
#include "cl_MTK_Side_Cluster_ISC_Impl.hpp"
// Includes to make double sided sets and new cells
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_Tracer.hpp"
#include "cl_Stopwatch.hpp" //CHR/src


#include "cl_MTK_Writer_Exodus.hpp"

namespace moris
{
    namespace mtk
    {

        Intersection_Detect::Intersection_Detect(
                std::shared_ptr<moris::mtk::Mesh_Manager> aMeshManager,
                moris::moris_index                        aMeshIndex,
                moris::ParameterList &                    aParameterList,
                moris::uint                               aNumBulkPhases): mMeshManager(aMeshManager),
                        mMeshIndex(aMeshIndex),
                        mNumBulkPhases(aNumBulkPhases),
                        mEntityLocaltoGlobalMap(4),
                        mNewNodeCoords(0)
                        {
            //get the periodic mesh set names
            std::string tMeshSideSetNames = aParameterList.get< std::string >( "periodic_side_set_pair" );

            //store the names in mMeshSideSetPairs
            string_to_cell_of_cell( aParameterList.get<std::string >( "periodic_side_set_pair" ), mMeshSideSetPairs );

            //Data to store Id's and Indices of the added Vertices and Cells
            mIntersectedMeshData.set_up_external_entity_data( mMeshManager->get_integration_mesh( mMeshIndex ) );

                        }

        Intersection_Detect::~Intersection_Detect()
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
        Intersection_Detect::perform()
        {

            // Initialize phase interaction table
            Matrix<IndexMat> tPhaseInteractionTable( mNumBulkPhases, mNumBulkPhases );

            for( uint Ii = 0; Ii < mNumBulkPhases ; Ii++ )
            {
                for(uint Ij = 0 ; Ij < mNumBulkPhases ; Ij++)
                {
                    tPhaseInteractionTable( Ii, Ij )= Ij + Ii *  mNumBulkPhases;
                }
            }

            //loop over the side set pairs
            for( uint tPairCount = 0 ; tPairCount < mMeshSideSetPairs.size() ; tPairCount++ )
            {
                //Keep track of time
                Tracer tTracer("MTK", "No Type", "Searching For Pairs");

                // permutation pair indicates the indicial relation between 3D and surface coordinates
                // namely depending on the x-y-z orientation , which indices must be picked
                std::pair< moris_id, moris_id > tPermutationOrder = this->permutation_pair( tPairCount );

                //obtain the integration mesh
                moris::mtk::Integration_Mesh* tIntegrationMesh = mMeshManager->get_integration_mesh( mMeshIndex );

                //prerequisite offset vector and sets on each periodic surface
                moris::Matrix<DDRMat >     tOffsetVector;
                moris::Cell< std::string > tFirstSideSetNames;
                moris::Cell< std::string > tSecondSideSetNames;

                //find the offset of two surfaces
                this->offset_vector( tOffsetVector, tFirstSideSetNames, tSecondSideSetNames, tPairCount );

                //initialize all the clusters on each side
                moris::Cell< mtk::Cluster const * > tSideClusters1;
                moris::Cell< mtk::Cluster const * > tSideClusters2;

                //cells on the side set are locally indexed and they are assigned a phase
                moris::Cell< moris_index > tLocalCellNumToColor1 ;

                //loop over the
                for (uint i = 0; i < tFirstSideSetNames.size() ; i++)
                {
                    //get the set on the first side
                    moris::mtk::Set* tSet = tIntegrationMesh->get_set_by_name( tFirstSideSetNames(i) );

                    //get color(phase) associated with the set
                    Matrix<IndexMat> tColor = tSet->get_set_colors();

                    // get clusters in the second set
                    moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet->get_clusters_on_set();

                    //Assign phases to the clusters
                    tLocalCellNumToColor1.resize( tLocalCellNumToColor1.size()+ tSetClusters.size(), tColor( 0 ) );

                    //append the clusters the list of the first side clusters
                    tSideClusters1.append( tSetClusters );
                }

                //clusters on the side set are locally indexed and they are assigned a phase to each one
                moris::Cell< moris_index > tLocalCellNumToColor2 ;

                //loop over the existing side sets on the 2nd side
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
                    tSideClusters2.append( tSetClusters );
                }

                //Initialize the identifier matrix
                moris::Matrix< IndexMat > tSideClusterIdentifier1;
                moris::Matrix< IndexMat > tSideClusterIdentifier2;

                //Initialize the background cell to side cluster map
                std::unordered_map< moris::moris_index, moris::Cell<moris_index> > tBackgroundCellToSideClusterMap1;
                std::unordered_map< moris::moris_index, moris::Cell<moris_index> > tBackgroundCellToSideClusterMap2;

                //populate the map and identifier
                this->generate_identifier( tSideClusters1, tPairCount, tSideClusterIdentifier1, tBackgroundCellToSideClusterMap1) ;
                this->generate_identifier( tSideClusters2, tPairCount, tSideClusterIdentifier2, tBackgroundCellToSideClusterMap2) ;

                //loop over the background cells in order to cut the matching pairs
                for ( const auto & tIter : tBackgroundCellToSideClusterMap1)
                {

                    //Find all the side clusters that lie in the left and right pair
                    moris::Cell<moris_index> const & tSideClustersIndices1  =  tIter.second;
                    moris::Cell<moris_index> const & tSideClustersIndices2  =  tBackgroundCellToSideClusterMap2[tIter.first];

                    //Consider different cases to minmize the effort of polygon clipping
                    //case 1 where both sides just have one cluster ( they can be trivial or non trivial)
                    if ( tSideClustersIndices1.size() == 1 and tSideClustersIndices2.size() == 1 )
                    {
                        //parametric coordinates are predefined as the corners and center 
                        moris::Cell< Matrix<DDRMat> > tParamCoordsCell2 (4, Matrix<DDRMat>(2,3) );
                        tParamCoordsCell2(0) = { { -1.0, +1.0,  0.0 }, { -1.0, -1.0,  0.0 } };
                        tParamCoordsCell2(1) = { { -1.0,  0.0, -1.0 }, { -1.0,  0.0, +1.0 } };
                        tParamCoordsCell2(2) = { { -1.0,  0.0,  1.0 }, { +1.0,  0.0,  1.0 } };
                        tParamCoordsCell2(3) = { {  1.0,  0.0, +1.0 }, {  1.0,  0.0, -1.0 } };

                        // extract the phase of each side cluster
                        moris_index tPhase2  = tLocalCellNumToColor2 ( tSideClustersIndices2(0) );
                        moris_index tPhase1  = tLocalCellNumToColor1 (tSideClustersIndices1(0) );

                        //phase interaction table value 
                        moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                        //create a dummy cell in order to specify that all of the cut cells need to be assigned
                        moris::Cell<moris_index> tDummyCel= {0,1,2,3};

                        // create the dbl sided cluster based on the cut cell and IP cell information
                        this->create_dbl_sided_cluster(tParamCoordsCell2,tDummyCel,
                                tSideClusters1(tSideClustersIndices1(0))->get_interpolation_cell(),
                                tSideClusters2(tSideClustersIndices2(0))->get_interpolation_cell(),
                                tPairCount,tPhaseToPhaseIndex);
                    }

                    //case 2 where left side has one cluster(trivial/non-trivial) and right has more than one
                    else if ( tSideClustersIndices1.size() == 1 and tSideClustersIndices2.size() > 1 )
                    {

                        // iterate through the clusters on the right in order to create dbl sided cluster separately
                        for( uint iCluster = 0 ; iCluster < tSideClustersIndices2.size() ; iCluster++ )
                        {

                            // obtain numbder of surfaces present in each side cluster
                            uint tPrimaryCellInClusterNum = tSideClusters2( tSideClustersIndices2(iCluster) )->get_num_primary_cells();

                            // Initialize the cell containing coordinates of the cut surfaces
                            moris::Cell < moris::Matrix <DDRMat> > tParamCoordsCell (tPrimaryCellInClusterNum, moris::Matrix <DDRMat>(2,3) );

                            // Iterate through each surface to get the coordinates on the side ordinal
                            for(uint iCell = 0 ; iCell <  tPrimaryCellInClusterNum ; iCell++ )
                            {
                                // get the coordinates on the side ordinal and transpose
                                moris::Matrix <DDRMat>  t3DcoordMatrix = trans ( tSideClusters2( tSideClustersIndices2(iCluster) )->get_cell_local_coords_on_side_wrt_interp_cell(iCell) );

                                //convert 3D coordinates to suraface coordinates
                                moris::Matrix <DDRMat> tSurfaceCoordMatrix(2,3);
                                tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder.first  );
                                tSurfaceCoordMatrix.get_row( 1 ) = t3DcoordMatrix.get_row(  tPermutationOrder.second );

                                // sort them counter clock wise in order to be treated the same as other coordinates
                                this->SortandRemove(tSurfaceCoordMatrix );

                                //add the surface trinagle to the list of coordinates
                                tParamCoordsCell( iCell ) =  tSurfaceCoordMatrix ;
                            }

                            //create a dummy cell to indicate all the triangle elements will be processed to create dbl sided set
                            moris::Cell<moris_index> tDummyCell(tPrimaryCellInClusterNum);

                            //fill in the dummy cell with dummy values to include all the tri elements
                            for( uint i = 0 ; i <tParamCoordsCell.size(); i++ )
                            {
                                tDummyCell(i) = i ;
                            }

                            // get the phase information of the clusters
                            moris_index tPhase2  = tLocalCellNumToColor2 ( tSideClustersIndices2(iCluster) );
                            moris_index tPhase1  = tLocalCellNumToColor1 (tSideClustersIndices1(0) );

                            //phase interaction table value
                            moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                            // create the dbl sided cluster
                            this->create_dbl_sided_cluster(tParamCoordsCell, tDummyCell,
                                    tSideClusters1(tSideClustersIndices1(0))->get_interpolation_cell(),
                                    tSideClusters2(tSideClustersIndices2(iCluster))->get_interpolation_cell(),
                                    tPairCount,tPhaseToPhaseIndex );

                        }
                    }

                    //case 3 where right side has one cluster(trivial/non-trivial) and left has more than one
                    else if ( tSideClustersIndices1.size() > 1 and tSideClustersIndices2.size() == 1 )
                    {

                        // iterate through the clusters on the right in order to create dbl sided cluster separately
                        for( uint iCluster = 0 ; iCluster < tSideClustersIndices1.size() ; iCluster++ )
                        {

                            // obtain numbder of surfaces present in each side cluster
                            uint tPrimaryCellInClusterNum = tSideClusters1( tSideClustersIndices1(iCluster) )->get_num_primary_cells();

                            // Initialize the cell containing coordinates of the cut surfaces
                            moris::Cell < moris::Matrix <DDRMat> > tParamCoordsCell (tPrimaryCellInClusterNum, moris::Matrix <DDRMat>(2,3) );

                            // Iterate through each surface to get the coordinates on the side ordinal
                            for(uint iCell = 0 ; iCell <  tPrimaryCellInClusterNum ; iCell++ )
                            {
                                // get the coordinates on the side ordinal and transpose
                                moris::Matrix <DDRMat>  t3DcoordMatrix = trans ( tSideClusters1( tSideClustersIndices1(iCluster) )->get_cell_local_coords_on_side_wrt_interp_cell(iCell) );

                                //convert 3D coordinates to suraface coordinates
                                moris::Matrix <DDRMat> tSurfaceCoordMatrix(2,3);
                                tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder.first  );
                                tSurfaceCoordMatrix.get_row( 1 ) = t3DcoordMatrix.get_row(  tPermutationOrder.second );

                                // sort them counter clock wise in order to be treated the same as other coordinates
                                this->SortandRemove(tSurfaceCoordMatrix );

                                //add the surface trinagle to the list of coordinates
                                tParamCoordsCell( iCell ) =  tSurfaceCoordMatrix ;
                            }

                            //create a dummy cell to indicate all the triangle elements will be processed to create dbl sided set
                            moris::Cell<moris_index> tDummyCell(tPrimaryCellInClusterNum);

                            //fill in the dummy cell with dummy values to include all the tri elements
                            for( uint i = 0 ; i <tParamCoordsCell.size(); i++ )
                            {
                                tDummyCell(i) = i ;
                            }

                            // get the phase information of the each side cluster on the sides
                            moris_index tPhase2  = tLocalCellNumToColor2 ( tSideClustersIndices2(0) );
                            moris_index tPhase1  = tLocalCellNumToColor1 (tSideClustersIndices1(iCluster) );

                            //phase interaction table value
                            moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                            // iterate through the clusters on the right in order to create dbl sided cluster separately
                            this->create_dbl_sided_cluster(tParamCoordsCell, tDummyCell,
                                    tSideClusters1(tSideClustersIndices1(iCluster))->get_interpolation_cell(),
                                    tSideClusters2(tSideClustersIndices2(0))->get_interpolation_cell(),
                                    tPairCount,tPhaseToPhaseIndex );
                        }
                    }

                    //case 4 where both sides have more than one cluster within
                    else
                    {
                        uint tIGCellNum = 0; 
                        //loop through the clusters to obtain number of IG surfaces of each cluster
                        for(uint iCluster = 0 ; iCluster < tSideClustersIndices1.size(); iCluster++)
                        {
                            tIGCellNum  += tSideClusters1( tSideClustersIndices1( iCluster ) )->get_num_primary_cells();
                        }

                        // Initialize the cut cell surfaces 
                        moris::Cell < moris::Matrix <DDRMat> > tParamCoordsCell1 (tIGCellNum, moris::Matrix <DDRMat>(2,3) );

                        // A local map determining that each IG cell/surfaces belongs to which local cluster
                        moris::Matrix< moris::IndexMat> tIGCellToSideClusterMap1( tIGCellNum , 1 ) ;

                        //iterative counter to assign the cut cell surfaces and IG cell to local cluster index map
                        uint iCounter = 0;

                        //loop through each cluster to assign coordinates of the surface 
                        for(uint iCluster = 0 ; iCluster < tSideClustersIndices1.size(); iCluster++)
                        {
                            // number of surfaces in each cluster
                            uint tPrimaryCellInClusterNum = tSideClusters1( tSideClustersIndices1(iCluster) )->get_num_primary_cells();

                            //loop through IG cells in the cluster to get their coordinates
                            for(uint iCell = 0 ; iCell < tPrimaryCellInClusterNum ; iCell++ )
                            {
                                // get the coordinates on the side ordinal and transpose
                                moris::Matrix <DDRMat>  t3DcoordMatrix = trans ( tSideClusters1( tSideClustersIndices1(iCluster) )->get_cell_local_coords_on_side_wrt_interp_cell(iCell) );

                                //convert 3D coordinates to suraface coordinates
                                moris::Matrix <DDRMat> tSurfaceCoordMatrix(2,3);
                                tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder.first  );
                                tSurfaceCoordMatrix.get_row( 1 ) = t3DcoordMatrix.get_row(  tPermutationOrder.second );

                                //add the surface trinagle to the list of coordinates
                                tParamCoordsCell1( iCounter ) = tSurfaceCoordMatrix;

                                //add the local side cluster index of the IG cell which it belongs to
                                tIGCellToSideClusterMap1( iCounter ) = iCluster;

                                //increase the counter by 1
                                iCounter++;
                            }
                        }

                        //initialize total number of IG cells/surfaces on side 2
                        tIGCellNum = 0 ;

                        //loop through the clusters to obtain number of IG surfaces of each cluster
                        for(uint iCluster = 0 ; iCluster < tSideClustersIndices2.size(); iCluster++)
                        {
                            tIGCellNum += tSideClusters2(tSideClustersIndices2(iCluster))->get_num_primary_cells();
                        }

                        // Initialize the cell of surfaces that will be intersected
                        moris::Cell < moris::Matrix <DDRMat> > tParamCoordsCell2 (tIGCellNum, moris::Matrix <DDRMat>(2,3) );

                        // A local map determining that each IG cell/surfaces belongs to which local cluster
                        moris::Matrix< moris::IndexMat> tIGCellToSideClusterMap2( tIGCellNum, 1) ;

                        //iterative counter to assign the cut cell surfaces and IG cell to local cluster index map
                        iCounter = 0;

                        //loop through each cluster to assign coordinates of the surface 
                        for(uint iCluster = 0 ; iCluster < tSideClustersIndices2.size(); iCluster++)
                        {
                            // number of surfaces in each cluster
                            uint tPrimaryCellInClusterNum = tSideClusters2( tSideClustersIndices2(iCluster) )->get_num_primary_cells();

                            //loop through IG cells in the cluster to get their coordinates
                            for(uint iCell = 0 ; iCell < tPrimaryCellInClusterNum ; iCell++ )
                            {
                                // get the coordinates on the side ordinal and transpose
                                moris::Matrix <DDRMat>  t3DcoordMatrix = trans ( tSideClusters2( tSideClustersIndices2(iCluster) )->get_cell_local_coords_on_side_wrt_interp_cell(iCell) );

                                //convert 3D coordinates to suraface coordinates
                                moris::Matrix <DDRMat> tSurfaceCoordMatrix(2,3);
                                tSurfaceCoordMatrix.get_row( 0 ) = t3DcoordMatrix.get_row( tPermutationOrder.first  );
                                tSurfaceCoordMatrix.get_row( 1 ) = t3DcoordMatrix.get_row(  tPermutationOrder.second );

                                //add the surface trinagle to the list of coordinates
                                tParamCoordsCell2( iCounter ) = tSurfaceCoordMatrix;

                                //add the local side cluster index of the IG cell which it belongs to
                                tIGCellToSideClusterMap2( iCounter ) = iCluster;

                                //increase the counter by 1
                                iCounter++;
                            }
                        }

                        //initialize the cut polygons and their respective indices that come from their parent element
                        moris::Cell < moris::Matrix <DDRMat> > tCutPolygons;
                        moris::Matrix< moris::IndexMat> tCutPolygonIdentifier;

                        // Polygon clipping algorthim
                        this->elementwise_bruteforce_search(
                                tParamCoordsCell1,tIGCellToSideClusterMap1,
                                tParamCoordsCell2,tIGCellToSideClusterMap2,
                                tCutPolygons, tCutPolygonIdentifier);

                        //a map from the identifier of each cut cell to all cut cells with the same identifier
                        std::unordered_map< moris_index , moris::Cell< moris_index > > tCutCellIdentifierToCutCell;

                        //populate the map
                        this->group_cut_cells( tCutPolygonIdentifier, tCutCellIdentifierToCutCell );

                        //unique multiplier to distingush between sub clusters
                        uint tMultiplier = std::max(   tParamCoordsCell1.size(),   tParamCoordsCell2.size() );

                        // iterate through grouped cut cell in order to create dbl sideded cluster with the correct 
                        //IP cell and phase information
                        for( auto & iCutCellGroup : tCutCellIdentifierToCutCell)
                        {
                            // indices of the cut polygon that exist in the particular dbl sided cluster
                            moris::Cell<moris_index> & tIndicesinCutCell  =  iCutCellGroup.second;

                            //local index of the side cluster on the right
                            uint j = (iCutCellGroup.first) % tMultiplier;

                            //local index of the side cluster on the right
                            uint i = std::floor (iCutCellGroup.first / tMultiplier) ;

                            //phase information of the side clusters
                            moris_index tPhase2  = tLocalCellNumToColor2 ( tSideClustersIndices2(j) );
                            moris_index tPhase1  = tLocalCellNumToColor1 ( tSideClustersIndices1(i) );

                            //phase interaction table value
                            moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                            //create the dbl sided cluster
                            this->create_dbl_sided_cluster(tCutPolygons,  tIndicesinCutCell,
                                    tSideClusters1( tSideClustersIndices1(i))->get_interpolation_cell(), tSideClusters2( tSideClustersIndices2(j))->get_interpolation_cell(),
                                    tPairCount, tPhaseToPhaseIndex) ;
                        }
                        //Find the intersection of two surfaces along with the relative IP cells
                    }
                }
            }

            //add the double sided set based on the interactions
            this->construct_add_dbl_sided_set(tPhaseInteractionTable);
        }


        //---------------------------------------------------------------------

        void
        Intersection_Detect::EdgeIntersect(
                moris::Matrix < DDRMat > const & aFirstTRICoords,
                moris::Matrix < DDRMat > const &  aSecondTRICoords,
                moris::Matrix < DDRMat > & aIntersectedPoints,
                moris::Matrix < DDUMat > & aIntersectVec) const
        {
            uint tNumIntersections = 0;

            aIntersectVec.set_size(1,3,0);

            for( uint i = 0 ; i < 3 ; i++ )
            {
                for( uint j = 0 ; j < 3 ; j++ )
                {
                    //form the matrix on RHS
                    moris::Matrix < DDRMat > b = aSecondTRICoords.get_column( j ) - aFirstTRICoords.get_column( i ) ;

                    //matrix on the LHS
                    moris::Matrix < DDRMat > A( 2,2 );
                    A.set_column( 0, aFirstTRICoords.get_column( (i+1) % 3  ) - aFirstTRICoords.get_column( i ) )  ;
                    A.set_column( 1, - aSecondTRICoords.get_column( (j+1) % 3 ) + aSecondTRICoords.get_column( j ) )   ;

                    //solve the system
                    if ( std::abs( A(0,0)*A(1,1) - A(0,1)*A(1,0) ) > 0.0000001 )
                    {
                        moris::Matrix < DDRMat > R = inv(A) * b;

                        moris::real eps = 0.01;

                        //Intersection Condition
                        if ( R(0) >= 0 and ( R(0) - 1 ) <= eps and R(1) >=0 and ( R(1) - 1) <= eps)
                        {
                            //grow the matrix to insert the new intersection
                            aIntersectedPoints.resize(2,tNumIntersections+1);

                            aIntersectedPoints.get_column(tNumIntersections) = aFirstTRICoords.get_column( i ) + R(0) * A.get_column( 0 ) ;

                            aIntersectVec( i ) = 1 ;

                            //increase number of intersection by 1
                            tNumIntersections++ ;
                        }
                    }
                }
            }

        }

        //-----------------------------------------------------------------------------------------

        void  Intersection_Detect::PointsXInY(
                moris::Matrix < moris::DDRMat >  const &  aFirstTRICoords,
                moris::Matrix < moris::DDRMat  > const &  aSecondTRICoords,
                moris::Matrix < moris::DDRMat  >       &  aIntersectedPoints) const
        {
            uint tNumIntersections = 0;

            //Interior points
            moris::Matrix < moris::DDRMat  > v0 = aSecondTRICoords.get_column(1) - aSecondTRICoords.get_column(0);
            moris::Matrix < moris::DDRMat  > v1 = aSecondTRICoords.get_column(2) - aSecondTRICoords.get_column(0);

            //Baricenteric Coordinates
            real d00 = dot( v0, v0 );
            real d01 = dot( v0, v1 );
            real d11 = dot( v1, v1 );

            real id = 1.0 / ( d00 * d11 - d01* d01 );

            for( uint i = 0 ; i < 3 ; i++)
            {
                moris::Matrix < moris::DDRMat  > v2 = aFirstTRICoords.get_column(i) - aSecondTRICoords.get_column(0);

                real d02 = dot( v0, v2 );
                real d12 = dot( v1, v2 );

                real u = ( d11 * d02 - d01 * d12 ) * id;
                real v = ( d00 * d12 - d01 * d02 ) * id;

                real eps = 0.001;

                if ( u >= -eps and v >= -eps and (u + v <= 1+eps) )
                {
                    aIntersectedPoints.resize( 2, tNumIntersections + 1 );

                    aIntersectedPoints.get_column( tNumIntersections ) = aFirstTRICoords.get_column(i) ;

                    tNumIntersections++ ;
                }

            }
        }

        //------------------------------------------------------------------------

        void  Intersection_Detect::SortandRemove
        ( moris::Matrix < moris::DDRMat  > & aIntersectedPoints ) const
        {
            real eps = 0.0001 ;

            uint tNumColumn = aIntersectedPoints.n_cols();

            if ( tNumColumn > 0 )
            {
                // C
                moris::Matrix < moris::DDRMat  > C(2,1);
                C.get_row(0) = sum( aIntersectedPoints.get_row( 0 ) )/tNumColumn;
                C.get_row(1) = sum( aIntersectedPoints.get_row( 1 ) )/tNumColumn;

                moris::Matrix < moris::DDRMat  > ao( 1, tNumColumn );

                //order polygon corners counter
                for( uint i = 0 ; i < tNumColumn ; i++ )
                {
                    moris::Matrix < moris::DDRMat  > d = aIntersectedPoints.get_column( i ) - C;

                    ao(i) = std::atan2( d( 1 ) , d( 0 ) );
                }

                moris::Matrix < moris::DDRMat  > aoSorted;

                moris::Matrix < moris::DDUMat  > IdMatrix( 1, tNumColumn );

                moris::sort( ao, aoSorted ,"ascend", 1 );

                //find the index matrix
                for( uint i = 0 ; i < tNumColumn; i++ )
                {
                    for (uint j = 0 ; j < tNumColumn ; j++ )
                    {
                        if( ao( j ) == aoSorted( i ) )
                        {
                            IdMatrix( i ) = j ;
                            break;
                        }
                    }
                }

                //sort the points based on IdMatrix
                moris::Matrix < moris::DDRMat  > tSortedIntersectedPoints(2,tNumColumn);

                for( uint i = 0; i < tNumColumn ; i++ )
                {
                    tSortedIntersectedPoints.get_column( i ) = aIntersectedPoints.get_column( IdMatrix( i ) );
                }

                //remove duplicates of the point
                uint i = 0;
                uint j = 1;

                while( j < tNumColumn )
                {
                    if ( norm( tSortedIntersectedPoints.get_column(i) -  tSortedIntersectedPoints.get_column(j)  ) > eps )
                    {
                        i++;

                        tSortedIntersectedPoints.get_column(i) = tSortedIntersectedPoints.get_column(j) ;

                        j++;

                    }
                    else
                    {
                        j++;
                    }
                }

                tSortedIntersectedPoints.resize(2,i+1);

                aIntersectedPoints = tSortedIntersectedPoints;
            }
        }

        //------------------------------------------------------------------------

        void
        Intersection_Detect::Intersect(
                moris::Matrix < moris::DDRMat  >  const & aFirstTRICoords,
                moris::Matrix < moris::DDRMat  > const &  aSecondTRICoords,
                moris::Matrix < moris::DDRMat  > & aIntersectedPoints,
                moris::Matrix < moris::DDUMat  > & aIntersectVec) const
        {
            this->EdgeIntersect(aFirstTRICoords,aSecondTRICoords,aIntersectedPoints, aIntersectVec);

            Matrix < DDRMat > P1;

            this->PointsXInY( aFirstTRICoords, aSecondTRICoords, P1 ) ;


            if ( P1.n_cols() > 1)
            {
                aIntersectVec.set_size( 1, 3, 1 );
            }

            aIntersectedPoints  = join_horiz( aIntersectedPoints, P1) ;

            this->PointsXInY( aSecondTRICoords, aFirstTRICoords, P1 ) ;

            aIntersectedPoints  = join_horiz( aIntersectedPoints, P1) ;

            this->SortandRemove (aIntersectedPoints ) ;
        }

        //------------------------------------------------------------------------------
        void
        Intersection_Detect::elementwise_bruteforce_search (
                moris::Cell< Matrix<moris::DDRMat > > const & aParamCoordsCell1,
                moris::Matrix<IndexMat > const & aIGCellToSideClusterMap1,
                moris::Cell< Matrix<moris::DDRMat > > const  & aParamCoordsCell2,
                moris::Matrix<IndexMat > const & aIGCellToSideClusterMap2,
                moris::Cell< Matrix<moris::DDRMat > >  & aIntersectedAreas,
                moris::Matrix < IndexMat > & aIntersectedAreasIdentifier) const
        {
            // multiplier used to assign an id to each cut surfaces based on the parent cells
            uint tMultiplier = std::max(aParamCoordsCell1.size(), aParamCoordsCell2.size() );

            // max size of the cut  polygons
            uint tMaxSize = aParamCoordsCell1.size() * aParamCoordsCell2.size() ;

            //reserve the space for max 
            aIntersectedAreas.reserve( tMaxSize );

            //set the max size in order to avoid resizing at each step
            aIntersectedAreasIdentifier.set_size( tMaxSize, 1 );

            // Initialize the output
            moris::Matrix < moris::DDUMat  >  tnc;

            //counter to assign for cut polygon and their identifier
            uint iCounter = 0 ;

            // Loop over second mesh
            for( uint iI = 0 ; iI < aParamCoordsCell1.size(); iI++ )
            {
                // Loop over the first mesh
                for( uint iJ = 0 ; iJ < aParamCoordsCell2.size(); iJ++)
                {
                    // initialize matrix as it needs to be refilled
                    moris::Matrix < moris::DDRMat  > tP;

                    // find the intersection of 2 element triangulation
                    this->Intersect( aParamCoordsCell1(iI), aParamCoordsCell2(iJ), tP ,tnc) ;

                    // If it is a polygon add to the output
                    if (tP.n_cols()>2)
                    {
                        //add the cut polygon to the cell and assign identifier to it
                        aIntersectedAreas.push_back(tP);
                        aIntersectedAreasIdentifier(iCounter) = aIGCellToSideClusterMap1(iI)*tMultiplier+ aIGCellToSideClusterMap2(iJ);

                        //increase the counter by 1
                        iCounter++;
                    }
                }
            }

            //resize the cell and matrix to the actual size
            aIntersectedAreas.shrink_to_fit();
            aIntersectedAreasIdentifier.resize(iCounter,1);
        }

        //------------------------------------------------------------------------------

        void
        Intersection_Detect::create_dbl_sided_cluster( moris::Cell< Matrix<DDRMat> > & tP,
                moris::Cell<moris_index> & aIndicesinCutCell,
                moris::mtk::Cell const & aInterpCell1,
                moris::mtk::Cell const & aInterpCell2,
                uint & aPairCount,
                moris_index  & aPhaseToPhase)
        {
            //requires that we give a unique coordinates and associated index to each of them
            Matrix< DDRMat > tUniqueIntersectedPoints;

            //join all the coordinates together
            for (size_t i = 0 ; i < aIndicesinCutCell.size(); i++ )
            {
                tUniqueIntersectedPoints = join_horiz( tUniqueIntersectedPoints, tP( aIndicesinCutCell(i) ) );
            }

            //remove duplicates of the point
            uint i = 0;
            uint j = 1;
            uint tNumCol = tUniqueIntersectedPoints.n_cols() ;

            while( j < tNumCol )
            {
                if ( norm( tUniqueIntersectedPoints.get_column(i) -  tUniqueIntersectedPoints.get_column(j)  ) >  0.00000000001 )
                {
                    i++;

                    tUniqueIntersectedPoints.get_column(i) = tUniqueIntersectedPoints.get_column(j) ;

                    j++;

                }
                else
                {
                    j++;
                }
            }

            tUniqueIntersectedPoints.resize(2,i+1);  

            //intersection cell from the algorithm comes as coordinates, convert them to unique coordinate indices
            moris::Cell< Matrix< IndexMat > > tPVertexIndex( aIndicesinCutCell.size() , Matrix< IndexMat >(1,3)) ;



            //loop over the cell of coordinates to assign indicies
            for( size_t i = 0 ; i < aIndicesinCutCell.size() ; i++ )
            {
                if(tP( aIndicesinCutCell(i) ).n_cols() != 3 )
                {
                    tPVertexIndex(i) = Matrix<IndexMat>(1, tP( aIndicesinCutCell(i) ).n_cols());
                }

                for( size_t j = 0; j < tP( aIndicesinCutCell(i) ).n_cols(); j++ )
                {
                    for( size_t k = 0 ; k < tUniqueIntersectedPoints.n_cols() ; k++ )
                    {
                        if ( norm( tP( aIndicesinCutCell(i) ).get_column(j) - tUniqueIntersectedPoints.get_column(k) ) < 0.00000000001 )
                        {
                            tPVertexIndex(i)(j) = k ;
                            break;
                        }
                    }
                }
            }

            // get the IP mesh in order to interpolate for physical coordinates
            moris::mtk::Interpolation_Mesh* tInterpMesh = mMeshManager->get_interpolation_mesh( mMeshIndex );

            //get the indices of the interpolation cell of the first side nodes
            moris::Matrix<moris::IndexMat> tElementNodes =
                    tInterpMesh->get_entity_connected_to_entity_loc_inds(
                            aInterpCell1.get_index(),
                            moris::EntityRank::ELEMENT,
                            moris::EntityRank::NODE);

            // coordinates of nodes attached to element to interpolate for physical coordinates
            moris::Matrix<moris::DDRMat> tIPCellCoordinates( tElementNodes.numel() , 3);

            //extract the corner coordinates of the IP cell
            for (uint  i = 0 ; i < tElementNodes.numel() ; i++ )
            {
                tIPCellCoordinates.get_row(i)  =  tInterpMesh->get_node_coordinate( tElementNodes( i ) ).get_row(0);

            }

            //numbder of nodes at the surface
            uint tNumSurfaceNodes = tUniqueIntersectedPoints.n_cols();

            //paramteric coordinates master side 3D
            Matrix <DDRMat> tParametricCoordinates3D(tNumSurfaceNodes,3);

            //permutation of the 2D coordinates to 3D
            std::pair<moris_id, moris_id> tCoordOrder = this->permutation_pair(aPairCount);

            //3rd coordinate that 
            Matrix <DDRMat> tDummyCoords(tNumSurfaceNodes,1,-1);

            //fill in the coordinate matrix
            tParametricCoordinates3D.get_column(3 - tCoordOrder.first - tCoordOrder.second ) = tDummyCoords.get_column(0);
            tParametricCoordinates3D.get_column( tCoordOrder.first ) = trans( tUniqueIntersectedPoints.get_row(0) );
            tParametricCoordinates3D.get_column( tCoordOrder.second) = trans( tUniqueIntersectedPoints.get_row(1) );

            //initialize the 3D physical coordinates
            Matrix <DDRMat> tPhysicalCoordinates3D(tNumSurfaceNodes,3);

            //interpolate the physical coordinates
            moris::mtk::Interpolation::trilinear_interpolation_multivalue(tIPCellCoordinates, tParametricCoordinates3D, tPhysicalCoordinates3D);

            //head vertex of the tetrahedron which lies in the center of the cube
            Matrix< DDRMat > tCoordinatesTop;
            Matrix< DDRMat > tHeadnodeCoords(1,3,0);
            moris::mtk::Interpolation::trilinear_interpolation(tIPCellCoordinates, tHeadnodeCoords, tCoordinatesTop);

            //create master and slave vertices
            moris::Cell<moris::mtk::Vertex* > tMasterVerticesCell;
            this->create_master_vertices(tMasterVerticesCell, tPhysicalCoordinates3D,aInterpCell1, aPairCount,tCoordinatesTop);

            //Slave Side Vertices
            //get the indices of the interpolation cell nodes
            tElementNodes =
                    tInterpMesh->get_entity_connected_to_entity_loc_inds(
                            aInterpCell2.get_index(),
                            moris::EntityRank::ELEMENT,
                            moris::EntityRank::NODE);

            //extract the corner coordinates of the IP cell
            for (uint  i = 0 ; i < tElementNodes.numel() ; i++ )
            {
                tIPCellCoordinates.get_row(i)  =  tInterpMesh->get_node_coordinate( tElementNodes( i ) ).get_row(0);

            }

            //Duplicate the parametric coordinates for the salve side
            Matrix <DDRMat> tSlaveParametricCoordinates3D = tParametricCoordinates3D.copy();

            // Use the same matrices and change the 3rd dimension coordinates
            tDummyCoords.fill( 1.0 );
            tSlaveParametricCoordinates3D.get_column( 3 - tCoordOrder.first - tCoordOrder.second ) = tDummyCoords.get_column( 0 );

            // get physical coordinates of the slave side
            moris::mtk::Interpolation::trilinear_interpolation_multivalue( tIPCellCoordinates, tSlaveParametricCoordinates3D, tPhysicalCoordinates3D );

            // physcal coordinates of the 
            moris::mtk::Interpolation::trilinear_interpolation( tIPCellCoordinates, tHeadnodeCoords, tCoordinatesTop );

            moris::Cell<moris::mtk::Vertex* > tSlaveVerticesCell;
            this->create_slave_vertices(tSlaveVerticesCell, tPhysicalCoordinates3D,aInterpCell2, aPairCount,tCoordinatesTop);

            //create constant version of the vertices on cluster
            moris::Cell< moris::mtk::Vertex const *> tMasterVerticesConst( tMasterVerticesCell.size() - 1,nullptr);
            moris::Cell< moris::mtk::Vertex const *> tSlaveVerticesConst( tMasterVerticesCell.size() - 1,nullptr);

            for(uint i = 0 ; i < tMasterVerticesCell.size() - 1 ; i++)
            {
                tMasterVerticesConst( i ) = tMasterVerticesCell( i ) ;
                tSlaveVerticesConst( i ) = tSlaveVerticesCell( i ) ;
            }

            //store the vertices to prevent memory leak
            mMasterVertices.append(tMasterVerticesConst);
            mSlaveVertices.append(tSlaveVerticesConst);


            tic tTimerCell;

            //initialize master,slave side IG cells
            moris::Cell<moris::mtk::Cell const *> tMasterIntegCells;
            moris::Cell<moris::mtk::Cell const *> tSlaveIntegCells;

            //Added surfaces if intersection area is a polygon
            uint tAddedSurafceNum = 0 ;

            //iterate through each intersection surface and make Ig cells
            for( size_t tClusterNum = 0 ; tClusterNum < aIndicesinCutCell.size() ; tClusterNum++ )
            {
                //if polygon break it further to tri by fan triangulation
                if ( tPVertexIndex(tClusterNum).n_cols() > 3 )
                {
                    //find number of newly aded surfaces if it is a polygon
                    uint tNumEdges = tPVertexIndex(tClusterNum).n_cols();

                    uint tNumTriangles = tNumEdges - 2 ;

                    tAddedSurafceNum += ( tNumTriangles - 1 ) ;

                    //loop over the subtriangles created by fan triangulation
                    for(uint tTRI = 0 ; tTRI < tNumTriangles ; tTRI++ )
                    {
                        // vertex Index of the cretae sub triangles
                        Matrix< IndexMat > tTmpVertexIndex = { {tPVertexIndex(tClusterNum)(0),tPVertexIndex(tClusterNum)(tTRI+1),tPVertexIndex(tClusterNum)(tTRI+2) } };

                        //master IG cell
                        moris::mtk::Cell const * tMasterIgCell = this->create_master_ig_cell( tMasterVerticesCell, tTmpVertexIndex, aInterpCell1,aPairCount ) ;
                        tMasterIntegCells.push_back( tMasterIgCell );

                        //slave IG cell
                        moris::mtk::Cell const * tSlaveIgCell = this->create_slave_ig_cell( tSlaveVerticesCell, tTmpVertexIndex, aInterpCell2,aPairCount ) ;
                        tSlaveIntegCells.push_back( tSlaveIgCell );
                    }
                }
                else
                {
                    //master IG cells
                    moris::mtk::Cell const * tMasterIgCell = this->create_master_ig_cell( tMasterVerticesCell, tPVertexIndex(tClusterNum), aInterpCell1, aPairCount ) ;
                    tMasterIntegCells.push_back( tMasterIgCell );

                    //slave IG cells
                    moris::mtk::Cell const * tSlaveIgCell = this->create_slave_ig_cell( tSlaveVerticesCell, tPVertexIndex(tClusterNum), aInterpCell2 ,aPairCount) ;
                    tSlaveIntegCells.push_back( tSlaveIgCell );
                }
            }

            //create side ordinal matrices, the entry is always 3
            moris::Matrix<moris::IndexMat>          tMasterIntegrationCellSideOrdinals( 1, aIndicesinCutCell.size() + tAddedSurafceNum, 3);
            moris::Matrix<moris::IndexMat>          tSlaveIntegrationCellSideOrdinals( 1, aIndicesinCutCell.size() + tAddedSurafceNum, 3);

            //create the master side cluster
            moris::mtk::Side_Cluster_ISC* tMasterSideCluster = new moris::mtk::Side_Cluster_ISC(false,
                    & aInterpCell1,
                    tMasterIntegCells,
                    tMasterIntegrationCellSideOrdinals,
                    tMasterVerticesConst,
                    tParametricCoordinates3D);

            //create the slave side cluster
            moris::mtk::Side_Cluster_ISC* tSlaveSideCluster  = new moris::mtk::Side_Cluster_ISC(false,
                    & aInterpCell2,
                    tSlaveIntegCells,
                    tSlaveIntegrationCellSideOrdinals,
                    tSlaveVerticesConst,
                    tSlaveParametricCoordinates3D);

            // create double side cluster
            mtk::Double_Side_Cluster* tDblSideCluster = new mtk::Double_Side_Cluster(
                    tMasterSideCluster,
                    tSlaveSideCluster,
                    tSlaveSideCluster->mVerticesInCluster);

            //store data to be used further
            mDoubleSidedClusters.push_back(tDblSideCluster);
            mDoubleSidedClustersIndex.push_back(aPhaseToPhase);

            //master
            mMasterSidedClusters.push_back(tMasterSideCluster);

            //slave
            mSlaveSidedClusters.push_back(tSlaveSideCluster);

            // Offset due to addition of new mtk cells at each loop
            uint tOffsetCells  = mMasterSideCells.size();

            // Append the integration cells
            mMasterSideCells.append(tMasterIntegCells);
            mSlaveSideCells.append(tSlaveIntegCells);

            //Determine the size of the elements local to global entity map
            uint tNumExistingElem = mEntityLocaltoGlobalMap(0).size() ;
            uint tNewElem =  tMasterIntegCells.size() + tSlaveIntegCells.size()+ mMasterSideCells.size()*2 ;

            // Resize number of elements
            mEntityLocaltoGlobalMap(3).resize(tNumExistingElem + tNewElem);

            // If it is requested to ouput the intersection mesh
            if ( mVisFlag )
            {

                // Create the maps for cell
                for(uint Ii = 0 ; Ii < tMasterIntegCells.size() ; Ii++)
                {
                    mMasterCellIndextoCellMap[ tMasterIntegCells(Ii)->get_index() ] = Ii+ tOffsetCells ;

                    mEntityLocaltoGlobalMap(3) ( Ii+ tOffsetCells)  =  tMasterIntegCells(Ii)->get_id();
                }

                //Determine the size of the local to global map for nodes
                uint tNumExistingNodes = mEntityLocaltoGlobalMap(0).size();
                uint tNewNodes = tMasterVerticesCell.size()+ tSlaveVerticesCell.size()+ mNewNodeCoords.size();

                // Resize the local to global map for the vertices
                mEntityLocaltoGlobalMap(0).resize(tNumExistingNodes + tNewNodes);

                // Offset due to addition of vertices at each loop
                uint tMasterOffset  = mNewNodeCoords.size();
                uint tSlaveOffset  = mNewNodeCoords.size() + tMasterVerticesCell.size();

                // Assign ids and indices to the maps
                for( uint Ii = 0 ; Ii < tMasterVerticesCell.size() ; Ii++ )
                {
                    // Assign Coordinates
                    Matrix<DDRMat> tCoords = tMasterVerticesCell(Ii)->get_coords();
                    mNewNodeCoords.push_back(tCoords);

                    mEntityLocaltoGlobalMap(0)(Ii + tMasterOffset)  = tMasterVerticesCell(Ii)->get_id();

                }

                // Assign ids and indices to the maps
                for( uint Ii = 0 ; Ii < tSlaveVerticesCell.size() ; Ii++ )
                {
                    // Assign Coordinates
                    mNewNodeCoords.push_back(tSlaveVerticesCell(Ii)->get_coords());

                    mEntityLocaltoGlobalMap(0)(Ii + tSlaveOffset) = tSlaveVerticesCell(Ii)->get_id();
                }
            }
        }

        //------------------------------------------------------------------------------------------------------------

        moris::mtk::Cell const *
        Intersection_Detect::create_master_ig_cell( moris::Cell<moris::mtk::Vertex *> tMasterVertices, Matrix<IndexMat> tVertexIndex ,  moris::mtk::Cell const & aMasterInterpCell,uint aPairCount)
        {
            // Sort the nodes such that it is consistent with Exodus convention
            Matrix<IndexMat> tVertexIndexCounterClockWise ;
            if( mMeshSideSetPairs(aPairCount)(0) == "5")
            {
                tVertexIndexCounterClockWise = tVertexIndex ;
            }
            else
            {
                tVertexIndexCounterClockWise = { {tVertexIndex(0) , tVertexIndex(2), tVertexIndex(1) } } ;

            }

            //cell of Ig vertices
            moris::Cell<moris::mtk::Vertex *> tIgCellVertices(4);

            //surface vertices
            for(uint i = 0 ; i < tVertexIndex.n_cols() ; i++ )
            {
                tIgCellVertices( i ) =  tMasterVertices( tVertexIndexCounterClockWise ( i ) ) ;
            }

            //head node of tet4
            tIgCellVertices( 3 ) =  tMasterVertices( tMasterVertices.size() - 1 ) ;

            //allocate id and index
            moris_index tCellIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::ELEMENT) ;
            moris_id tCellId= mIntersectedMeshData.allocate_entity_ids_external_entity_data(1, EntityRank::ELEMENT);

            //create cell info and geometry
            mtk::Cell_Info_Factory tCellInfoFactory;
            std::shared_ptr<moris::mtk::Cell_Info> tLinearCellInfo = tCellInfoFactory.create_cell_info_sp(Geometry_Type::TET,mtk::Interpolation_Order::LINEAR);

            //create the IG Cell
            moris::mtk::Cell* tIgCell = new Cell_ISC(
                    tCellId,
                    tCellIndex,
                    aMasterInterpCell.get_owner(),
                    tLinearCellInfo,
                    tIgCellVertices);

            mIntersectedMeshData.update_first_available_index_external_data( tCellIndex+1, EntityRank::ELEMENT);

            return tIgCell;
        }

        //------------------------------------------------------------------------------------------------------------

        moris::mtk::Cell const *
        Intersection_Detect::create_slave_ig_cell( moris::Cell<moris::mtk::Vertex *> tSlaveVertices, Matrix<IndexMat> tVertexIndex ,  moris::mtk::Cell const & aSlaveInterpCell, uint aPairCount)
        {
            // Sort the nodes such that it is consistent with Exodus convention
            Matrix<IndexMat> tVertexIndexCounterClockWise ;

            if( mMeshSideSetPairs(aPairCount)(0) == "5")
            {
                tVertexIndexCounterClockWise = { {tVertexIndex(0) , tVertexIndex(2), tVertexIndex(1) } } ;
            }
            else
            {
                tVertexIndexCounterClockWise = tVertexIndex ;

            }

            //cell of Ig vertices
            moris::Cell<moris::mtk::Vertex *> tIgCellVertices(4);

            //surface vertices
            for(uint i = 0 ; i < tVertexIndex.n_cols() ; i++ )
            {
                tIgCellVertices( i ) = tSlaveVertices( tVertexIndexCounterClockWise ( i ) ) ;
            }

            //head node of tet4
            tIgCellVertices( 3 ) = tSlaveVertices( tSlaveVertices.size() - 1 ) ;

            //allocate id and index
            moris_index tCellIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::ELEMENT) ;
            moris_id tCellId = mIntersectedMeshData.allocate_entity_ids_external_entity_data( 1, EntityRank::ELEMENT);

            // linear cell info
            mtk::Cell_Info_Factory tCellInfoFactory;
            std::shared_ptr<moris::mtk::Cell_Info> tLinearCellInfo = tCellInfoFactory.create_cell_info_sp(Geometry_Type::TET,mtk::Interpolation_Order::LINEAR);

            //construct the IG cell
            moris::mtk::Cell* tIgCell = new Cell_ISC(
                    tCellId,
                    tCellIndex,
                    aSlaveInterpCell.get_owner(),
                    tLinearCellInfo,
                    tIgCellVertices);

            mIntersectedMeshData.update_first_available_index_external_data( tCellIndex+1, EntityRank::ELEMENT);

            return tIgCell;
        }


        //------------------------------------------------------------------------------------------------------------
        //name the cluster set
        void Intersection_Detect::construct_add_dbl_sided_set (moris::Matrix < IndexMat > const & tPhaseInteractionTable)
        {
            Tracer tTracer("MTK", "Double Sided Set", "Construction of Double Sided Sets");

            //initialize all the double sided sets
            mDoubleSideSets.resize( tPhaseInteractionTable.numel() );
            mMasterSideSets.resize( tPhaseInteractionTable.numel() );

            //get the integration mesh
            moris::mtk::Integration_Mesh*  tIntegMesh  = mMeshManager->get_integration_mesh( mMeshIndex );

            uint tSpatialDim = tIntegMesh->get_spatial_dim();

            //assign double sided clusters to the right double sided set
            for( uint i = 0 ; i < mDoubleSidedClusters.size(); i++ )
            {
                mDoubleSideSets(mDoubleSidedClustersIndex(i)).push_back(mDoubleSidedClusters(i));

                mMasterSideSets(mDoubleSidedClustersIndex(i)).push_back(mMasterSidedClusters(i));
            }

            //add double sided set to the mesh and name them
            for( size_t aBulkPhase1 = 0 ; aBulkPhase1 <  tPhaseInteractionTable.n_cols() ; aBulkPhase1++)
            {
                for( size_t aBulkPhase2 = 0 ; aBulkPhase2 <  tPhaseInteractionTable.n_cols() ; aBulkPhase2++ )
                {
                    uint tPhaseToPhaseIndex = tPhaseInteractionTable( aBulkPhase1, aBulkPhase2 );

                    //name the double sided set
                    std::string  tDoubleSideSetName     = "P"+std::to_string(aBulkPhase1)+std::to_string(aBulkPhase2);

                    //check if double sided set is empty
                    if(  mDoubleSideSets(tPhaseToPhaseIndex).size() > 0)
                    {
                        //colors of the double sided set,  needs to be modified
                        moris::Matrix<moris::IndexMat>      tColors         = {{0}};

                        // Construct the double side set
                        moris::mtk::Set* tDblSideSet = new moris::mtk::Double_Side_Set(tDoubleSideSetName,  mDoubleSideSets(tPhaseToPhaseIndex), tColors, tSpatialDim);

                        // add double sided periodic boundary condition to the integration mesh
                        tIntegMesh->add_double_side_set( tDblSideSet );
                    }

                    //add it to be stored in the intersection mesh
                    mSideSideSetLabelToOrd[tDoubleSideSetName] = tPhaseToPhaseIndex ;

                    //add double sided name to the list of labels
                    mSideSetLabels.push_back(tDoubleSideSetName);

                }
            }
            //assign an name and index to the bulk elements (cells on the master side)
            mBlockSetLabelToOrd["Dummy"] = 0;
            mBlockSetLabels.push_back("Dummy");
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect::offset_vector(moris::Matrix<DDRMat > &tOffsetVector,  moris::Cell< std::string > & tFirstSideSetNames,moris::Cell< std::string > &tSecondSideSetNames,uint aPairCount )
        {
            // Integration mesh
            moris::mtk::Integration_Mesh*  tIntegrationMesh  = mMeshManager->get_integration_mesh( mMeshIndex );

            //get all the sidesets
            moris::Cell<std::string> tSetNames = tIntegrationMesh->get_set_names(EntityRank::FACE);

            //loop through all the set names and finds sets attached to input side numbers
            std::string tPairFirstName = "SideSet_";
            tPairFirstName.append(mMeshSideSetPairs(aPairCount)(0));

            //loop through all the set names and finds sets attached to input side numbers
            std::string tPairSecondName = "SideSet_";
            tPairSecondName.append(mMeshSideSetPairs(aPairCount)(1));

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

            //populate the first side clusters
            for (uint i = 0; i < tFirstSideSetNames.size() ; i++)
            {
                moris::mtk::Set* tSet = tIntegrationMesh->get_set_by_name( tFirstSideSetNames(i) );

                // get clusters in the second set
                moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet->get_clusters_on_set();

                tFirstSideClusters.append(tSetClusters);

                if( tFirstSideClusters.size() > 0 )
                {
                    break;
                }
            }

            //populate the second side clusters
            for (uint i = 0; i < tSecondSideSetNames.size() ; i++)
            {
                moris::mtk::Set* tSet = tIntegrationMesh->get_set_by_name( tSecondSideSetNames(i) );

                // get clusters in the second set
                moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet->get_clusters_on_set();

                tSecondSideClusters.append(tSetClusters);

                if( tSecondSideClusters.size() > 0 )
                {
                    break;
                }
            }

            moris::Cell<moris::mtk::Cell const *> const & tSecondCells = tSecondSideClusters( 0 )->get_primary_cells_in_cluster();

            // get the side ordinals of the integration cells in the cluster
            moris::Matrix<moris::IndexMat> tSecondCellOrds = tSecondSideClusters( 0 )->get_cell_side_ordinals();

            //vertices on the second one
            moris::Cell<moris::mtk::Vertex const *> tSecondVertices = tSecondCells( 0 )->get_vertices_on_side_ordinal( tSecondCellOrds( 0 ) );

            moris::Matrix< moris::DDRMat > tSecondVertexCoords = tSecondVertices(0)->get_coords();

            moris::Matrix<moris::DDRMat> tSecondSideNormal = tSecondCells(0)->compute_outward_side_normal(tSecondCellOrds(0));

            //get the Integration cell cluster for the side cluster
            moris::Cell<moris::mtk::Cell const *> const & tFirstCells = tFirstSideClusters( 0 )->get_primary_cells_in_cluster();

            // get the side ordinals of the integration cells in the cluster
            moris::Matrix<moris::IndexMat> tFirstCellOrds = tFirstSideClusters( 0 )->get_cell_side_ordinals();

            //get the vertex information
            moris::Cell<moris::mtk::Vertex const *> tFirstVertices = tFirstCells( 0 )->get_vertices_on_side_ordinal( tFirstCellOrds( 0 ) );

            moris::Matrix< moris::DDRMat > tFirstVertexCoords = tFirstVertices(0)->get_coords();

            //form the output
            tOffsetVector = dot(tSecondVertexCoords-tFirstVertexCoords, tSecondSideNormal )*trans(tSecondSideNormal);
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect::create_master_vertices(moris::Cell< moris::mtk::Vertex *> & aMasterVertices, Matrix < DDRMat> & tPhysicalCoordinates3D,moris::mtk::Cell const & aMasterInterpCell, uint aPairCount, Matrix < DDRMat> & tNewNodeCoordinates)
        {
            uint tNumSurfaceNodes = tPhysicalCoordinates3D.n_rows();

            //initialize the vertices of the cell ( +1 for the head of the tet4 )
            aMasterVertices.resize( tNumSurfaceNodes+1,nullptr ) ;

            // allocate ids for all the vertices
            moris_id tVertexId = mIntersectedMeshData.allocate_entity_ids_external_entity_data( tNumSurfaceNodes+1,  EntityRank::NODE);

            for( uint i = 0 ; i < tNumSurfaceNodes ; i++)
            {
                //update the vertex index
                moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;

                //construct the new vertex
                aMasterVertices ( i ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tPhysicalCoordinates3D.get_row(i) ) ;

                // Update the index and id
                mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );
                tVertexId++;
            }

            moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;
            mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );

            //add the top vertex of tet4
            aMasterVertices( tNumSurfaceNodes ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewNodeCoordinates) ;

        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect::create_slave_vertices(moris::Cell< moris::mtk::Vertex *> & aSlaveVertices, Matrix < DDRMat> & tPhysicalCoordinates3D,moris::mtk::Cell const & aSlaveInterpCell, uint aPairCount,Matrix < DDRMat> & tNewNodeCoordinates)
        {
            uint tNumSurfaceNodes = tPhysicalCoordinates3D.n_rows();

            //initialize the vertices of the cell ( +1 for the head of the tet4 )
            aSlaveVertices.resize( tNumSurfaceNodes+1,nullptr ) ;

            // allocate ids for all the vertices
            moris_id tVertexId = mIntersectedMeshData.allocate_entity_ids_external_entity_data( tNumSurfaceNodes+1,  EntityRank::NODE);

            for( uint i = 0 ; i < tNumSurfaceNodes ; i++)
            {
                //update the vertex index
                moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;

                //construct the new vertex
                aSlaveVertices ( i ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex,tPhysicalCoordinates3D.get_row(i)) ;

                mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );
                tVertexId++;
            }

            moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;
            mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );

            //add the top vertex of tet4
            aSlaveVertices( tNumSurfaceNodes ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewNodeCoordinates) ;
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect::rotation_matrix(moris::Matrix< DDRMat > &aRotation, moris::Matrix< DDRMat > &aInverseRotation, uint aPairCount ) const
        {
            if( mMeshSideSetPairs(aPairCount)(0) == "4")
            {
                aRotation        = { { 1, 0, 0 }, { 0, 0, 1 }, { 0, 1, 0 } };
                aInverseRotation = { { 1, 0, 0 }, { 0, 0, 1 }, { 0, 1, 0 } };
            }
            else if ( mMeshSideSetPairs(aPairCount)(0) == "1" )
            {
                aRotation        = { { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, 1 } };
                aInverseRotation = { { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, 1 } };
            }
            else
            {
                aRotation        = { { 0, 0, 1 }, { 1, 0, 0 }, { 0, 1, 0 } };
                aInverseRotation = { { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 0 } };
            }
        }

        // ----------------------------------------------------------------------------

        std::pair< moris_id, moris_id >
        Intersection_Detect::permutation_pair(uint const & aPairCount) const
        {
            std::pair< moris_id, moris_id > tPermuationOrder;

            if( mMeshSideSetPairs(aPairCount)(0) == "4")
            {
                tPermuationOrder = std::make_pair ( 2 , 1 ) ;
            }
            else if ( mMeshSideSetPairs(aPairCount)(0) == "1" )
            {
                tPermuationOrder = std::make_pair( 0 , 2 ) ;
            }
            else
            {
                tPermuationOrder =  std::make_pair( 0 , 1 ) ;
            }
            return tPermuationOrder;
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect::build_input_data_for_intersection(
                mtk::Cluster const *                  aSideCluster,
                moris::Matrix< DDRMat > &             aInverseRotation,
                moris::Matrix< DDRMat >  &            aUniqueParamCoords,
                moris::Matrix< moris::IdMat > &       aLocalElemToNode)
        {
            // Integration cells in the cluster
            moris::Cell< moris::mtk::Cell const *>  const & tIntegCells = aSideCluster->get_primary_cells_in_cluster();

            // get the side ordinals of the integration cells in the cluster
            moris::Matrix<moris::IndexMat> tCellOrdinals = aSideCluster->get_cell_side_ordinals();

            // Element to node map
            moris::Matrix<moris::IdMat> tElementToNode( tIntegCells.size(), 3 );

            //parametric coordinate matrix for side 2, all the coordinates
            moris::Matrix <DDRMat> tParamCoordsMatrix( 2, 3 * tIntegCells.size() );

            //loop over the cells on side 2
            for( uint k = 0 ; k < tIntegCells.size() ; k++ )
            {
                //extract vertices laying on the side
                moris::Cell<moris::mtk::Vertex const *> tVertices = tIntegCells( k )->get_vertices_on_side_ordinal( tCellOrdinals( k ) );

                //loop over the vertices
                for ( uint j = 0 ; j < tVertices.size() ; j++ )
                {
                    //element to node map globally
                    tElementToNode( k, j ) = tVertices( j )->get_id();

                    //parametric coordinates of the vertices
                    Matrix<DDRMat> tTmpCoordMatrix = aSideCluster->get_vertex_local_coordinate_wrt_interp_cell( tVertices( j ) ).get_row( 0 );

                    //convert to desired form
                    tTmpCoordMatrix = tTmpCoordMatrix * aInverseRotation ;

                    //convert to 2d
                    tParamCoordsMatrix( 0, k * 3 + j ) = tTmpCoordMatrix( 1 );
                    tParamCoordsMatrix( 1, k * 3 + j ) = tTmpCoordMatrix( 2 );

                }
            }

            //find the unique ids to form the local to global map ( in order to work with the intersection algo)
            moris::Matrix<moris::IdMat> tUniqueElementIds;
            unique( tElementToNode, tUniqueElementIds );

            //construct a map of local indices to unique global id's
            std::unordered_map< moris::moris_id, moris::moris_index > tGlobaltoLocalmap;

            for(moris::uint iB =0; iB < tUniqueElementIds.numel() ; iB++)
            {

                tGlobaltoLocalmap[tUniqueElementIds(iB)] = iB+1;
            }

            //construct the input of unique parametric coordinates
            aUniqueParamCoords.resize(2,tUniqueElementIds.numel() );

            //loop over the unique ids
            for (uint i = 0 ; i < tUniqueElementIds.numel() ; i++ )
            {
                //populate the unique parametric coordinate matrix
                for( uint j = 0 ; j < tIntegCells.size() ; j++ )
                {
                    bool tFound = false;
                    //loop over the cells
                    for ( uint k = 0 ; k < 3 ; k++ )
                    {
                        //if the coordinate was found in global element to node map, we have the coordinate stored
                        if(  tUniqueElementIds(i) == tElementToNode(j,k) )
                        {
                            aUniqueParamCoords.get_column(i) = tParamCoordsMatrix.get_column(j*3+k);

                            tFound = true ;
                            break;
                        }
                    }
                    //no further loop if a coordinate is found
                    if( tFound )
                    {
                        break;
                    }
                }
            }

            // assemble element to node map based on the local indices (input for intersection algo)
            aLocalElemToNode.resize( tElementToNode.n_rows(), tElementToNode.n_cols() );

            for( uint tRowNum = 0 ; tRowNum < tElementToNode.n_rows() ; tRowNum++ )
            {
                for( uint tColNum = 0 ; tColNum < tElementToNode.n_cols() ; tColNum++ )
                {

                    aLocalElemToNode( tRowNum, tColNum ) = tGlobaltoLocalmap.at( tElementToNode( tRowNum, tColNum ) );
                }
            }
        }

        // ----------------------------------------------------------------------------

        bool
        Intersection_Detect::clusters_align(
                mtk::Cluster const *                          aFirstSideCluster,
                mtk::Cluster const *                          aSecondSideCluster,
                uint                                          aPairCount,
                moris::Matrix< DDRMat > const &               aOffsetMatrix) const
        {
            //right interpolation cell
            moris::mtk::Cell const & tInterpCell1 = aFirstSideCluster->get_interpolation_cell();

            //left interpolation cell
            moris::mtk::Cell const & tInterpCell2 = aSecondSideCluster->get_interpolation_cell();

            //coordinates to compare
            Matrix<DDRMat>  const & tInterpCoords1 = tInterpCell1.get_cell_geometric_coords_on_side_ordinal( std::stoi( mMeshSideSetPairs( aPairCount )( 0 ) ) - 1 );

            Matrix<DDRMat>  const & tInterpCoords2 = tInterpCell2.get_cell_geometric_coords_on_side_ordinal( std::stoi( mMeshSideSetPairs( aPairCount )( 1 ) ) - 1 );

            //subtract the offset to compare the coords
            Matrix<DDRMat>   const & tInterpCoords2Shifted = tInterpCoords2 - aOffsetMatrix;

            bool tMatched = false;

            // matching interpolation cell criteria (This condition might be necessary but not enough)
            if( std::abs( norm( tInterpCoords2Shifted ) - norm( tInterpCoords1 ) ) < 0.0000000001 )
            {
                //count nodes that are matching
                uint s = 0 ;

                //loop over first nodes
                for ( uint i = 0 ; i < 4 ; i++ )
                {
                    //loop over second nodes
                    for( uint j = 0 ; j < 4 ; j++)
                    {
                        if ( norm(tInterpCoords2Shifted.get_row(i)  - tInterpCoords1.get_row(j) ) < 0.00000000001 )
                        {
                            s++;
                            break;
                        }
                    }
                }
                if( s == 4 )
                {
                    tMatched = true;
                }
            }

            return tMatched;
        }

        // ----------------------------------------------------------------------------

        moris::Cell< moris::Matrix < moris::DDRMat > >
        Intersection_Detect::elementwise_intersection_search(
                moris::Matrix < moris::DDRMat >   const &  aFirstTRICoords,
                moris::Matrix < moris::DDRMat >   const &  aSecondTRICoords,
                moris::Matrix < moris::IdMat >    const &  aFirstTRINodeIndex,
                moris::Matrix < moris::IdMat >    const &  aSecondTRINodeIndex,
                moris::Matrix < moris::IndexMat > const &  aFirstTRIConnect,
                moris::Matrix < moris::IndexMat > const &  aSecondTRIConnect)
        {

            //###########################################################
            // The same notation as the original MATLAB algorithm is used
            //###########################################################

            moris::Cell< moris::Matrix < moris::DDRMat > > allIntersectionTriangles ;

            //list of triangles of second side to treat
            moris::Cell< uint > bl;

            //list of triangles first side to start with
            moris::Cell< uint > bil;

            //find the first matching triangles
            bool tExit = false;
            for( uint i = 0 ; i < aSecondTRINodeIndex.n_rows(); i++ )
            {
                for( uint j = 0 ; j < aFirstTRINodeIndex.n_rows(); j++)
                {
                    Matrix < IdMat > Tbbc = aSecondTRINodeIndex( {i,i}, {0, 2});
                    Matrix < IdMat > Taac = aFirstTRINodeIndex( {j,j}, {0, 2});

                    Matrix < DDRMat > NbTbbc( 2, 3 );
                    Matrix < DDRMat > NaTaac( 2, 3 );

                    for ( uint Ii = 0 ; Ii < 3 ; Ii++ )
                    {
                        NbTbbc.get_column( Ii ) = aSecondTRICoords.get_column( Tbbc( Ii ) - 1 );
                        NaTaac.get_column( Ii ) = aFirstTRICoords.get_column( Taac( Ii ) - 1 );
                    }

                    moris::Matrix < moris::DDRMat  > tP;

                    moris::Matrix < moris::DDUMat  >  tnc;

                    this->Intersect( NbTbbc, NaTaac, tP, tnc ) ;

                    if (tP.n_cols()>2)
                    {
                        bl.push_back( i+1 );

                        bil.push_back(j+1);

                        tExit = true;

                        break;
                    }
                }
                if( tExit )
                {
                    break;
                }
            }

            /* end of this part */

            if( tExit == false)
            {
                bl.push_back( 1 );

                bil.push_back(1);
            }

            //flag for triangles in second side treated
            moris::Cell< uint > bd( aSecondTRIConnect.n_rows()+1 , 0);

            //mark first triangle in b list.
            bd(0) = 1 ;

            //guard, to treat boundaries
            bd( aSecondTRIConnect.n_rows() ) = 1;

            uint bc;

            moris::Cell< uint > n( 3, 0 );

            while ( bl.size() > 0 )
            {
                bc = bl(0) - 1;
                // current triangle of second side
                bl.erase(0);
                bl.shrink_to_fit();

                //triangle of the first side to start with
                moris::Cell< uint > al(1 , bil(0));
                bil.erase( 0 );
                bil.shrink_to_fit();

                // same as for bd
                moris::Cell< uint > ad( aFirstTRIConnect.n_rows()+1 , 0.0);
                for (uint Ii = 0 ; Ii < al.size(); Ii++ )
                {
                    ad ( al (Ii) - 1) = 1.0 ;
                }

                ad( aFirstTRIConnect.n_rows() ) = 1;

                //triangles intersecting with neighbors
                n.resize( 3, 0 ) ;

                while ( al.size() > 0 )
                {
                    //take next candidate
                    uint ac = al(0) - 1;
                    al.erase( 0 );
                    al.shrink_to_fit();

                    Matrix < IdMat > Tbbc = aSecondTRINodeIndex( {bc,bc}, {0, 2});
                    Matrix < IdMat > Taac = aFirstTRINodeIndex( {ac,ac}, {0, 2});

                    Matrix < DDRMat > NbTbbc( 2, 3 );
                    Matrix < DDRMat > NaTaac( 2, 3 );

                    for ( uint Ii = 0 ; Ii < 3 ; Ii++ )
                    {
                        NbTbbc.get_column( Ii ) = aSecondTRICoords.get_column( Tbbc( Ii ) - 1 );
                        NaTaac.get_column( Ii ) = aFirstTRICoords.get_column( Taac( Ii ) - 1 );
                    }

                    moris::Matrix < moris::DDRMat  > tP;

                    moris::Matrix < moris::DDUMat  >  tnc;

                    this->Intersect( NbTbbc, NaTaac, tP, tnc ) ;


                    if ( tP.n_cols() != 0 )
                    {
                        // if intersection is more dimension than a line
                        if ( tP.n_cols() > 2 )
                        {
                            allIntersectionTriangles.push_back(tP) ;
                        }

                        Matrix < IdMat> Taac4to6 =  aFirstTRIConnect( {ac,ac}, {0, 2});
                        Matrix < DDUMat> tmp( 1, Taac4to6.n_cols() );


                        for (uint Ii = 0 ; Ii < Taac4to6.n_cols() ; Ii++ )
                        {
                            tmp(Ii) = ad( Taac4to6(Ii) - 1);
                        }


                        Matrix < DDUMat> comp( 1, tmp.n_cols() ,0.0 );

                        Matrix < DDNIMat> tn = find( tmp == comp ) ;


                        if ( tn.numel() != 0)
                        {

                            Matrix < IdMat> t(1,tn.numel()) ;

                            for ( uint Ii = 0 ; Ii < t.n_cols() ; Ii++ )
                            {
                                t ( Ii ) = aFirstTRIConnect(ac, tn(Ii) );
                            }

                            //add neighbours
                            for ( uint Ii = 0; Ii < t.n_cols() ; Ii++ )
                            {
                                al.push_back( t( Ii ) );
                            }


                            for ( uint Ii = 0; Ii < t.n_cols() ; Ii++ )
                            {
                                ad( t( Ii ) - 1) = 1.0 ;
                            }
                        }

                        // n(find(nc>0))  = ac;
                        Matrix < DDUMat > comp2 ( 1, tnc.n_cols(), 0 );

                        Matrix < DDNIMat > findnc = find(  comp2 < tnc) ;

                        for (uint Ii = 0 ; Ii < findnc.numel() ; Ii++ )
                        {
                            n( findnc( Ii )  ) = ac + 1 ;
                        }

                    }
                }

                // find non-treated neighbors
                // tmp=find(bd(Tb(bc,4:6))==0)
                Matrix < IdMat> Tbbc4to6 =  aSecondTRIConnect( {bc,bc}, {0, 2});

                Matrix < DDUMat> bdTbbc( 1, Tbbc4to6.n_cols() ) ;

                for (uint Ii = 0 ; Ii < Tbbc4to6.n_cols() ; Ii++ )
                {
                    bdTbbc( Ii ) = bd( Tbbc4to6( Ii ) - 1);
                }

                Matrix < DDUMat> zerovec( 1, bdTbbc.n_cols() ,0.0 );

                Matrix < DDNIMat > tmp = find( bdTbbc == zerovec );


                if (tmp.numel() != 0 )
                {
                    // take those which intersect
                    // idx=find(n(tmp)>0)
                    Matrix < DDUMat> ntmp( 1, tmp.n_rows() );

                    for (uint Ii = 0 ; Ii < ntmp.n_cols() ; Ii++ )
                    {
                        ntmp( Ii ) = n( tmp ( Ii ) );
                    }

                    Matrix < DDUMat> zerovec2( 1, ntmp.n_cols() ,0.0 );

                    Matrix < DDNIMat > idx = find ( zerovec2 < ntmp );

                    if ( idx.numel() != 0 )
                    {

                        // t=Tb(bc,3+tmp(idx))
                        Matrix < DDUMat> tmpidx(1, idx.n_rows() ) ;

                        for( uint Ii = 0 ; Ii < tmpidx.n_cols() ; Ii++ )
                        {
                            tmpidx (Ii) = tmp ( idx ( Ii ) ) ;
                        }

                        Matrix < IdMat > t(1,tmpidx.n_cols() ) ;

                        for ( uint Ii = 0 ; Ii < t.n_cols() ; Ii++ )
                        {
                            t( Ii ) =  aSecondTRIConnect (bc, tmpidx (Ii)) ;
                        }

                        //and add the ones that intersect
                        for ( uint Ii = 0 ; Ii < t.n_cols() ; Ii++ )
                        {
                            bl.push_back( t( Ii ) ) ;
                        }

                        //with starting candidates Ta
                        // bil=[bil n(tmp(idx))]
                        if ( tmpidx.numel() != 0 )
                        {
                            Matrix < DDUMat> ntmpidx(1, idx.numel() ) ;

                            for(uint Ii = 0 ; Ii < ntmpidx.n_cols() ; Ii++ )
                            {
                                ntmpidx( Ii ) = n ( tmpidx (Ii) );
                            }

                            for(uint Ii = 0 ; Ii < ntmpidx.n_cols() ; Ii++ )
                            {
                                bil.push_back( ntmpidx( Ii ) );
                            }

                            for ( uint Ii = 0 ; Ii < t.n_cols() ; Ii++ )
                            {
                                bd( t( Ii ) - 1 ) = 1 ;
                            }

                        }
                    }
                }
            }

            //return the cell of intersected
            return allIntersectionTriangles;
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect::generate_identifier( moris::Cell< mtk::Cluster const * > & aSideClusters,
                uint &                                 aPairCount,
                Matrix<IndexMat>                    &   aIdentifierMatrix,
                std::unordered_map< moris::moris_index, moris::Cell<moris_index> > & aBackgroundCellToSideClusterMap) const
        {
            //reserve memory for output matrix
            uint tOutputSize = aSideClusters.size();
            aIdentifierMatrix.set_size(1, tOutputSize);

            moris::mtk::Integration_Mesh*  tIntegrationMesh  = mMeshManager->get_integration_mesh( mMeshIndex );

            //multiplier for the identification number
            uint tMultiplier = tIntegrationMesh->get_num_entities(moris::EntityRank::NODE);

            uint tFirst;
            uint tSecond ;
            // depending on the side pick relevant ijk
            if( mMeshSideSetPairs(aPairCount)(0) == "4")
            {
                //pick j and k
                tFirst = 1;
                tSecond = 2;
            }
            else if ( mMeshSideSetPairs(aPairCount)(0) == "1" )
            {
                //pick i and k
                tFirst = 0;
                tSecond = 2;
            }
            else
            {
                //pick i and j
                tFirst = 0;
                tSecond = 1;
            }

            //loop over the clusters
            for(uint iI = 0 ; iI < tOutputSize; iI++ )
            {
                //access the ijk postion
                const luint * tIJK = aSideClusters(iI)->get_interpolation_cell().get_ijk();

                //calculate the identifier
                moris_index tIdentifier = tIJK[tFirst] * tMultiplier + tIJK[tSecond] ;

                //assign it to the matrix
                aIdentifierMatrix(iI) = tIdentifier;

                //add the cluster index to the relevant bg cell
                aBackgroundCellToSideClusterMap[tIdentifier].push_back(iI);
            }

        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect::group_cut_cells(moris::Matrix< IndexMat> const & aCutTrianglesIdentifier,
                std::unordered_map< moris_index , moris::Cell< moris_index > > & aCutCellIdentifierToCutCellIndex ) const
        {
            //loop over all the cut cells
            for(uint iI = 0 ; iI < aCutTrianglesIdentifier.n_rows() ; iI++ )
            {
                aCutCellIdentifierToCutCellIndex[aCutTrianglesIdentifier(iI)].push_back(iI);
            }
        }
    }

}
