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
                //find the rotation and inverse rotation
                //rotation goes from 2d to 3d
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

                //tFirstSideSetNames = { "SideSet_2_c_p1"};
                //tSecondSideSetNames = { "SideSet_4_c_p1"};
                //find the offset matrix
                moris::Matrix < DDRMat >  tVecTransform( 4, 1, 1 );
                moris::Matrix< moris::DDRMat > const & tOffsetMatrix = tVecTransform * tOffsetVector ;

                //initialize all the clusters on each side
                moris::Cell< mtk::Cluster const * > tSideClusters1;
                moris::Cell< mtk::Cluster const * > tSideClusters2;

                //cells on the side set are locally indexed and they are assigned a phase
                moris::Cell< moris_index > tLocalCellNumToColor1 ;

                //global cell id to local cell index map
                std::unordered_map< moris::moris_id,moris::moris_index > tGlobaltoLocalmap1Phase;

                //populate the first side clusters and assign them phases
                uint tLocalCellIndex = 0 ;

                //loop over the
                for (uint i = 0; i < tFirstSideSetNames.size() ; i++)
                {
                    //get the set on the first side
                    moris::mtk::Set* tSet = tIntegrationMesh->get_set_by_name( tFirstSideSetNames(i) );

                    //get color(phase) associated with the set
                    Matrix<IndexMat> tColor = tSet->get_set_colors();

                    // get clusters in the second set
                    moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet->get_clusters_on_set();

                    //loop over clusters on the set and assign the first cell in the cluster a color(phase)
                    for( uint j = 0 ; j < tSetClusters.size(); j++ )
                    {
                        //primary cells in cluster
                        moris::Cell< moris::mtk::Cell const *> const & tPrimaryCellsinCluster = tSetClusters( j )->get_primary_cells_in_cluster();

                        //build the map of global id to local index for side cells
                        tGlobaltoLocalmap1Phase[ tPrimaryCellsinCluster( 0 )->get_id() ] = tLocalCellIndex;
                        tLocalCellIndex++;

                        //store the phase of cell
                        tLocalCellNumToColor1.push_back( tColor( 0 ) );
                    }

                    //append the clusters the list of the first side clusters
                    tSideClusters1.append( tSetClusters );
                }

                //cells on the side set are locally indexed and they are assigned a phase to each one
                moris::Cell< moris_index > tLocalCellNumToColor2 ;

                //global cell id to local cell index map
                std::unordered_map<moris::moris_id,moris::moris_index> tGlobaltoLocalmap2Phase;

                //populate the first side clusters
                tLocalCellIndex = 0 ;

                //loop over the existing side sets on the 2nd side
                for (uint i = 0; i < tSecondSideSetNames.size() ; i++)
                {
                    //get the set on the first side
                    moris::mtk::Set* tSet = tIntegrationMesh->get_set_by_name( tSecondSideSetNames( i ) );

                    // get clusters in the second set
                    moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet->get_clusters_on_set();

                    //get color(phase) associated with the set
                    Matrix<IndexMat> tColor = tSet->get_set_colors();

                    for( uint j = 0 ; j < tSetClusters.size() ; j++ )
                    {
                        //primary cells in cluster
                        moris::Cell< moris::mtk::Cell const *> const & tPrimaryCellsinCluster = tSetClusters( j )->get_primary_cells_in_cluster();

                        //build the map of global id to local index for side cells
                        tGlobaltoLocalmap2Phase[ tPrimaryCellsinCluster( 0 )->get_id() ] = tLocalCellIndex;
                        tLocalCellIndex++;

                        //store the phase of cell
                        tLocalCellNumToColor2.push_back( tColor( 0 ) );
                    }

                    //append the clusters the list of the first side clusters
                    tSideClusters2.append( tSetClusters );
                }

                //loop over all the clusters on the first side
                for( uint tClusterNum1  = 0 ; tClusterNum1 < tSideClusters1.size(); tClusterNum1++ )
                {
                    //loop over all the clusters on the second side
                    for( uint tClusterNum2 = 0 ; tClusterNum2 < tSideClusters2.size() ; tClusterNum2++ )
                    {
                        // Two interpolation cells matched
                        if( this->clusters_align( tSideClusters1( tClusterNum1 ), tSideClusters2( tClusterNum2 ), tPairCount, tOffsetMatrix  ) )
                        {
                            // right interpolation cell
                            moris::mtk::Cell const & tInterpCell1 = tSideClusters1( tClusterNum1 )->get_interpolation_cell();

                            //left interpolation cell
                            moris::mtk::Cell const & tInterpCell2 = tSideClusters2( tClusterNum2 )->get_interpolation_cell();

                            //the case where left and right clusters are trivial
                            if ( tSideClusters1( tClusterNum1 )->is_trivial() and tSideClusters2( tClusterNum2 )->is_trivial() )
                            {
                                moris::Cell< Matrix<DDRMat> >tParamCoordsCell2 (4);
                                tParamCoordsCell2(0) = { { -1.0, +1.0,  0.0 }, { -1.0, -1.0,  0.0 } };
                                tParamCoordsCell2(1) = { { -1.0,  0.0, -1.0 }, { -1.0,  0.0, +1.0 } };
                                tParamCoordsCell2(2) = { { -1.0,  0.0,  1.0 }, { +1.0,  0.0,  1.0 } };
                                tParamCoordsCell2(3) = { {  1.0,  0.0, +1.0 }, {  1.0,  0.0, -1.0 } };

                                //get the Integration cell cluster for the side cluster
                                moris::Cell<moris::mtk::Cell const *> const & tCells2 = tSideClusters2( tClusterNum2 )->get_primary_cells_in_cluster();

                                //get the Integration cell cluster for the side cluster
                                moris::Cell<moris::mtk::Cell const *> const & tCells1 = tSideClusters1( tClusterNum1 )->get_primary_cells_in_cluster();

                                //obtain the phases of each side
                                moris_index tPhase2  = tLocalCellNumToColor2 ( tGlobaltoLocalmap2Phase[ tCells2( 0 )->get_id() ] );
                                moris_index tPhase1  = tLocalCellNumToColor1 ( tGlobaltoLocalmap1Phase[ tCells1( 0 )->get_id() ] );

                                //phase interaction table value
                                moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                                this->make_new_pairs( tParamCoordsCell2, tInterpCell1, tInterpCell2, tPairCount, tPhaseToPhaseIndex);

                            }

                            //case where left side is trivial
                            else if ( tSideClusters1( tClusterNum1 )->is_trivial() == true and tSideClusters2( tClusterNum2 )->is_trivial() == false )
                            {
                                //get the Integration cell cluster for the side cluster
                                moris::Cell<moris::mtk::Cell const *> const & tCells2 = tSideClusters2( tClusterNum2 )->get_primary_cells_in_cluster();

                                // get the side ordinals of the integration cells in the cluster
                                moris::Matrix< moris::IndexMat > tCellOrds2 = tSideClusters2( tClusterNum2 )->get_cell_side_ordinals();

                                //cell to store parametric coordinates of all triangle
                                moris::Cell < moris::Matrix <DDRMat> > tParamCoordsCell2;

                                for( uint i = 0 ; i < tCells2.size() ; i++ )
                                {
                                    //find vertices on the side
                                    moris::Cell<moris::mtk::Vertex const *> tVertices = tCells2( i )->get_vertices_on_side_ordinal( tCellOrds2(i) );

                                    //parametric coordinates of each triangle
                                    moris::Matrix <DDRMat> tParamCoordsMatrix2(2,3);

                                    //loop over the vertices on the side
                                    for ( uint j = 0 ; j < tVertices.size() ; j++ )
                                    {
                                        //temporary matrix of 3d parametric coordinates
                                        Matrix<DDRMat> tTmpCoordMatrix = tSideClusters2( tClusterNum2 )->get_vertex_local_coordinate_wrt_interp_cell( tVertices( j ) ).get_row( 0 );

                                        //organize 3d parametric coordinates in the desired format
                                        tTmpCoordMatrix = tTmpCoordMatrix * tInverseRotation ;

                                        //convert 3d parametric coordinates to 2d
                                        tParamCoordsMatrix2( 0, j ) = tTmpCoordMatrix(1);
                                        tParamCoordsMatrix2( 1, j ) = tTmpCoordMatrix(2);
                                    }

                                    //sort them counter clockwise
                                    this->SortandRemove(tParamCoordsMatrix2 );

                                    //store each triangle in the the cell
                                    tParamCoordsCell2.push_back( tParamCoordsMatrix2 );
                                }

                                //get the Integration cell cluster for the side cluster
                                moris::Cell<moris::mtk::Cell const *> const & tCells1 = tSideClusters1( tClusterNum1 )->get_primary_cells_in_cluster();

                                //obtain the phases of each side
                                moris_index tPhase2  = tLocalCellNumToColor2 ( tGlobaltoLocalmap2Phase[ tCells2( 0 )->get_id() ] );
                                moris_index tPhase1  = tLocalCellNumToColor1 ( tGlobaltoLocalmap1Phase[ tCells1( 0 )->get_id() ] );

                                //phase interaction table value
                                moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                                //creates the pair of the single sided clusters
                                this->make_new_pairs( tParamCoordsCell2, tInterpCell1, tInterpCell2, tPairCount, tPhaseToPhaseIndex );
                            }

                            //case where left side is trivial
                            else if ( tSideClusters1( tClusterNum1 )->is_trivial() == false and tSideClusters2( tClusterNum2 )->is_trivial() == true )
                            {
                                //get the Integration cell cluster for the side cluster
                                moris::Cell<moris::mtk::Cell const *> const & tCells1 = tSideClusters1( tClusterNum1 )->get_primary_cells_in_cluster();

                                // get the side ordinals of the integration cells in the cluster
                                moris::Matrix< moris::IndexMat > tCellOrds1 = tSideClusters1( tClusterNum1 )->get_cell_side_ordinals();

                                //cell to store parametric coordinates of all triangle
                                moris::Cell< moris::Matrix <DDRMat> > tParamCoordsCell1;

                                for( uint i = 0 ; i < tCells1.size() ; i++ )
                                {
                                    //find vertices on the 1st side cluster
                                    moris::Cell<moris::mtk::Vertex const *> tVertices = tCells1(i)->get_vertices_on_side_ordinal( tCellOrds1(i) );

                                    //parametric coordinates of each triangle
                                    moris::Matrix <DDRMat> tParamCoordsMatrix1( 2, 3 );

                                    //loop over the vertices on the side
                                    for ( uint j = 0 ; j < tVertices.size() ; j++ )
                                    {
                                        //temporary matrix of 3d parametric coordinates
                                        Matrix<DDRMat> tTmpCoordMatrix = tSideClusters1( tClusterNum1 )->get_vertex_local_coordinate_wrt_interp_cell( tVertices( j ) ).get_row( 0 );

                                        //convert to desired form
                                        tTmpCoordMatrix  = tTmpCoordMatrix  * tInverseRotation ;

                                        //convert 3d parametric coordinates to 2d
                                        tParamCoordsMatrix1( 0, j ) = tTmpCoordMatrix (1);
                                        tParamCoordsMatrix1( 1, j ) = tTmpCoordMatrix (2);

                                    }
                                    //sort them counter clockwise
                                    this->SortandRemove(tParamCoordsMatrix1 );

                                    //store each triangle in the the cell
                                    tParamCoordsCell1.push_back(tParamCoordsMatrix1);
                                }

                                //get the Integration cell cluster for the side cluster
                                moris::Cell<moris::mtk::Cell const *> const & tCells2 = tSideClusters2( tClusterNum2 )->get_primary_cells_in_cluster();

                                //obtain the phases of each side
                                moris_index tPhase2  = tLocalCellNumToColor2 ( tGlobaltoLocalmap2Phase[ tCells2( 0 )->get_id() ] );
                                moris_index tPhase1  = tLocalCellNumToColor1 ( tGlobaltoLocalmap1Phase[ tCells1( 0 )->get_id() ] );

                                //phase interaction table value
                                moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                                this->make_new_pairs( tParamCoordsCell1, tInterpCell1, tInterpCell2, tPairCount, tPhaseToPhaseIndex );
                            }

                            //where both cases are not trivial
                            else
                            {
                                /* the following part supplied the input for the intersection algorithim
                                 * There are six inputs needs to be constructed
                                 * 1)unique coordinates of the points in each side
                                 * 2)element to node map for each side
                                 * 3)connectivity info for each face
                                 */
                                // Initialize the data for side 1
                                moris::Matrix<DDRMat>                       aUniqueParamCoords2;
                                moris::Matrix<moris::IdMat>                  aLocalElemToNode2;

                                // Build the input data
                                this->build_input_data_for_intersection(
                                        tSideClusters2( tClusterNum2 ),
                                        tInverseRotation,
                                        aUniqueParamCoords2,
                                        aLocalElemToNode2);

                                // Initialize the data for side 1
                                moris::Matrix<DDRMat>                       aUniqueParamCoords1;
                                moris::Matrix<moris::IdMat>                  aLocalElemToNode1;

                                // Build the input data
                                this->build_input_data_for_intersection(
                                        tSideClusters1( tClusterNum1 ),
                                        tInverseRotation,
                                        aUniqueParamCoords1,
                                        aLocalElemToNode1);

                                //Intersection Algorithm that returns a cell of intersected area
                                moris::Cell< Matrix < DDRMat > > allIntersectionTriangles = this->elementwise_bruteforce_search(
                                        aUniqueParamCoords1, aUniqueParamCoords2,
                                        aLocalElemToNode1, aLocalElemToNode2);

                                //if there is any intersection
                                if( allIntersectionTriangles.size() > 0 )
                                {
                                    //get the Integration cell cluster for the side cluster
                                    moris::Cell<moris::mtk::Cell const *> const & tCells2 = tSideClusters2( tClusterNum2 )->get_primary_cells_in_cluster();

                                    //get the Integration cell cluster for the side cluster
                                    moris::Cell<moris::mtk::Cell const *> const & tCells1 = tSideClusters1( tClusterNum1 )->get_primary_cells_in_cluster();

                                    //obtain the phases of each side
                                    moris_index tPhase2  = tLocalCellNumToColor2 ( tGlobaltoLocalmap2Phase[ tCells2(0)->get_id() ] );
                                    moris_index tPhase1  = tLocalCellNumToColor1 ( tGlobaltoLocalmap1Phase[ tCells1(0)->get_id() ] );

                                    //phase interaction table value
                                    moris_index tPhaseToPhaseIndex = tPhaseInteractionTable( tPhase1, tPhase2 );

                                    this->make_new_pairs( allIntersectionTriangles, tInterpCell1, tInterpCell2, tPairCount,  tPhaseToPhaseIndex);
                                }
                            }
                        } /* end if when they are matched */

                    }/* iteration over the second (slave) side */

                }/* iteration over the first (master) side */
            }
            //add the double sided set based on the interactions
            this->constrcuct_add_dbl_sided_set(tPhaseInteractionTable);
        }


        //---------------------------------------------------------------------

        void
        Intersection_Detect::EdgeIntersect(
                moris::Matrix < DDRMat > const & aFirstTRICoords,
                moris::Matrix < DDRMat > const &  aSecondTRICoords,
                moris::Matrix < DDRMat > & aIntersectedPoints,
                moris::Matrix < DDUMat > & aIntersectVec)
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
                moris::Matrix < moris::DDRMat  >       &  aIntersectedPoints)
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
        ( moris::Matrix < moris::DDRMat  > & aIntersectedPoints )
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
                moris::Matrix < moris::DDUMat  > & aIntersectVec)
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

        moris::Cell< moris::Matrix < moris::DDRMat > >
        Intersection_Detect::elementwise_bruteforce_search (
                moris::Matrix < moris::DDRMat >   const &  aFirstTRICoords,
                moris::Matrix < moris::DDRMat >   const &  aSecondTRICoords,
                moris::Matrix < moris::IdMat >    const &  aFirstTRINodeIndex,
                moris::Matrix < moris::IdMat >    const &  aSecondTRINodeIndex)
        {

            // Initialize the output
            moris::Cell< moris::Matrix < moris::DDRMat > > allIntersectionTriangles ;

            // Loop over second mesh
            for( uint i = 0 ; i < aSecondTRINodeIndex.n_rows(); i++ )
            {
                // Loop over the first mesh
                for( uint j = 0 ; j < aFirstTRINodeIndex.n_rows(); j++)
                {
                    // Get the node id of the examined triangles
                    Matrix < IdMat > Tbbc = aSecondTRINodeIndex( {i,i}, {0, 2});
                    Matrix < IdMat > Taac = aFirstTRINodeIndex( {j,j}, {0, 2});

                    // Initialize the coordinate matrices
                    Matrix < DDRMat > NbTbbc( 2, 3 );
                    Matrix < DDRMat > NaTaac( 2, 3 );

                    // Assign coordinates
                    for ( uint i = 0 ; i < 3 ; i++ )
                    {
                        NbTbbc.get_column( i ) = aSecondTRICoords.get_column( Tbbc( i ) - 1 );
                        NaTaac.get_column( i ) = aFirstTRICoords.get_column( Taac( i ) - 1 );
                    }

                    // Initialize the polygon clipping outputs
                    moris::Matrix < moris::DDRMat  > tP;

                    moris::Matrix < moris::DDUMat  >  tnc;

                    // find the intersection of 2 element triangulation
                    this->Intersect( NbTbbc, NaTaac, tP, tnc ) ;

                    // If it is a polygon add to the output
                    if (tP.n_cols()>2)
                    {
                        allIntersectionTriangles.push_back(tP) ;
                    }
                }
            }

            //return the output
            return allIntersectionTriangles;
        }

        //------------------------------------------------------------------------------

        void
        Intersection_Detect::make_new_pairs( moris::Cell< Matrix<DDRMat> > tP,
                moris::mtk::Cell const & aInterpCell1,
                moris::mtk::Cell const & aInterpCell2,
                uint aPairCount,
                moris_index aPhaseToPhase)
        {


            //convert to desired form
            moris::Matrix<DDRMat >     tRotation;
            moris::Matrix<DDRMat >     tInverseRotation;
            this->rotation_matrix( tRotation, tInverseRotation, aPairCount );


            //make a custom map for every coordinate system
            //requires that we give a unique coordinates and associated index to each of them
            Matrix< DDRMat > tIntersectedPoints;

            //join all the coordiantes together
            for (size_t i = 0 ; i < tP.size(); i++ )
            {
                tIntersectedPoints = join_horiz( tIntersectedPoints, tP( i ) );
            }

            //make a cell of unique coordinates
            moris::Cell< Matrix<DDRMat> > tUniqueIntersectedPointsCell;

            //put the first element
            tUniqueIntersectedPointsCell.push_back(tIntersectedPoints.get_column(0));

            //iterate through other coordinates and add them to unique coordinate cell if not already there
            for(uint i = 1 ; i < tIntersectedPoints.n_cols() ; i++ )
            {
                uint tRepeatCount = 0 ;
                for( uint j = 0 ; j < tUniqueIntersectedPointsCell.size() ; j++ )
                {
                    if( norm( tIntersectedPoints.get_column(i) - tUniqueIntersectedPointsCell(j).get_column(0) ) < 0.00000000001 )
                    {
                        tRepeatCount++;
                        break;
                    }
                }

                if( tRepeatCount == 0 )
                {
                    tUniqueIntersectedPointsCell.push_back(tIntersectedPoints.get_column(i));
                }
            }

            //transfer coordinates of cells to a matrix
            Matrix<DDRMat> tUniqueIntersectedPoints( 2 ,tUniqueIntersectedPointsCell.size() );

            for( uint k = 0 ; k <tUniqueIntersectedPointsCell.size() ; k++ )
            {
                tUniqueIntersectedPoints.get_column(k) =  tUniqueIntersectedPointsCell(k).get_column(0);
            }

            //intersection cell from the algorithm comes as coordinates, convert them to unique coordinate indices
            moris::Cell< Matrix< IndexMat > > tPVertexIndex( tP.size() ) ;

            for( size_t i = 0 ; i < tP.size() ; i++ )
            {
                tPVertexIndex(i) = Matrix<IndexMat>(1, tP(i).n_cols());

                for( size_t j = 0; j < tP(i).n_cols(); j++ )
                {
                    for( size_t k = 0 ; k < tUniqueIntersectedPoints.n_cols() ; k++ )
                    {
                        if ( norm( tP(i).get_column(j) - tUniqueIntersectedPoints.get_column(k) ) < 0.00000000001 )
                        {
                            tPVertexIndex(i)(j) = k ;
                        }
                    }
                }
            }

            //create master and slave vertices
            moris::Cell<moris::mtk::Vertex* > tMasterVerticesCell = this->create_master_vertices(tUniqueIntersectedPoints,aInterpCell1, aPairCount);
            moris::Cell<moris::mtk::Vertex* > tSlaveVerticesCell = this->create_slave_vertices(tUniqueIntersectedPoints,aInterpCell2, aPairCount);

            //obtain the integration mesh
            moris::mtk::Integration_Mesh* tIntegrationMesh = mMeshManager->get_integration_mesh( mMeshIndex );

            // Get number of nodes
            uint tNumExistingNodes = tIntegrationMesh->get_num_entities(moris::EntityRank::NODE);
            uint tNewNodes = tMasterVerticesCell.size()+ tSlaveVerticesCell.size()+ mNewNodeCoords.size();

            // Resize the local to global map for the vertices
            mEntityLocaltoGlobalMap(0).resize(tNumExistingNodes + tNewNodes);

            // Offset due to addition of vertices at each loop
            uint tMasterOffset  = mNewNodeCoords.size();
            uint tSlaveOffset  = mNewNodeCoords.size() + tMasterVerticesCell.size();

            // Assign ids and indices to the maps
            for( uint Ii = 0 ; Ii < tMasterVerticesCell.size() ; Ii++ )
            {

                mNodeIndexToCoordsMap[ tMasterVerticesCell(Ii)->get_index() ] = Ii + tMasterOffset ;

                // Assign Coordinates
                Matrix<DDRMat> tCoords = tMasterVerticesCell(Ii)->get_coords();
                mNewNodeCoords.push_back(tCoords);

                mEntityLocaltoGlobalMap(0)(tMasterVerticesCell(Ii)->get_index()) = tMasterVerticesCell(Ii)->get_id();

            }

            for( uint Ii = 0 ; Ii < tSlaveVerticesCell.size() ; Ii++ )
            {
                mNodeIndexToCoordsMap[ tSlaveVerticesCell(Ii)->get_index() ] = Ii + tSlaveOffset;

                // Assign Coordinates
                mNewNodeCoords.push_back(tSlaveVerticesCell(Ii)->get_coords());

                mEntityLocaltoGlobalMap(0)(tSlaveVerticesCell(Ii)->get_index()) = tSlaveVerticesCell(Ii)->get_id();
            }


            //initialize master,slave side IG cells
            moris::Cell<moris::mtk::Cell const *> tMasterIntegCells;
            moris::Cell<moris::mtk::Cell const *> tSlaveIntegCells;

            //Added surfaces if intersection area is a polygon
            uint tAddedSurafceNum = 0 ;

            //iterate through each intersection surface and make Ig cells
            for( size_t tClusterNum = 0 ; tClusterNum < tP.size() ; tClusterNum++ )
            {
                //if polygon break it further to tri by fan triangulation
                if ( tPVertexIndex(tClusterNum).n_cols() > 3 )
                {
                    uint tNumEdges = tPVertexIndex(tClusterNum).n_cols();

                    uint tNumTriangles = tNumEdges - 2 ;

                    tAddedSurafceNum += ( tNumTriangles - 1 ) ;

                    for(uint tTRI = 0 ; tTRI < tNumTriangles ; tTRI++ )
                    {
                        Matrix< IndexMat > tTmpVertexIndex = { {tPVertexIndex(tClusterNum)(0),tPVertexIndex(tClusterNum)(tTRI+1),tPVertexIndex(tClusterNum)(tTRI+2) } };

                        moris::mtk::Cell const * tMasterIgCell = this->create_master_ig_cell( tMasterVerticesCell, tTmpVertexIndex, aInterpCell1,aPairCount ) ;
                        tMasterIntegCells.push_back( tMasterIgCell );

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
            moris::Matrix<moris::IndexMat>          tMasterIntegrationCellSideOrdinals( 1, tP.size() + tAddedSurafceNum, 3);
            moris::Matrix<moris::IndexMat>          tSlaveIntegrationCellSideOrdinals( 1, tP.size() + tAddedSurafceNum, 3);

            //parametric coordiate matrix for slave and master vertices
            moris::Matrix< DDRMat > tMasterParamCoords( tUniqueIntersectedPoints.n_cols(), 3 );
            moris::Matrix< DDRMat > tSlaveParamCoords( tUniqueIntersectedPoints.n_cols(), 3 ) ;

            for(uint  i = 0 ; i < tUniqueIntersectedPoints.n_cols() ; i++ )
            {
                moris::Matrix< DDRMat > tMasterParamCoord = { { -1,  tUniqueIntersectedPoints( 0, i ), tUniqueIntersectedPoints( 1, i ) } };
                moris::Matrix< DDRMat > tSlaveParamCoord = { { +1,  tUniqueIntersectedPoints( 0, i ), tUniqueIntersectedPoints( 1, i ) } };
                tMasterParamCoord = tMasterParamCoord*tRotation;
                tSlaveParamCoord  = tSlaveParamCoord*tRotation;
                tMasterParamCoords.get_row( i ) = tMasterParamCoord.get_row(0);
                tSlaveParamCoords.get_row( i )  = tSlaveParamCoord.get_row(0);
            }

            //create constant version of the vertices on cluster
            moris::Cell< moris::mtk::Vertex const *> tMasterVerticesConst;
            moris::Cell< moris::mtk::Vertex const *> tSlaveVerticesConst;

            for(uint i = 0 ; i < tMasterVerticesCell.size() - 1 ; i++)
            {
                tMasterVerticesConst.push_back( tMasterVerticesCell( i ) );
                tSlaveVerticesConst.push_back( tSlaveVerticesCell( i ) );
            }

            //store the vertices to prevent memory leak
            mMasterVertices.append(tMasterVerticesConst);
            mSlaveVertices.append(tSlaveVerticesConst);

            //create the master side cluster
            moris::mtk::Side_Cluster_ISC* tMasterSideCluster = new moris::mtk::Side_Cluster_ISC(false,
                    & aInterpCell1,
                    tMasterIntegCells,
                    tMasterIntegrationCellSideOrdinals,
                    tMasterVerticesConst,
                    tMasterParamCoords);

            //create the slave side cluster
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

            //store data to be used further
            mDoubleSidedClusters.push_back(tDblSideCluster);
            mDoubleSidedClustersIndex.push_back(aPhaseToPhase);

            //master
            mMasterSidedClusters.push_back(tMasterSideCluster);

            //slave
            mSlaveSidedClusters.push_back(tSlaveSideCluster);

            // Get number of elements
            uint tNumExistingElem = tIntegrationMesh->get_num_entities(moris::EntityRank::ELEMENT);//obtain the integration mesh
            uint tNewElem =  tMasterIntegCells.size() + tSlaveIntegCells.size()+ mMasterSideCells.size()*2 ;

            // Resize number of elements
            mEntityLocaltoGlobalMap(3).resize(tNumExistingElem + tNewElem);

            // Offset due to addition of new mtk cells at each loop
            uint tOffsetCells  = mMasterSideCells.size();

            // Append the integration cells
            mMasterSideCells.append(tMasterIntegCells);
            mSlaveSideCells.append(tSlaveIntegCells);

            // Create the maps for cell
            for(uint Ii = 0 ; Ii < tMasterIntegCells.size() ; Ii++)
            {
                mMasterCellIndextoCellMap[ tMasterIntegCells(Ii)->get_index() ] = Ii+ tOffsetCells ;

                mEntityLocaltoGlobalMap(3)( tMasterIntegCells(Ii)->get_index() ) =  tMasterIntegCells(Ii)->get_id();

                mEntityLocaltoGlobalMap(3)( tSlaveIntegCells(Ii)->get_index() ) =  tMasterIntegCells(Ii)->get_id();
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
            moris::Cell<moris::mtk::Vertex *> tIgCellVertices;

            //surface vertices
            for(uint i = 0 ; i < tVertexIndex.n_cols() ; i++ )
            {
                tIgCellVertices.push_back( tMasterVertices( tVertexIndexCounterClockWise ( i ) ) );
            }

            //head node of tet4
            tIgCellVertices.push_back( tMasterVertices( tMasterVertices.size() - 1 ) );

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
            moris::Cell<moris::mtk::Vertex *> tIgCellVertices;

            //surface vertices
            for(uint i = 0 ; i < tVertexIndex.n_cols() ; i++ )
            {
                tIgCellVertices.push_back( tSlaveVertices( tVertexIndexCounterClockWise ( i ) ) );
            }

            //head node of tet4
            tIgCellVertices.push_back( tSlaveVertices( tSlaveVertices.size() - 1 ) );

            //allocate id and index
            moris_index tCellIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::ELEMENT) ;
            moris_id tCellId = mIntersectedMeshData.allocate_entity_ids_external_entity_data(1, EntityRank::ELEMENT);

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
        void Intersection_Detect::constrcuct_add_dbl_sided_set (moris::Matrix < IndexMat > tPhaseInteractionTable)
        {
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

                    mSideSetLabels.push_back(tDoubleSideSetName);

                }
            }

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
            }

            //populate the second side clusters
            for (uint i = 0; i < tSecondSideSetNames.size() ; i++)
            {
                moris::mtk::Set* tSet = tIntegrationMesh->get_set_by_name( tSecondSideSetNames(i) );

                // get clusters in the second set
                moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet->get_clusters_on_set();

                tSecondSideClusters.append(tSetClusters);
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

        moris::Cell< moris::mtk::Vertex *>
        Intersection_Detect::create_master_vertices(Matrix < DDRMat> tUniqueIntersectedPoints,moris::mtk::Cell const & aMasterInterpCell, uint aPairCount )
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
            moris::Matrix<moris::DDRMat> tCoordinates( tElementNodes.numel() , 3);

            for (uint  i = 0 ; i < tElementNodes.numel() ; i++ )
            {
                moris::Matrix<moris::DDRMat> tmp =  tInterpMesh->get_node_coordinate( tElementNodes( i ) );
                tCoordinates.get_row(i) = tmp.get_row(0);
            }

            for( uint i = 0 ; i < tNumSurfaceNodes ; i++)
            {
                //form the parametric coordinates of the nodes
                moris::Matrix<DDRMat> tParamCoordsRelativeToElem = {{-1.0,tUniqueIntersectedPoints(0,i),tUniqueIntersectedPoints(1,i)}};

                tParamCoordsRelativeToElem = tParamCoordsRelativeToElem* tRotation ;

                //obtain physical coordinates of the node
                Matrix <DDRMat > tNewCoords;
                moris::mtk::Interpolation::trilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem, tNewCoords);

                //update the vertex index
                moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;

                //construct the new vertex
                tVertcies ( i ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex,tNewCoords ) ;

                // Update the index and id
                mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );
                tVertexId++;
            }

            //center node coordinates
            moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
            Matrix <DDRMat > tParamCoordsRelativeToElem = {{0.0, 0.0, 0.0}};
            moris::mtk::Interpolation::trilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem, tNewNodeCoordinates);


            moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;
            mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );

            //add the top vertex of tet4
            tVertcies( tNumSurfaceNodes ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewNodeCoordinates) ;

            return tVertcies;
        }

        // ----------------------------------------------------------------------------

        moris::Cell< moris::mtk::Vertex *>
        Intersection_Detect::create_slave_vertices(Matrix < DDRMat> tUniqueIntersectedPoints,moris::mtk::Cell const & aSlaveInterpCell, uint aPairCount )
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
            moris::Matrix<moris::DDRMat> tCoordinates( tElementNodes.numel() , 3);

            for (uint  i = 0 ; i < tElementNodes.numel() ; i++ )
            {
                moris::Matrix<moris::DDRMat> tmp =  tInterpMesh->get_node_coordinate( tElementNodes( i ) );
                tCoordinates.get_row(i) = tmp.get_row(0);
            }

            for( uint i = 0 ; i < tNumSurfaceNodes ; i++)
            {
                //form the parametric coordinates of the nodes
                moris::Matrix<DDRMat> tParamCoordsRelativeToElem = {{+1.0,tUniqueIntersectedPoints(0,i),tUniqueIntersectedPoints(1,i)}};

                //convert to desired form
                tParamCoordsRelativeToElem = tParamCoordsRelativeToElem * tRotation;

                //obtain physical coordinates of the node
                Matrix <DDRMat > tNewCoords;
                moris::mtk::Interpolation::trilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem, tNewCoords);

                moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;

                //construct the new vertex
                tVertcies ( i ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex,tNewCoords ) ;

                mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );
                tVertexId++;
            }

            // center node coordinates
            moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
            Matrix <DDRMat > tParamCoordsRelativeToElem = {{0.0, 0.0, 0.0}};
            moris::mtk::Interpolation::trilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem, tNewNodeCoordinates);

            moris_index tVertexIndex = mIntersectedMeshData.get_first_available_index_external_data( EntityRank::NODE) ;
            mIntersectedMeshData.update_first_available_index_external_data( tVertexIndex+1, EntityRank::NODE );

            //add the top vertex of tet4
            tVertcies( tNumSurfaceNodes ) = new moris::mtk::Vertex_ISC( tVertexId, tVertexIndex, tNewNodeCoordinates) ;

            return tVertcies;
        }

        // ----------------------------------------------------------------------------

        void
        Intersection_Detect::rotation_matrix(moris::Matrix< DDRMat > &aRotation, moris::Matrix< DDRMat > &aInverseRotation, uint aPairCount )
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
                moris::Matrix< DDRMat > const &               aOffsetMatrix)
        {
            // right interpolation cell
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

                    for ( uint k = 0 ; k < 3 ; k++ )
                    {
                        NbTbbc.get_column( k ) = aSecondTRICoords.get_column( Tbbc( k ) - 1 );
                        NaTaac.get_column( k ) = aFirstTRICoords.get_column( Taac( k ) - 1 );
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
                for (uint i = 0 ; i < al.size(); i++ )
                {
                    ad ( al (i) - 1) = 1.0 ;
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

                    for ( uint i = 0 ; i < 3 ; i++ )
                    {
                        NbTbbc.get_column( i ) = aSecondTRICoords.get_column( Tbbc( i ) - 1 );
                        NaTaac.get_column( i ) = aFirstTRICoords.get_column( Taac( i ) - 1 );
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

                        for (uint i = 0 ; i < Taac4to6.n_cols() ; i++ )
                        {
                            tmp(i) = ad( Taac4to6(i) - 1);
                        }

                        Matrix < DDUMat> comp( 1, tmp.n_cols() ,0.0 );

                        Matrix < DDNIMat> tn = find( tmp == comp ) ;

                        if ( tn.numel() != 0)
                        {

                            Matrix < IdMat> t(1,tn.numel()) ;

                            for ( uint i = 0 ; i < t.n_cols() ; i++ )
                            {
                                t ( i ) = aFirstTRIConnect(ac, tn(i) );
                            }

                            //add neighbours
                            for ( uint i = 0; i < t.n_cols() ; i++ )
                            {
                                al.push_back( t( i ) );
                            }

                            for ( uint i = 0; i < t.n_cols() ; i++ )
                            {
                                ad( t( i ) - 1) = 1.0 ;
                            }
                        }

                        // n(find(nc>0))  = ac;
                        Matrix < DDUMat > comp2 ( 1, tnc.n_cols(), 0 );

                        Matrix < DDNIMat > findnc = find(  comp2 < tnc) ;

                        for (uint i = 0 ; i < findnc.numel() ; i++ )
                        {
                            n( findnc( i )  ) = ac + 1 ;
                        }
                    }
                }

                // find non-treated neighbors
                // tmp=find(bd(Tb(bc,4:6))==0)
                Matrix < IdMat> Tbbc4to6 =  aSecondTRIConnect( {bc,bc}, {0, 2});

                Matrix < DDUMat> bdTbbc( 1, Tbbc4to6.n_cols() ) ;

                for (uint i = 0 ; i < Tbbc4to6.n_cols() ; i++ )
                {
                    bdTbbc( i ) = bd( Tbbc4to6( i ) - 1);
                }

                Matrix < DDUMat> zerovec( 1, bdTbbc.n_cols() ,0.0 );

                Matrix < DDNIMat > tmp = find( bdTbbc == zerovec );

                if (tmp.numel() != 0 )
                {
                    // take those which intersect
                    // idx=find(n(tmp)>0)
                    Matrix < DDUMat> ntmp( 1, tmp.n_rows() );

                    for (uint i = 0 ; i < ntmp.n_cols() ; i++ )
                    {
                        ntmp( i ) = n( tmp ( i ) );
                    }

                    Matrix < DDUMat> zerovec2( 1, ntmp.n_cols() ,0.0 );

                    Matrix < DDNIMat > idx = find ( zerovec2 < ntmp );

                    if ( idx.numel() != 0 )
                    {

                        // t=Tb(bc,3+tmp(idx))
                        Matrix < DDUMat> tmpidx(1, idx.n_rows() ) ;

                        for( uint i = 0 ; i < tmpidx.n_cols() ; i++ )
                        {
                            tmpidx (i) = tmp ( idx ( i ) ) ;
                        }

                        Matrix < IdMat > t(1,tmpidx.n_cols() ) ;

                        for ( uint i = 0 ; i < t.n_cols() ; i++ )
                        {
                            t( i ) =  aSecondTRIConnect (bc, tmpidx (i)) ;
                        }

                        //and add the ones that intersect
                        for ( uint i = 0 ; i < t.n_cols() ; i++ )
                        {
                            bl.push_back( t( i ) ) ;
                        }

                        //with starting candidates Ta
                        // bil=[bil n(tmp(idx))]
                        if ( tmpidx.numel() != 0 )
                        {
                            Matrix < DDUMat> ntmpidx(1, idx.numel() ) ;

                            for(uint i = 0 ; i < ntmpidx.n_cols() ; i++ )
                            {
                                ntmpidx( i ) = n ( tmpidx (i) );
                            }

                            for(uint i = 0 ; i < ntmpidx.n_cols() ; i++ )
                            {
                                bil.push_back( ntmpidx( i ) );
                            }

                            for ( uint i = 0 ; i < t.n_cols() ; i++ )
                            {
                                bd( t( i ) - 1 ) = 1 ;
                            }
                        }
                    }
                }
            }

            //return the cell of intersected
            return allIntersectionTriangles;
        }
    }

}
