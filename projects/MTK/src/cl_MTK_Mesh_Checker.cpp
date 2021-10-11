/*
 * cl_MTK_Mesh_Checker.cpp
 *
 *  Created on: Jan 15, 2020
 *      Author: doble
 */

#include "cl_MTK_Mesh_Checker.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_all_true.hpp"
#include "fn_norm.hpp"
#include "fn_unique.hpp"
#include "op_minus.hpp"

#include <unordered_map>

#include "cl_MPI_Tools.hpp"

namespace moris
{

namespace mtk
{

    moris_index
    Serialized_Mesh_Data::get_vertex_index(
        const moris_index& aVertexId,
        const moris_index& aProcIndex ) const
    {
        MORIS_ERROR( par_rank() == 0, "ONLY VALID ON PROC 0" );
        auto tIter = mCollectVertexMaps( aProcIndex ).find( aVertexId );
        MORIS_ERROR( tIter != mCollectVertexMaps( aProcIndex ).end(), "Id not in map" );

        return tIter->second;
    }

    //
    moris::Matrix< moris::DDRMat >
    Serialized_Mesh_Data::get_vertex_coords(
        const moris_index& aVertexId,
        const moris_index& aProcIndex ) const
    {
        MORIS_ERROR( par_rank() == 0, "ONLY VALID ON PROC 0" );

        moris_index tIndex = get_vertex_index( aVertexId, aProcIndex );


        return mCollectVertexCoords( aProcIndex ).get_row( tIndex );
    }

    //--------------------------------------------------------------------------------
    Mesh_Checker::Mesh_Checker() {}
    //--------------------------------------------------------------------------------
    Mesh_Checker::Mesh_Checker(
        moris_index         aMeshIndex,
        Interpolation_Mesh* aIpMesh,
        Integration_Mesh*   aIgMesh ) :
        mMeshIndex( aMeshIndex ),
        mIpMesh( aIpMesh ),
        mIgMesh( aIgMesh )
    {
        this->serialize_mesh();

        // core mesh gather
        this->gather_serialized_mesh( &mSerializedIpMesh );
        this->gather_serialized_mesh( &mSerializedIgMesh );

        // ip specific gathering
        this->gather_serialized_ip_mesh( &mSerializedIpMesh );

        // setup maps for each proc
        // vertex maps
        this->setup_vertex_maps( &mSerializedIpMesh );
        this->setup_vertex_maps( &mSerializedIgMesh );

        // Cell maps
        this->setup_cell_maps( &mSerializedIpMesh );
        this->setup_cell_maps( &mSerializedIgMesh );
    }
    //--------------------------------------------------------------------------------
    Mesh_Checker::~Mesh_Checker() {}
    //--------------------------------------------------------------------------------
    bool
    Mesh_Checker::perform()
    {
        // verify the coordinates
        mIpVertexDiag = this->verify_vertex_coordinates( &mSerializedIpMesh, true );
        mIgVertexDiag = this->verify_vertex_coordinates( &mSerializedIgMesh, true );

        // verify the ownership
        mIpVertexOwnerDiag = this->verify_vertex_ownership( &mSerializedIpMesh );
        mIgVertexOwnerDiag = this->verify_vertex_ownership( &mSerializedIgMesh );

        // verify the ownership
        mIpCellOwnerDiag = this->verify_cell_ownership( &mSerializedIpMesh );
        mIgCellOwnerDiag = this->verify_cell_ownership( &mSerializedIgMesh );


        return false;
    }

    //--------------------------------------------------------------------------------
    bool
    Mesh_Checker::verify_double_side_sets( Integration_Mesh const* aIgMesh )
    {
        moris::uint tNumDoubleSideSets = aIgMesh->get_num_double_side_set();

        // iterate through the double side sets
        for ( moris::uint iSet = 0; iSet < tNumDoubleSideSets; iSet++ )
        {
            moris::Cell< Cluster const* > tClustersInSet = aIgMesh->get_double_side_set_cluster( (moris_index)iSet );

            // iterate through clusters in set
            for ( auto iCl : tClustersInSet )
            {
                // verify the master
                bool tMasterValid = this->verify_side_cluster( iCl, Master_Slave::MASTER );

                if ( !tMasterValid )
                {
                    MORIS_LOG_ERROR( "\nInvalid master cluster" );
                }

                // verify the slave
                bool tSlaveValid = this->verify_side_cluster( iCl, Master_Slave::SLAVE );

                if ( !tSlaveValid )
                {
                    MORIS_LOG_ERROR( "\nInvalid master cluster" );
                }

                if ( !tMasterValid || !tSlaveValid )
                {
                    return false;
                }
            }
        }

        return true;
    }
    //--------------------------------------------------------------------------------
    bool
    Mesh_Checker::verify_side_cluster( Cluster const* aCluster,
        enum Master_Slave                             aMasterSlave )
    {
        // check if trivial
        bool tTrivial = aCluster->is_trivial( aMasterSlave );

        if ( tTrivial )
        {
            // get the interpolation cell
            mtk::Cell const& tIpCell = aCluster->get_interpolation_cell( aMasterSlave );

            // get the integration primary cells
            moris::Cell< moris::mtk::Cell const* > const& tIgCells = aCluster->get_primary_cells_in_cluster( aMasterSlave );

            if ( tIgCells.size() != 1 )
            {
                MORIS_LOG_ERROR( "\nTrivial cluster needs to have 1 primary integration cell." );
                return false;
            }

            // get the side ordinals
            moris::Matrix< moris::IndexMat > tSideOrd = aCluster->get_cell_side_ordinals( aMasterSlave );

            if ( tSideOrd.numel() != 1 )
            {
                MORIS_LOG_ERROR( "\nTrivial cluster needs to have exactly 1 side ordinal." );
                return false;
            }

            moris::Cell< moris::mtk::Vertex const* > tIpVertsOnSide = tIpCell.get_vertices_on_side_ordinal( tSideOrd( 0 ) );
            moris::Cell< moris::mtk::Vertex const* > tIgVertsOnSide = tIgCells( 0 )->get_vertices_on_side_ordinal( tSideOrd( 0 ) );

            if ( tIpVertsOnSide.size() != tIgVertsOnSide.size() )
            {
                MORIS_LOG_ERROR( "\nTrivial cluster ip and ig vertex on side mismatch." );
                return false;
            }


            // check the coordinates
            for ( moris::uint iV = 0; iV < tIpVertsOnSide.size(); iV++ )
            {
                if ( moris::norm( tIpVertsOnSide( iV )->get_coords() - tIgVertsOnSide( iV )->get_coords() ) > 1e-8 )
                {
                    return false;
                }
            }
        }
        else
        {

            // get the interpolation cell
            mtk::Cell const& tIpCell = aCluster->get_interpolation_cell( aMasterSlave );

            // create cell info
            Cell_Info_Factory      tCellInfoFactory;
            moris::mtk::Cell_Info* tCellInfo = tCellInfoFactory.create_cell_info( tIpCell.get_geometry_type(), tIpCell.get_interpolation_order() );

            // get the vertices in the cluster
            moris::Cell< moris::mtk::Vertex const* > tVertsInCluster = aCluster->get_vertices_in_cluster( aMasterSlave );

            // ip verts
            Matrix< DDRMat > tIpCoords = tIpCell.get_vertex_coords();

            for ( size_t i = 0; i < tVertsInCluster.size(); i++ )
            {
                // local coordinates
                moris::Matrix< moris::DDRMat > tLocalCoords = aCluster->get_vertex_local_coordinate_wrt_interp_cell( tVertsInCluster( i ), aMasterSlave );

                // Evalute the basis function at the point
                moris::Matrix< moris::DDRMat > tN;
                tCellInfo->eval_N( tLocalCoords, tN );

                // Evaluate the nodes global coordinate from the basis weights
                moris::Matrix< moris::DDRMat > tInterpedNodeCoord = tN * tIpCoords;

                // Verify the interpolated coordinate is equal to the node coordinate row
                if ( norm( tInterpedNodeCoord - tVertsInCluster( i )->get_coords() ) > 1e-8 )
                {
                    return false;
                }
            }

            delete tCellInfo;
        }

        return true;
    }
    //--------------------------------------------------------------------------------
    bool
    Mesh_Checker::verify_vertex_coordinates(
        Serialized_Mesh_Data* aSerializedMesh,
        bool                  aStackedVertexFlag )
    {
        moris::uint tValidInd = 1;

        if ( par_rank() == 0 )
        {

            Matrix< IndexMat > tAllNodeMaps   = this->concatenate_cell_of_mats( aSerializedMesh->mCollectVertexIds, 1 );
            Matrix< DDRMat >   tAllNodeCoords = this->concatenate_cell_of_mats( aSerializedMesh->mCollectVertexCoords, 1 );


            // check coordinates
            std::unordered_map< moris_id, moris_index > tNodeIdToIndMap;
            for ( moris::uint i = 0; i < tAllNodeMaps.numel(); i++ )
            {
                auto tIter = tNodeIdToIndMap.find( tAllNodeMaps( i ) );
                if ( tIter == tNodeIdToIndMap.end() )
                {
                    tNodeIdToIndMap[tAllNodeMaps( i )] = i;
                }

                else
                {
                    moris_index tIndex = tIter->second;
                    if ( moris::norm( tAllNodeCoords.get_row( i ) - tAllNodeCoords.get_row( tIndex ) ) > 1e-6 )
                    {
                        tValidInd = 0;
                    }
                }
            }
        }


        broadcast( tValidInd );

        barrier();

        if ( tValidInd == 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    bool
    Mesh_Checker::verify_vertex_ownership( Serialized_Mesh_Data* aSerializedMesh )
    {
        moris::uint tValidInd = 1;

        if ( par_rank() == 0 )
        {
            Matrix< IndexMat > tVertexOwners( 1, aSerializedMesh->mVertexSerialIndexToId.size(), MORIS_INDEX_MAX );
            // iterate through processors
            for ( moris::uint iP = 0; iP < aSerializedMesh->mCollectVertexIds.size(); iP++ )
            {

                // iterate through vertices on proc
                for ( moris::uint iV = 0; iV < aSerializedMesh->mCollectVertexIds( iP ).numel(); iV++ )
                {
                    moris_id tVertexId = aSerializedMesh->mCollectVertexIds( iP )( iV );

                    moris_id tSerialVertInd = this->get_vertex_serial_index_from_id( tVertexId, aSerializedMesh );
                    moris_id tProcVertInd   = this->get_vertex_proc_index_from_id( tVertexId, iP, aSerializedMesh );


                    if ( tVertexOwners( tSerialVertInd ) == MORIS_INDEX_MAX )
                    {
                        tVertexOwners( tSerialVertInd ) = aSerializedMesh->mCollectVertexOwners( iP )( tProcVertInd );
                    }
                    else
                    {
                        // ownership issue
                        if ( aSerializedMesh->mCollectVertexOwners( iP )( tProcVertInd ) != tVertexOwners( tSerialVertInd ) )
                        {
                            std::cout << "Issue with Vertex Ownership, Vertex ID = " << std::setw( 12 ) << tVertexId << " | Owner on " << std::setw( 4 ) << iP << " = ";
                            std::cout << aSerializedMesh->mCollectVertexOwners( iP )( tProcVertInd ) << " | expected owner = " << std::setw( 4 ) << tVertexOwners( tSerialVertInd );
                            std::cout << std::endl;
                            moris::print( aSerializedMesh->get_vertex_coords( tVertexId, tVertexOwners( tSerialVertInd ) ), "Coords on expected owner" );
                            moris::print( aSerializedMesh->get_vertex_coords( tVertexId, iP ), "Coords on current" );

                            tValidInd = 0;
                        }
                    }
                }
            }

            // iterate to see if a proc thinks someone owns the vert but they dont
            for ( moris::uint iV = 0; iV < tVertexOwners.numel(); iV++ )
            {
                moris_index tOwner       = tVertexOwners( iV );
                moris_index tID          = aSerializedMesh->mVertexSerialIndexToId( iV );
                moris_id    tProcVertInd = this->get_vertex_proc_index_from_id( tID, tOwner, aSerializedMesh );

                if ( aSerializedMesh->mCollectVertexOwners( tOwner )( tProcVertInd ) != tOwner )
                {
                    std::cout << "Vertex Ownership issue 2" << std::endl;
                    tValidInd = 0;
                }
            }
        }

        broadcast( tValidInd );

        barrier();

        if ( tValidInd == 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    bool
    Mesh_Checker::verify_cell_ownership( Serialized_Mesh_Data* aSerializedMesh )
    {
        moris::uint tValidInd = 1;

        if ( par_rank() == 0 )
        {
            Matrix< IndexMat > tCellOwners( 1, aSerializedMesh->mCellSerialIndexToId.size(), MORIS_INDEX_MAX );
            // iterate through processors
            for ( moris::uint iP = 0; iP < aSerializedMesh->mCollectCellIds.size(); iP++ )
            {
                // iterate through vertices on proc
                for ( moris::uint iV = 0; iV < aSerializedMesh->mCollectCellIds( iP ).numel(); iV++ )
                {
                    moris_id tCellId = aSerializedMesh->mCollectCellIds( iP )( iV );


                    moris_id tSerialVertInd = this->get_cell_serial_index_from_id( tCellId, aSerializedMesh );
                    moris_id tProcVertInd   = this->get_cell_proc_index_from_id( tCellId, iP, aSerializedMesh );

                    if ( tCellOwners( tSerialVertInd ) == MORIS_INDEX_MAX )
                    {
                        tCellOwners( tSerialVertInd ) = aSerializedMesh->mCollectCellOwners( iP )( tProcVertInd );
                    }
                    else
                    {
                        // ownership issue
                        if ( aSerializedMesh->mCollectCellOwners( iP )( tProcVertInd ) != tCellOwners( tSerialVertInd ) )
                        {
                            std::cout << "Issue with Cell Ownership, Cell ID = " << tCellId << " | Owner on " << iP << " = ";
                            std::cout << aSerializedMesh->mCollectCellOwners( iP )( tProcVertInd ) << " | expected owner = " << tCellOwners( tSerialVertInd ) << std::endl;
                            tValidInd = 0;
                        }
                    }
                }
            }

            // iterate to see if a proc thinks someone owns the vert but they dont
            for ( moris::uint iV = 0; iV < tCellOwners.numel(); iV++ )
            {
                moris_index tOwner       = tCellOwners( iV );
                moris_index tID          = aSerializedMesh->mCellSerialIndexToId( iV );
                moris_id    tProcVertInd = this->get_cell_proc_index_from_id( tID, tOwner, aSerializedMesh );

                if ( aSerializedMesh->mCollectCellOwners( tOwner )( tProcVertInd ) != tOwner )
                {
                    std::cout << "Cell Ownership issue 2" << std::endl;
                    tValidInd = 0;
                }
            }
        }


        broadcast( tValidInd );

        barrier();

        if ( tValidInd == 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    //--------------------------------------------------------------------------------
    void
    Mesh_Checker::print_diagnostics()
    {
        if ( par_rank() == 0 )
        {
            uint tFirstColWidth  = 15;
            uint tSecondColWidth = 11;

            std::cout << "Mesh Checking Diagnostics:" << std::endl;
            std::cout << std::setw( tFirstColWidth ) << "" << std::setw( tSecondColWidth ) << "Interp Mesh" << std::setw( tSecondColWidth ) << " | Integ Mesh" << std::endl;
            std::cout << std::setw( tFirstColWidth ) << "Vertex Ids: " << std::setw( tSecondColWidth ) << this->bool_to_string( mIpVertexDiag ) << " |" << std::setw( tSecondColWidth ) << this->bool_to_string( mIgVertexDiag ) << std::endl;
            std::cout << std::setw( tFirstColWidth ) << "Vertex Coords: " << std::setw( tSecondColWidth ) << this->bool_to_string( mIpVertexDiag ) << " |" << std::setw( tSecondColWidth ) << this->bool_to_string( mIgVertexDiag ) << std::endl;
            std::cout << std::setw( tFirstColWidth ) << "Vertex Owners: " << std::setw( tSecondColWidth ) << this->bool_to_string( mIpVertexOwnerDiag ) << " |" << std::setw( tSecondColWidth ) << this->bool_to_string( mIgVertexOwnerDiag ) << std::endl;
            std::cout << std::setw( tFirstColWidth ) << "Cell Owners: " << std::setw( tSecondColWidth ) << this->bool_to_string( mIpCellOwnerDiag ) << " |" << std::setw( tSecondColWidth ) << this->bool_to_string( mIgCellOwnerDiag ) << std::endl;
        }
    }
    //--------------------------------------------------------------------------------

    moris_index
    Mesh_Checker::get_vertex_proc_index_from_id( moris_id aId,
        moris_index                                       aProc,
        Serialized_Mesh_Data*                             aSerializedMesh )
    {
        MORIS_ERROR( aSerializedMesh->mCollectVertexMaps( aProc ).find( aId ) != aSerializedMesh->mCollectVertexMaps( aProc ).end(), "Vertex id not in the map" );
        return aSerializedMesh->mCollectVertexMaps( aProc )[aId];
    }

    moris_index
    Mesh_Checker::get_vertex_serial_index_from_id( moris_id aId, Serialized_Mesh_Data* aSerializedMesh )
    {
        MORIS_ERROR( aSerializedMesh->mSerialVertexMap.find( aId ) != aSerializedMesh->mSerialVertexMap.end(), "Vertex id already in the map" );
        return aSerializedMesh->mSerialVertexMap[aId];
    }

    moris_index
    Mesh_Checker::get_cell_proc_index_from_id( moris_id aId,
        moris_index                                     aProc,
        Serialized_Mesh_Data*                           aSerializedMesh )
    {
        MORIS_ERROR( aSerializedMesh->mCollectCellMaps( aProc ).find( aId ) != aSerializedMesh->mCollectCellMaps( aProc ).end(), "Vertex id not in the map" );
        return aSerializedMesh->mCollectCellMaps( aProc )[aId];
    }

    moris_index
    Mesh_Checker::get_cell_serial_index_from_id( moris_id aId, Serialized_Mesh_Data* aSerializedMesh )
    {
        MORIS_ERROR( aSerializedMesh->mSerialCellMap.find( aId ) != aSerializedMesh->mSerialCellMap.end(), "Vertex id already in the map" );
        return aSerializedMesh->mSerialCellMap[aId];
    }

    //--------------------------------------------------------------------------------


    void
    Mesh_Checker::serialize_mesh()
    {
        this->serialize_mesh_core();
    }
    //--------------------------------------------------------------------------------
    void
    Mesh_Checker::serialize_mesh_core()
    {
        // set spatial dimensions
        mSerializedIpMesh.mSpatialDim = mIpMesh->get_spatial_dim();
        mSerializedIgMesh.mSpatialDim = mIgMesh->get_spatial_dim();

        // serialize vertex
        this->serialize_vertices( mIpMesh, &mSerializedIpMesh );
        this->serialize_vertices( mIgMesh, &mSerializedIgMesh );

        // serialize cells
        this->serialize_cells( mIpMesh, &mSerializedIpMesh );
        this->serialize_cells( mIgMesh, &mSerializedIgMesh );

        // for the interpolation matrix serialize the data for vertex interp
        this->serialize_vertex_t_matrices( mIpMesh, &mSerializedIpMesh );
    }
    //--------------------------------------------------------------------------------
    void
    Mesh_Checker::serialize_vertices(
        Mesh*                 aMesh,
        Serialized_Mesh_Data* aSerialMesh )
    {
        size_t tNumNodes = aMesh->get_num_entities( (moris::EntityRank)EntityRank::NODE );

        aSerialMesh->mVertexIds.resize( tNumNodes, 1 );
        aSerialMesh->mVertexOwners.resize( tNumNodes, 1 );
        aSerialMesh->mVertexCoordinates.resize( tNumNodes, aMesh->get_spatial_dim() );

        for ( size_t i = 0; i < tNumNodes; i++ )
        {
            moris::mtk::Vertex* tVertex = &aMesh->get_mtk_vertex( i );
            // add cordinate to the data
            aSerialMesh->mVertexCoordinates.set_row( i, tVertex->get_coords() );

            // add id to the data
            aSerialMesh->mVertexIds( i ) = tVertex->get_id();

            // get the owner
            aSerialMesh->mVertexOwners( i ) = tVertex->get_owner();
        }
    }

    void
    Mesh_Checker::serialize_vertex_t_matrices(
        Interpolation_Mesh*   aMesh,
        Serialized_Mesh_Data* aSerializedMesh )
    {
        size_t tNumNodes = aMesh->get_num_entities( (moris::EntityRank)EntityRank::NODE );

        for ( size_t i = 0; i < tNumNodes; i++ )
        {
            Vertex& tVertex = aMesh->get_mtk_vertex( i );
            if ( tVertex.has_interpolation( mMeshIndex ) )
            {
                // vertex interpolation
                Vertex_Interpolation* tInterpolation = tVertex.get_interpolation( mMeshIndex );

                aSerializedMesh->mVertexTMatrixBasisOwners.push_back( tInterpolation->get_owners() );
                aSerializedMesh->mVertexTMatrixBasisIds.push_back( tInterpolation->get_ids() );
                aSerializedMesh->mVertexTMatrixWeights.push_back( *tInterpolation->get_weights() );
            }

            else
            {
                aSerializedMesh->mVertexTMatrixBasisOwners.push_back( Matrix< IndexMat >( 1, 1, MORIS_INDEX_MAX ) );
                aSerializedMesh->mVertexTMatrixBasisIds.push_back( Matrix< IndexMat >( 1, 1, MORIS_INDEX_MAX ) );
                aSerializedMesh->mVertexTMatrixWeights.push_back( Matrix< DDRMat >( 1, 1, 0.0 ) );
            }
        }
    }

    void
    Mesh_Checker::serialize_cells( Mesh* aMesh,
        Serialized_Mesh_Data*            aSerializedMesh )
    {
        size_t tNumCells = aMesh->get_num_entities( (moris::EntityRank)EntityRank::ELEMENT );

        aSerializedMesh->mCellIds.resize( tNumCells, 1 );
        aSerializedMesh->mCellOwners.resize( tNumCells, 1 );

        for ( size_t i = 0; i < tNumCells; i++ )
        {
            Cell const& tCell = aMesh->get_mtk_cell( i );
            // add id to the data
            aSerializedMesh->mCellIds( i ) = tCell.get_id();

            // get the owner
            aSerializedMesh->mCellOwners( i ) = tCell.get_owner();
        }
    }

    //--------------------------------------------------------------------------------
    void
    Mesh_Checker::gather_serialized_mesh( Serialized_Mesh_Data* aSerializedMesh )
    {

        moris_index tTag = 600;

        // gather vertex
        moris::all_gather_vector( aSerializedMesh->mVertexIds, aSerializedMesh->mCollectVertexIds, tTag, 1 );
        moris::all_gather_vector( aSerializedMesh->mVertexOwners, aSerializedMesh->mCollectVertexOwners, tTag + 1, 1 );
        moris::all_gather_vector( aSerializedMesh->mVertexCoordinates, aSerializedMesh->mCollectVertexCoords, tTag + 2, 1 );

        // gather cells
        moris::all_gather_vector( aSerializedMesh->mCellIds, aSerializedMesh->mCollectCellIds, tTag + 3, 1 );
        moris::all_gather_vector( aSerializedMesh->mCellOwners, aSerializedMesh->mCollectCellOwners, tTag + 4, 1 );
    }

    void
    Mesh_Checker::gather_serialized_ip_mesh( Serialized_Mesh_Data* aSerializedMesh )
    {
        moris_index tTag = 650;


        // get the basis weights
        Matrix< DDRMat >   tVertexTMatrixWeightsData;
        Matrix< IndexMat > tVertexTMatrixWeightsOffsets;
        this->cell_of_mats_to_flattened_mat( aSerializedMesh->mVertexTMatrixWeights, tVertexTMatrixWeightsData, tVertexTMatrixWeightsOffsets );

        moris::Cell< Matrix< DDRMat > >   tGatheredFlattenedWeightsData;
        moris::Cell< Matrix< IndexMat > > tGatheredFlattenedWeightsOffsets;
        moris::all_gather_vector( tVertexTMatrixWeightsData, tGatheredFlattenedWeightsData, tTag, 0 );
        moris::all_gather_vector( tVertexTMatrixWeightsOffsets, tGatheredFlattenedWeightsOffsets, tTag + 1, 0 );

        // get the basis ids
        Matrix< IndexMat > tVertexTMatrixIdsData;
        Matrix< IndexMat > tVertexTMatrixIdsOffsets;
        this->cell_of_mats_to_flattened_mat( aSerializedMesh->mVertexTMatrixBasisIds, tVertexTMatrixIdsData, tVertexTMatrixIdsOffsets );


        // gather the vertex t-matrices
        moris::Cell< Matrix< IndexMat > > tGatheredFlattenedBasisIdData;
        moris::Cell< Matrix< IndexMat > > tGatheredFlattenedBasisIdOffsets;
        moris::all_gather_vector( tVertexTMatrixIdsData, tGatheredFlattenedBasisIdData, tTag + 2, 0 );
        moris::all_gather_vector( tVertexTMatrixIdsOffsets, tGatheredFlattenedBasisIdOffsets, tTag + 3, 0 );

        //
        if ( par_rank() == 0 )
        {
            aSerializedMesh->mCollectVertexTMatrixWeights.resize( tGatheredFlattenedWeightsData.size() );
            aSerializedMesh->mCollectVertexTMatrixBasisIds.resize( tGatheredFlattenedBasisIdData.size() );

            // iterate through the data and construct cell of mats
            for ( moris::uint i = 0; i < tGatheredFlattenedWeightsData.size(); i++ )
            {
                this->flattened_mat_to_cell_of_mats( tGatheredFlattenedWeightsData( i ), tGatheredFlattenedWeightsOffsets( i ), aSerializedMesh->mCollectVertexTMatrixWeights( i ) );
                this->flattened_mat_to_cell_of_mats( tGatheredFlattenedBasisIdData( i ), tGatheredFlattenedBasisIdOffsets( i ), aSerializedMesh->mCollectVertexTMatrixBasisIds( i ) );
            }
        }


        //    moris::Cell<Matrix<DDRMat>> tTest;
        //    this->serial_mat_to_cell_of_mats(tVertexTMatrixWeightsData,tVertexTMatrixWeightsOffsets,tTest);
        //
        //    for(moris::uint i = 0; i < aSerializedMesh->mVertexTMatrixWeights.size(); i++)
        //    {
        //        MORIS_ASSERT(moris::norm(aSerializedMesh->mVertexTMatrixWeights(i)-tTest(i))<1e-8,"Mismatch");
        //    }
        //
        //    moris::print(tVertexTMatrixWeightsData,"tVertexTMatrixWeightsData");
        //    moris::print(tVertexTMatrixWeightsOffsets,"tVertexTMatrixWeightsOffsets");
    }

    //--------------------------------------------------------------------------------

    void
    Mesh_Checker::setup_vertex_maps( Serialized_Mesh_Data* aSerializedMesh )
    {
        moris_index tUniqueVertCount = 0;

        aSerializedMesh->mCollectVertexMaps.resize( aSerializedMesh->mCollectVertexIds.size() );

        // iterate through processors
        for ( moris::uint iP = 0; iP < aSerializedMesh->mCollectVertexIds.size(); iP++ )
        {
            Matrix< IndexMat > tUniqueCheck;
            unique( aSerializedMesh->mCollectVertexIds( iP ), tUniqueCheck );

            MORIS_ERROR( aSerializedMesh->mCollectVertexIds( iP ).numel() == tUniqueCheck.numel(), "Unique vertex issue from a proc" );
            // iterate through vertices on proc
            for ( moris::uint iV = 0; iV < aSerializedMesh->mCollectVertexIds( iP ).numel(); iV++ )
            {
                moris_id tVertexId = aSerializedMesh->mCollectVertexIds( iP )( iV );

                MORIS_ERROR( aSerializedMesh->mCollectVertexMaps( iP ).find( tVertexId ) == aSerializedMesh->mCollectVertexMaps( iP ).end(), "Vertex id already in the map" );
                aSerializedMesh->mCollectVertexMaps( iP )[tVertexId] = iV;

                if ( aSerializedMesh->mSerialVertexMap.find( tVertexId ) == aSerializedMesh->mSerialVertexMap.end() )
                {
                    aSerializedMesh->mSerialVertexMap[tVertexId] = tUniqueVertCount;
                    aSerializedMesh->mVertexSerialIndexToId.push_back( tVertexId );
                    tUniqueVertCount++;
                }
            }
        }
    }

    void
    Mesh_Checker::setup_cell_maps( Serialized_Mesh_Data* aSerializedMesh )
    {
        moris_index tUniqueCellCount = 0;

        aSerializedMesh->mCollectCellMaps.resize( aSerializedMesh->mCollectCellIds.size() );

        // iterate through processors
        for ( moris::uint iP = 0; iP < aSerializedMesh->mCollectCellIds.size(); iP++ )
        {
            Matrix< IndexMat > tUniqueCheck;
            unique( aSerializedMesh->mCollectCellIds( iP ), tUniqueCheck );
            MORIS_ERROR( aSerializedMesh->mCollectCellIds( iP ).numel() == tUniqueCheck.numel(), "Unique vertex issue from a proc" );

            // iterate through vertices on proc
            for ( moris::uint iC = 0; iC < aSerializedMesh->mCollectCellIds( iP ).numel(); iC++ )
            {
                moris_id tCellId = aSerializedMesh->mCollectCellIds( iP )( iC );

                MORIS_ERROR( aSerializedMesh->mCollectCellMaps( iP ).find( tCellId ) == aSerializedMesh->mCollectCellMaps( iP ).end(), "Vertex id already in the map" );
                aSerializedMesh->mCollectCellMaps( iP )[tCellId] = iC;

                if ( aSerializedMesh->mSerialCellMap.find( tCellId ) == aSerializedMesh->mSerialCellMap.end() )
                {
                    aSerializedMesh->mSerialCellMap[tCellId] = tUniqueCellCount;
                    aSerializedMesh->mCellSerialIndexToId.push_back( tCellId );
                    tUniqueCellCount++;
                }
            }
        }
    }

    //--------------------------------------------------------------------------------

    std::string
    Mesh_Checker::bool_to_string( bool aBool )
    {
        if ( aBool )
        {
            return "Pass";
        }
        else
        {
            return "Fail";
        }
    }
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

}// namespace mtk
}// namespace moris
