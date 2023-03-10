/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Cut_Integration_Mesh.cpp
 *
 */

#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "typedefs.hpp"

#include "fn_stringify_matrix.hpp"

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_XTK_Model.hpp"
// #include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_MPI_Tools.hpp"
#include "cl_XTK_Cell_No_CM.hpp"
#include "fn_determine_cell_topology.hpp"
using namespace moris;

namespace xtk
{
    // ----------------------------------------------------------------------------------

    IG_Cell_Group::IG_Cell_Group(
            moris_index aNumCellsInGroup )
            : mIgCellGroup( aNumCellsInGroup, nullptr )
    {
    }

    // ----------------------------------------------------------------------------------

    IG_Cell_Group::IG_Cell_Group()
    {
    }

    // ----------------------------------------------------------------------------------

    void
    IG_Cell_Group::add_Cell( moris::mtk::Cell* aCell )
    {
        moris_index tNewCellOrdinal = (moris_index)mIgCellGroup.size();

        if ( mIgCellIndexToCellOrdinal.find( aCell->get_index() ) != mIgCellIndexToCellOrdinal.end() )
        {
            std::cout << "New cell index = " << aCell->get_index() << " | mIgCellIndexToCellOrdinal.find( aCell->get_index() ) = " << mIgCellIndexToCellOrdinal.find( aCell->get_index() )->first << std::endl;
        }

        MORIS_ASSERT( mIgCellIndexToCellOrdinal.find( aCell->get_index() ) == mIgCellIndexToCellOrdinal.end(), "Duplicate vertex in group" );
        mIgCellGroup.push_back( aCell );
        mIgCellIndexToCellOrdinal[ aCell->get_index() ] = tNewCellOrdinal;
    }

    // ----------------------------------------------------------------------------------

    moris_index
    IG_Cell_Group::get_cell_group_ordinal( moris_index aCell )
    {
        auto tIter = mIgCellIndexToCellOrdinal.find( aCell );
        MORIS_ERROR( tIter != mIgCellIndexToCellOrdinal.end(), "Provided cell not in group" );
        return tIter->second;
    }

    // ----------------------------------------------------------------------------------

    void
    IG_Cell_Group::remove_cell( moris_index aCell )
    {

        moris_index tOrdinal = this->get_cell_group_ordinal( aCell );

        mIgCellGroup.erase( tOrdinal );
        mIgCellIndexToCellOrdinal.erase( aCell );

        // rewrite index to ordinal map
        std::unordered_map< moris_index, moris_index > tTempMap;

        for ( std::pair< moris_index, moris_index > el : mIgCellIndexToCellOrdinal )
        {
            moris_index key = el.first;
            moris_index val = el.second;
            if ( el.first >= aCell )
            {
                key = el.first - 1;
            }
            if ( el.second >= tOrdinal )
            {
                val = el.second - 1;
            }

            tTempMap[ key ] = val;
        }
        mIgCellIndexToCellOrdinal = tTempMap;
    }

    // ----------------------------------------------------------------------------------

    bool
    IG_Cell_Group::cell_is_in_group( moris_index aCell )
    {
        return mIgCellIndexToCellOrdinal.find( aCell ) != mIgCellIndexToCellOrdinal.end();
    }

    // ----------------------------------------------------------------------------------

    void
    IG_Cell_Group::shift_indices( moris_index aCell )
    {
        std::unordered_map< moris_index, moris_index > tTempMap;

        for ( std::pair< moris_index, moris_index > el : mIgCellIndexToCellOrdinal )
        {
            moris_index key = el.first;
            moris_index val = el.second;
            if ( el.first >= aCell )
            {
                key = el.first - 1;
            }

            tTempMap[ key ] = val;
        }
        mIgCellIndexToCellOrdinal = tTempMap;
    }

    // ----------------------------------------------------------------------------------

    IG_Cell_Side_Group::IG_Cell_Side_Group( moris_index aEstimatedNumCells )
    {
        mIgCells.reserve( aEstimatedNumCells );
        mIgCellSideOrdinals.reserve( aEstimatedNumCells );
    }

    // ----------------------------------------------------------------------------------

    IG_Cell_Double_Side_Group::IG_Cell_Double_Side_Group( moris_index aEstimatedNumCells )
    {
        mLeaderIgCells.reserve( aEstimatedNumCells );
        mLeaderIgCellSideOrdinals.reserve( aEstimatedNumCells );
        mFollowerIgCells.reserve( aEstimatedNumCells );
        mFollowerIgCellSideOrdinals.reserve( aEstimatedNumCells );
    }

    // ----------------------------------------------------------------------------------

    IG_Vertex_Group::IG_Vertex_Group( moris_index aNumVerticesInGroup )
            : mIgVertexGroup( 0, nullptr )
    {
        mIgVertexGroup.reserve( aNumVerticesInGroup );
    }

    // ----------------------------------------------------------------------------------

    std::size_t
    IG_Vertex_Group::size()
    {
        return mIgVertexGroup.size();
    }

    // ----------------------------------------------------------------------------------

    void
    IG_Vertex_Group::reserve( std::size_t aReserveSize )
    {
        mIgVertexGroup.reserve( aReserveSize );
        mIgVertexLocalCoords.reserve( aReserveSize );
    }

    // ----------------------------------------------------------------------------------

    void
    IG_Vertex_Group::add_vertex(
            moris::mtk::Vertex const *          aVertex,
            std::shared_ptr< Matrix< DDRMat > > aVertexLocalCoord )
    {
        moris_index tNewVertexOrdinal = (moris_index)mIgVertexGroup.size();

        if ( mIgVertexIndexToVertexOrdinal.find( aVertex->get_index() ) != mIgVertexIndexToVertexOrdinal.end() )
        {
            std::cout << "New vertex index = " << aVertex->get_index() << " | mIgVertexIndexToVertexOrdinal.find( aVertex->get_index() ) = " << mIgVertexIndexToVertexOrdinal.find( aVertex->get_index() )->first << std::endl;
            this->print();
        }

        MORIS_ASSERT( mIgVertexIndexToVertexOrdinal.find( aVertex->get_index() ) == mIgVertexIndexToVertexOrdinal.end(), "Duplicate vertex in group" );
        MORIS_ASSERT( mIgVertexGroup.size() == mIgVertexLocalCoords.size(), "Size issue" );
        mIgVertexGroup.push_back( aVertex );
        mIgVertexLocalCoords.push_back( aVertexLocalCoord );
        mIgVertexIndexToVertexOrdinal[ aVertex->get_index() ] = tNewVertexOrdinal;
    }

    // ----------------------------------------------------------------------------------

    void
    IG_Vertex_Group::add_vertex_local_coord_pointers()
    {
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Vertex const *
    IG_Vertex_Group::get_vertex( moris_index aGroupVertexOrdinal )
    {
        return mIgVertexGroup( aGroupVertexOrdinal );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    IG_Vertex_Group::get_vertex_group_ordinal( moris_index aVertex )
    {
        auto tIter = mIgVertexIndexToVertexOrdinal.find( aVertex );
        MORIS_ERROR( tIter != mIgVertexIndexToVertexOrdinal.end(), "Provided vertex not in group" );
        return tIter->second;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Matrix< DDRMat > >
    IG_Vertex_Group::get_vertex_local_coords( moris_index aVertex )
    {
        return mIgVertexLocalCoords( this->get_vertex_group_ordinal( aVertex ) );
    }

    // ----------------------------------------------------------------------------------

    bool
    IG_Vertex_Group::vertex_is_in_group( moris_index aVertex )
    {
        return mIgVertexIndexToVertexOrdinal.find( aVertex ) != mIgVertexIndexToVertexOrdinal.end();
    }

    // ----------------------------------------------------------------------------------

    uint
    IG_Vertex_Group::get_vertex_local_coords_dim() const
    {
        if ( mIgVertexLocalCoords.size() > 0 )
        {
            return mIgVertexLocalCoords( 0 )->n_cols();
        }
        else
        {
            return 0;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    IG_Vertex_Group::print()
    {
        // iterate through vertices
        for ( uint iV = 0; iV < this->size(); iV++ )
        {
            std::cout << std::setw( 8 ) << iV << " : ";
            std::cout << "Vertex Id: " << std::setw( 8 ) << mIgVertexGroup( iV )->get_id();
            std::cout << " | Vertex Index: " << std::setw( 8 ) << mIgVertexGroup( iV )->get_index();

            Matrix< DDRMat > tVertexCoords = mIgVertexGroup( iV )->get_coords();

            std::cout << " | Coords:";
            for ( uint iSpatial = 0; iSpatial < tVertexCoords.numel(); iSpatial++ )
            {
                std::cout << " " << std::setw( 16 ) << std::scientific << tVertexCoords( iSpatial );
            }
            std::cout << std::endl;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    IG_Vertex_Group::remove_vertex( moris_index aVertex )
    {
        moris_index tOrdinal = this->get_vertex_group_ordinal( aVertex );

        mIgVertexGroup.erase( tOrdinal );
        mIgVertexLocalCoords.erase( tOrdinal );
        mIgVertexIndexToVertexOrdinal.erase( aVertex );

        // rewrite index to ordinal map
        std::unordered_map< moris_index, moris_index > tTempMap;

        for ( std::pair< moris_index, moris_index > el : mIgVertexIndexToVertexOrdinal )
        {
            moris_index key = el.first;
            moris_index val = el.second;
            if ( el.first >= aVertex )
            {
                key = el.first - 1;
            }
            if ( el.second >= tOrdinal )
            {
                val = el.second - 1;
            }

            tTempMap[ key ] = val;
        }

        mIgVertexIndexToVertexOrdinal = tTempMap;
    }

    // ----------------------------------------------------------------------------------

    void
    IG_Vertex_Group::shift_indices( moris_index aVertex )
    {
        std::unordered_map< moris_index, moris_index > tTempMap;

        for ( std::pair< moris_index, moris_index > el : mIgVertexIndexToVertexOrdinal )
        {
            moris_index key = el.first;
            moris_index val = el.second;
            if ( el.first >= aVertex )
            {
                key = el.first - 1;
            }

            tTempMap[ key ] = val;
        }

        mIgVertexIndexToVertexOrdinal = tTempMap;
    }

    // ----------------------------------------------------------------------------------

    Cut_Integration_Mesh::Cut_Integration_Mesh(
            moris::mtk::Mesh* aBackgroundMesh,
            xtk::Model*       aXTKModel )
            : mBackgroundMesh( aBackgroundMesh )
            , mXTKModel( aXTKModel )
    {
        //
        mSpatialDimension = aBackgroundMesh->get_spatial_dim();

        // setup the integration cells
        uint tNumBackgroundCells    = mBackgroundMesh->get_num_elems();
        uint tNumBackgroundVertices = mBackgroundMesh->get_num_nodes();

        // cell setup
        mIntegrationCells.resize( tNumBackgroundCells, nullptr );
        mIntegrationCellIndexToId.resize( tNumBackgroundCells, MORIS_INDEX_MAX );
        mIntegrationCellToCellGroupIndex.resize( tNumBackgroundCells, 0 );
        mParentCellCellGroupIndex.resize( tNumBackgroundCells, MORIS_INDEX_MAX );

        mIntegrationCellBulkPhase.resize( tNumBackgroundCells, MORIS_INDEX_MAX );

        for ( uint iCell = 0; iCell < tNumBackgroundCells; iCell++ )
        {
            mIntegrationCells( iCell ) = &mBackgroundMesh->get_mtk_cell( (moris_index)iCell );

            // verify that we are not doubling up vertices in the id map
            MORIS_ERROR( mIntegrationCellIdToIndexMap.find( mIntegrationCells( iCell )->get_id() ) == mIntegrationCellIdToIndexMap.end(), "Provided Cell Id is already in the integration vertex map: Vertex Id =%uon process %u", mIntegrationCells( iCell )->get_id(), par_rank() );

            // add the vertex to id to index map
            mIntegrationCellIdToIndexMap[ mIntegrationCells( iCell )->get_id() ] = (moris_index)iCell;

            mParentCellCellGroupIndex( iCell ) = iCell;
        }

        // vertex setup
        mIntegrationVertices.resize( tNumBackgroundVertices, nullptr );
        mIntegrationVertexIndexToId.resize( tNumBackgroundVertices, MORIS_INDEX_MAX );
        mVertexCoordinates.resize( tNumBackgroundVertices, nullptr );
        mIgVertexParentEntityIndex.resize( tNumBackgroundVertices );
        mIgVertexParentEntityRank.resize( tNumBackgroundVertices, 0 );

        for ( uint iV = 0; iV < tNumBackgroundVertices; iV++ )
        {
            // get a vertex pointer into our data
            mIntegrationVertices( iV ) = &mBackgroundMesh->get_mtk_vertex( (moris_index)iV );

            // this vertex's parent is itself
            mIgVertexParentEntityIndex( iV ) = iV;

            // check the vertex index lines up correctly
            MORIS_ERROR( mIntegrationVertices( iV )->get_index() == (moris_index)iV, "Vertex index mismatch" );

            // add the vertex to the local to global vertex map
            mIntegrationVertexIndexToId( mIntegrationVertices( iV )->get_index() ) = mIntegrationVertices( iV )->get_id();

            // store the coordinate
            mVertexCoordinates( mIntegrationVertices( iV )->get_index() ) = std::make_shared< moris::Matrix< DDRMat > >( mIntegrationVertices( iV )->get_coords() );

            // verify that we are not doubling up vertices in the id map
            MORIS_ERROR( mIntegrationVertexIdToIndexMap.find( mIntegrationVertices( iV )->get_id() ) == mIntegrationVertexIdToIndexMap.end(), "Provided Vertex Id is already in the integration vertex map: Vertex Id =%uon process %u", mIntegrationVertices( iV )->get_id(), par_rank() );

            // add the vertex to id to index map
            mIntegrationVertexIdToIndexMap[ mIntegrationVertices( iV )->get_id() ] = (moris_index)iV;
        }

        this->setup_comm_map();

        // max vertex and cell ids (to allocate new ones later)
        mGlobalMaxVertexId = mBackgroundMesh->get_max_entity_id( EntityRank::NODE ) + 1;
        mGlobalMaxCellId   = mBackgroundMesh->get_max_entity_id( EntityRank::ELEMENT ) + 1;

        mFirstControlledCellIndex   = mIntegrationCells.size();
        mFirstControlledVertexIndex = mIntegrationVertices.size();

        // setup base cell blocks
        this->create_base_cell_blocks();
        if ( mXTKModel->mDiagnostics )
        {
            std::string tCellDiagFile    = mXTKModel->get_diagnostic_file_name( std::string( "Cells_On_Init" ) );
            std::string tVertDiagFile    = mXTKModel->get_diagnostic_file_name( std::string( "Vertex_On_Init" ) );
            std::string tGENVertDiagFile = mXTKModel->get_diagnostic_file_name( std::string( "Vertex_On_Init_GEN" ) );
            std::string tGroupDiagFile   = mXTKModel->get_diagnostic_file_name( std::string( "Groups_On_Init" ) );
            this->print_cells( false, tCellDiagFile );
            this->print_vertices( false, tVertDiagFile );
            this->print_groupings( tGroupDiagFile );
            mXTKModel->get_geom_engine()->print_gen_vertices( tGENVertDiagFile, this );
        }
    }

    // ----------------------------------------------------------------------------------

    Cut_Integration_Mesh::~Cut_Integration_Mesh()
    {
        // delete the objects holding B-spline information
        this->delete_Bspline_mesh_info();
    }

    //-------------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::delete_Bspline_mesh_info()
    {
        // delete the objects themselves
        for ( auto iBspMesh : mBsplineMeshInfos )
        {
            delete iBspMesh;
        }

        // clear memory for cell
        mBsplineMeshInfos.clear();
    }

    // ----------------------------------------------------------------------------------

    uint
    Cut_Integration_Mesh::get_spatial_dim() const
    {
        return mSpatialDimension;
    }

    // ----------------------------------------------------------------------------------

    MeshType
    Cut_Integration_Mesh::get_mesh_type() const
    {
        return MeshType::XTK;
    }

    // ----------------------------------------------------------------------------------

    uint
    Cut_Integration_Mesh::get_num_sets() const
    {
        return 0;
    }

    // ----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cut_Integration_Mesh::get_node_coordinate( moris_index aNodeIndex ) const
    {
        return ( *mVertexCoordinates( aNodeIndex ) );
    }

    // ----------------------------------------------------------------------------------

    uint
    Cut_Integration_Mesh::get_node_owner( moris_index aNodeIndex ) const
    {
        return mIntegrationVertices( aNodeIndex )->get_owner();
    }

    // ----------------------------------------------------------------------------------

    uint
    Cut_Integration_Mesh::get_element_owner( moris_index aElementIndex ) const
    {
        return mIntegrationCells( aElementIndex )->get_owner();
    }

    // ----------------------------------------------------------------------------------

    Matrix< IdMat >
    Cut_Integration_Mesh::get_communication_table() const
    {
        return mCommunicationMap;
    }

    // ----------------------------------------------------------------------------------

    std::map< moris_id, moris_index >
    Cut_Integration_Mesh::get_communication_map()
    {
        // construct the index map if it hasn't already, or the communication table was updated
        if( !mCommMapHasBeenConstructed )
        {
            mCommunicationIndexMap.clear();
            for ( uint iProc = 0; iProc < mCommunicationMap.numel(); iProc++ )
            {
                mCommunicationIndexMap[ mCommunicationMap( iProc ) ] = (moris_index)iProc;
            }
        }

        // mark the index map as constructed
        mCommMapHasBeenConstructed = true;

        // return the index map
        return mCommunicationIndexMap;
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::add_proc_to_comm_table( moris_index aProcRank )
    {
        moris_index tIndex = mCommunicationMap.numel();

        for ( uint i = 0; i < mCommunicationMap.numel(); i++ )
        {
            MORIS_ERROR( mCommunicationMap( i ) != aProcRank, "Processor rank already in communication table" );
        }

        mCommunicationMap.resize( 1, mCommunicationMap.numel() + 1 );
        mCommunicationMap( tIndex ) = aProcRank;

        // mark the index map as needing to be re-constructed
        mCommMapHasBeenConstructed = false;
    }

    // ----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Cut_Integration_Mesh::get_element_indices_in_block_set( uint aSetIndex )
    {
        // get cells in set
        std::shared_ptr< IG_Cell_Group > tCellsInBlock = mBlockSetCellGroup( aSetIndex );

        if ( tCellsInBlock != nullptr )
        {
            Matrix< IndexMat > tCellIndices( 1, tCellsInBlock->mIgCellGroup.size() );

            for ( uint iCell = 0; iCell < tCellsInBlock->mIgCellGroup.size(); iCell++ )
            {
                tCellIndices( iCell ) = tCellsInBlock->mIgCellGroup( iCell )->get_index();
            }

            return tCellIndices;
        }

        return Matrix< IndexMat >( 0, 0 );
    }

    // ----------------------------------------------------------------------------------

    enum CellTopology
    Cut_Integration_Mesh::get_blockset_topology( const std::string& aSetName )
    {
        // get index
        moris_index tBlockIndex = this->get_block_set_index( aSetName );

        return mBlockCellTopo( tBlockIndex );
    }

    // ----------------------------------------------------------------------------------

    enum CellShape
    Cut_Integration_Mesh::get_IG_blockset_shape( const std::string& aSetName )
    {
        return CellShape::INVALID;
    }

    // ----------------------------------------------------------------------------------

    enum CellShape
    Cut_Integration_Mesh::get_IP_blockset_shape( const std::string& aSetName )
    {
        return CellShape::INVALID;
    }

    // ----------------------------------------------------------------------------------

    moris_id
    Cut_Integration_Mesh::get_glb_entity_id_from_entity_loc_index(
            moris_index       aEntityIndex,
            enum EntityRank   aEntityRank,
            const moris_index aDiscretizationIndex ) const
    {
        MORIS_ERROR( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT, "Only supported for nodes and cells" );
        if ( aEntityRank == EntityRank::NODE )
        {
            return mIntegrationVertices( aEntityIndex )->get_index();
        }
        else if ( aEntityRank == EntityRank::ELEMENT )
        {
            return mIntegrationCells( aEntityIndex )->get_index();
        }

        else
        {
            return 0;
        }
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_loc_entity_ind_from_entity_glb_id(
            moris_id        aEntityId,
            enum EntityRank aEntityRank ) const
    {
        // warning element map is set up after integration mesh has been constructed
        MORIS_ERROR( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT, "Only a node map and element map is implemented in XTK" );

        if ( aEntityRank == EntityRank::NODE )
        {
            auto tIter = mIntegrationVertexIdToIndexMap.find( aEntityId );

            MORIS_ERROR( tIter != mIntegrationVertexIdToIndexMap.end(),
                    "Cut_Integration_Mesh::get_loc_entity_ind_from_entity_glb_id() - "
                    "Provided Entity Id is not in the map, Has the map been initialized?: aEntityId =%u EntityRank = %u on process %u",
                    aEntityId,
                    (uint)aEntityRank,
                    par_rank() );
            return tIter->second;
        }
        else if ( aEntityRank == EntityRank::ELEMENT )
        {
            auto tIter = mIntegrationCellIdToIndexMap.find( aEntityId );

            MORIS_ERROR( tIter != mIntegrationCellIdToIndexMap.end(),
                    "Cut_Integration_Mesh::get_loc_entity_ind_from_entity_glb_id() - "
                    "Provided Entity Id is not in the map, Has the map been initialized?: aEntityId =%u EntityRank = %u on process %u",
                    aEntityId,
                    (uint)aEntityRank,
                    par_rank() );
            return tIter->second;
        }

        else
        {
            return 0;
        }
    }

    // ----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Cut_Integration_Mesh::get_entity_connected_to_entity_loc_inds(
            moris_index       aEntityIndex,
            enum EntityRank   aInputEntityRank,
            enum EntityRank   aOutputEntityRank,
            const moris_index aDiscretizationIndex ) const
    {
        MORIS_ERROR( aInputEntityRank == EntityRank::ELEMENT && aOutputEntityRank == EntityRank::NODE,
                "Only support element to node connectivity" );

        return this->get_mtk_cell( aEntityIndex ).get_vertex_inds();
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< std::string >
    Cut_Integration_Mesh::get_set_names( enum EntityRank aSetEntityRank ) const
    {
        switch ( aSetEntityRank )
        {
            case EntityRank::NODE:
            {
                return {};
                break;
            }
            case EntityRank::EDGE:
            {
                MORIS_ASSERT( this->get_facet_rank() == EntityRank::EDGE, "side sets are defined on edges in 2d" );
                return mSideSetLabels;
                break;
            }
            case EntityRank::FACE:
            {
                MORIS_ASSERT( this->get_facet_rank() == EntityRank::FACE, "side sets are defined on faces in 3d" );
                return mSideSetLabels;
                break;
            }
            case EntityRank::ELEMENT:
            {
                return mBlockSetNames;
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Currently only supporting block, node and side sets in XTK enriched integration meshes" );
            }
                return moris::Cell< std::string >( 0 );
                break;
        }
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_block_set_index( std::string aBlockSetLabel ) const
    {
        auto tIter = mBlockSetLabelToOrd.find( aBlockSetLabel );

        MORIS_ERROR( tIter != mBlockSetLabelToOrd.end(), "block set set label not found" );

        return tIter->second;
    }

    // ----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Cut_Integration_Mesh::get_block_entity_loc_inds( std::string aSetName ) const
    {
        // get index
        moris_index tBlockIndex = this->get_block_set_index( aSetName );

        // get cells in set
        std::shared_ptr< IG_Cell_Group > tCellsInBlock = mBlockSetCellGroup( tBlockIndex );

        Matrix< IndexMat > tCellIndices( 1, tCellsInBlock->mIgCellGroup.size() );

        for ( uint iCell = 0; iCell < tCellsInBlock->mIgCellGroup.size(); iCell++ )
        {
            tCellIndices( iCell ) = tCellsInBlock->mIgCellGroup( iCell )->get_index();
        }

        return tCellIndices;
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_side_set_index( std::string aSideSetLabel ) const
    {
        auto tIter = mSideSideSetLabelToOrd.find( aSideSetLabel );

        MORIS_ERROR( tIter != mSideSideSetLabelToOrd.end(), "side side set label not found" );

        return tIter->second;
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::get_sideset_elems_loc_inds_and_ords(
            const std::string&  aSetName,
            Matrix< IndexMat >& aElemIndices,
            Matrix< IndexMat >& aSidesetOrdinals ) const
    {
        // get the index
        moris_index tSideSetIndex = this->get_side_set_index( aSetName );

        std::shared_ptr< IG_Cell_Side_Group > tCellSidesInSet = mSideSetCellSides( tSideSetIndex );

        // iterate through side clusters and count number of sides in set

        if ( tCellSidesInSet != nullptr )
        {
            uint tNumSides = tCellSidesInSet->mIgCells.size();

            // size outputs
            aElemIndices.resize( 1, tNumSides );
            aSidesetOrdinals.resize( 1, tNumSides );

            for ( uint iSide = 0; iSide < tNumSides; iSide++ )
            {
                aElemIndices( iSide )     = tCellSidesInSet->mIgCells( iSide )->get_index();
                aSidesetOrdinals( iSide ) = tCellSidesInSet->mIgCellSideOrdinals( iSide );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Cut_Integration_Mesh::get_set_entity_loc_inds(
            enum EntityRank aSetEntityRank,
            std::string     aSetName ) const
    {
        switch ( aSetEntityRank )
        {
            case EntityRank::NODE:
            {
                return Matrix< IndexMat >( 0, 0 );
                break;
            }
            case EntityRank::EDGE:
            {
                // MORIS_ASSERT(this->get_facet_rank() == EntityRank::EDGE,"side sets are defined on edges in 2d");
                return Matrix< IndexMat >( 0, 0 );
                break;
            }
            case EntityRank::FACE:
            {
                // MORIS_ASSERT(this->get_facet_rank() == EntityRank::FACE,"side sets are defined on faces in 3d");
                return Matrix< IndexMat >( 0, 0 );
                break;
            }
            case EntityRank::ELEMENT:
            {
                return this->get_block_entity_loc_inds( aSetName );
                return Matrix< IndexMat >( 0, 0 );
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Currently only supporting block, node and side sets in XTK enriched integration meshes" );
                return Matrix< IndexMat >( 0, 0 );
                break;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< std::shared_ptr< Matrix< DDRMat > > >*
    Cut_Integration_Mesh::get_all_vertex_coordinates_loc_inds()
    {
        return &mVertexCoordinates;
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< std::shared_ptr< IG_Cell_Group > >&
    Cut_Integration_Mesh::get_all_cell_groups()
    {
        return mIntegrationCellGroups;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Child_Mesh_Experimental >
    Cut_Integration_Mesh::get_child_mesh( moris_index aChildMeshIndex )
    {
        MORIS_ERROR( mChildMeshes.size() > (uint)aChildMeshIndex, "Child mesh index out of bounds" );
        return mChildMeshes( aChildMeshIndex );
    }

    // ----------------------------------------------------------------------------

    uint
    Cut_Integration_Mesh::get_num_child_meshes() const
    {
        return mChildMeshes.size();
    }

    // ----------------------------------------------------------------------------

    mtk::Vertex&
    Cut_Integration_Mesh::get_mtk_vertex( moris_index aVertexIndex )
    {
        MORIS_ASSERT( aVertexIndex < (moris_index)mIntegrationVertices.size(),
                "Vertex index out of bounds" );

        return *mIntegrationVertices( aVertexIndex );
    }

    // ----------------------------------------------------------------------------

    mtk::Vertex const &
    Cut_Integration_Mesh::get_mtk_vertex( moris_index aVertexIndex ) const
    {
        MORIS_ASSERT( aVertexIndex < (moris_index)mIntegrationVertices.size(),
                "Vertex index out of bounds" );

        return *mIntegrationVertices( aVertexIndex );
    }

    // ----------------------------------------------------------------------------------

    std::unordered_map< moris_id, moris_index >
    Cut_Integration_Mesh::get_vertex_glb_id_to_loc_vertex_ind_map() const
    {
        return mIntegrationVertexIdToIndexMap;
    }

    // ----------------------------------------------------------------------------------

    bool
    Cut_Integration_Mesh::vertex_exists( moris_index tId ) const
    {
        return mIntegrationVertexIdToIndexMap.find( tId ) != mIntegrationVertexIdToIndexMap.end();
    }

    // ----------------------------------------------------------------------------------

    mtk::Cell const &
    Cut_Integration_Mesh::get_mtk_cell( moris_index aElementIndex ) const
    {
        return *mIntegrationCells( aElementIndex );
    }
    mtk::Cell&
    Cut_Integration_Mesh::get_mtk_cell( moris_index aElementIndex )
    {
        return *mIntegrationCells( aElementIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Vertex*
    Cut_Integration_Mesh::get_mtk_vertex_pointer( moris_index aVertexIndex )
    {
        return mIntegrationVertices( aVertexIndex );
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< IG_Vertex_Group >
    Cut_Integration_Mesh::get_vertex_group( moris_index aVertexGroupIndex )
    {
        return mIntegrationVertexGroups( aVertexGroupIndex );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_parent_cell_group_index( moris_index aParentCellIndex )
    {
        return mParentCellCellGroupIndex( aParentCellIndex );
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::replace_controlled_ig_cell(
            moris_index                              aCellIndex,
            moris_id                                 aCellId,
            std::shared_ptr< moris::mtk::Cell_Info > aCellInfo,
            moris::Cell< moris::mtk::Vertex* >&      aVertexPointers )
    {
        MORIS_ERROR( aCellIndex >= mFirstControlledCellIndex, "Cannot set integration cell that I do not control." );

        // controlled index
        moris_index tIndexInControlledCells = aCellIndex - mFirstControlledCellIndex;

        mControlledIgCells( tIndexInControlledCells )->set_id( aCellId );
        mControlledIgCells( tIndexInControlledCells )->set_mtk_cell_info( aCellInfo );
        mControlledIgCells( tIndexInControlledCells )->set_vertex_pointers( aVertexPointers );
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::set_integration_cell(
            moris_index                            aCellIndex,
            std::shared_ptr< xtk::Cell_XTK_No_CM > aNewCell )
    {
        MORIS_ERROR( aCellIndex >= mFirstControlledCellIndex, "Cannot set integration cell that I do not control." );

        // controlled index
        moris_index tIndexInControlledCells = aCellIndex - mFirstControlledCellIndex;

        mControlledIgCells( tIndexInControlledCells ) = aNewCell;
        mIntegrationCells( aCellIndex )               = mControlledIgCells( tIndexInControlledCells ).get();
    }

    void
    Cut_Integration_Mesh::add_integration_cell(
            moris_index                            aCellIndex,
            std::shared_ptr< xtk::Cell_XTK_No_CM > aNewCell )
    {
        MORIS_ERROR( (uint)aCellIndex >= mIntegrationCells.size(), "Index mismatch between adding cell and current data." );

        mControlledIgCells.push_back( aNewCell );
        mIntegrationCells.push_back( aNewCell.get() );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_integration_cell_controlled_index(
            moris_index aCellIndex )
    {
        MORIS_ERROR( aCellIndex >= mFirstControlledCellIndex, "Cell index provided I do not control" );

        return aCellIndex - mFirstControlledCellIndex;
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_integration_vertex_controlled_index(
            moris_index aVertIndex )
    {
        MORIS_ERROR( aVertIndex >= mFirstControlledVertexIndex,
                "Cut_Integration_Mesh::get_integration_vertex_controlled_index() - Vertex index provided I do not control" );

        if ( aVertIndex >= mFirstControlledVertexIndex )
        {
            return aVertIndex - mFirstControlledVertexIndex;
        }
        else
        {
            return -1;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::add_cell_to_cell_group(
            moris_index aCellIndex,
            moris_index aCellGroupIndex )
    {
        MORIS_ERROR( aCellGroupIndex < (moris_index)mIntegrationCellGroups.size(),
                "Cut_Integration_Mesh::add_cell_to_cell_group() - Child Mesh Index out of bounds." );
        MORIS_ERROR( aCellIndex < (moris_index)mIntegrationCells.size(),
                "Cut_Integration_Mesh::add_cell_to_cell_group() - Cell Index out of bounds." );

        // mIntegrationCellGroups( aCellGroupIndex )->mIgCellGroup.push_back( mIntegrationCells( aCellIndex ) );
        mIntegrationCellGroups( aCellGroupIndex )->add_Cell( mIntegrationCells( aCellIndex ) );
        mIntegrationCellToCellGroupIndex( aCellIndex ).push_back( aCellGroupIndex );
    }

    // ----------------------------------------------------------------------------------

    moris_id
    Cut_Integration_Mesh::allocate_entity_ids(
            moris::size_t   aNumIdsToAllocate,
            enum EntityRank aEntityRank )
    {
        MORIS_ERROR( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT, "Only node and element ids can be allocated with xtk." );
        int tProcRank = moris::par_rank();
        int tProcSize = moris::par_size();

        // size_t is defined as uint here because of aNumRequested
        // Initialize gathered information outputs (information which will be scattered across processors)
        moris::Cell< moris_id > aGatheredInfo;
        moris::Cell< moris_id > tFirstId( 1 );
        moris::Cell< moris_id > tNumIdsRequested( 1 );

        tNumIdsRequested( 0 ) = (moris_id)aNumIdsToAllocate;

        moris::gather( tNumIdsRequested, aGatheredInfo );

        moris::Cell< moris_id > tProcFirstID( tProcSize );

        if ( tProcRank == 0 )
        {
            // Loop over entities print the number of entities requested by each processor
            for ( int iProc = 0; iProc < tProcSize; ++iProc )
            {
                if ( aEntityRank == EntityRank::NODE )
                {
                    // Give each processor their desired amount of IDs
                    tProcFirstID( iProc ) = mGlobalMaxVertexId;

                    // Increment the first available node ID
                    mGlobalMaxVertexId = mGlobalMaxVertexId + aGatheredInfo( iProc );
                }

                else if ( aEntityRank == EntityRank::ELEMENT )
                {
                    // Give each processor their desired amount of IDs
                    tProcFirstID( iProc ) = mGlobalMaxCellId;

                    // Increment the first available node ID
                    mGlobalMaxCellId = mGlobalMaxCellId + aGatheredInfo( iProc );
                }
            }
        }

        moris::scatter( tProcFirstID, tFirstId );

        return tFirstId( 0 );
    }

    // ----------------------------------------------------------------------------------

    moris_id
    Cut_Integration_Mesh::allocate_subphase_ids( moris::size_t aNumIdsToAllocate )
    {
        int tProcRank = moris::par_rank();
        int tProcSize = moris::par_size();

        // size_t is defined as uint here because of aNumRequested
        // Initialize gathered information outputs (information which will be scattered across processors)
        moris::Cell< moris_id > aGatheredInfo;
        moris::Cell< moris_id > tFirstId( 1 );
        moris::Cell< moris_id > tNumIdsRequested( 1 );

        // put current processors ID request size into the Cell that will be shared across procs
        tNumIdsRequested( 0 ) = (moris_id)aNumIdsToAllocate;

        // hand ID range size request to root processor
        moris::gather( tNumIdsRequested, aGatheredInfo );

        // initialize list holding the first ID in range for each processor
        moris::Cell< moris_id > tProcFirstID( tProcSize );

        // Subphase IDs up to the number of IP cells have already been used. Hence, the first free ID is:
        moris_index tFirstSubphaseId = mBackgroundMesh->get_max_entity_id( EntityRank::ELEMENT ) + 1;

        // Manage information on the root processor
        if ( tProcRank == 0 )
        {
            // Loop over entities print the number of entities requested by each processor
            for ( int iProc = 0; iProc < tProcSize; ++iProc )
            {
                // Give each processor their desired amount of IDs
                tProcFirstID( iProc ) = tFirstSubphaseId;

                // Increment the first available node ID
                tFirstSubphaseId = tFirstSubphaseId + aGatheredInfo( iProc );
            }
        }

        // on proc 0: split up the list of first IDs for every proc and send it to every other proc
        // on all procs: receive the assigned first SP ID as tFirstId
        moris::scatter( tProcFirstID, tFirstId );

        // return the first SP ID assigned
        return tFirstId( 0 );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_first_available_index( enum EntityRank aEntityRank ) const
    {
        MORIS_ERROR( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT, "Only can handle this question for nodes and elements" );

        if ( aEntityRank == EntityRank::NODE )
        {
            return mIntegrationVertices.size();
        }

        if ( aEntityRank == EntityRank::ELEMENT )
        {
            return mIntegrationCells.size();
        }

        return MORIS_INDEX_MAX;
    }

    // ----------------------------------------------------------------------------------

    uint
    Cut_Integration_Mesh::get_num_ig_cell_groups()
    {
        return mIntegrationCellGroups.size();
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< IG_Cell_Group >
    Cut_Integration_Mesh::get_ig_cell_group( moris_index aGroupIndex )
    {
        return mIntegrationCellGroups( aGroupIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const &
    Cut_Integration_Mesh::get_ig_cell_group_memberships( moris_index aIgCellIndex )
    {
        return mIntegrationCellToCellGroupIndex( aIgCellIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Cell*
    Cut_Integration_Mesh::get_ig_cell_group_parent_cell( moris_index aGroupIndex )
    {
        return mIntegrationCellGroupsParentCell( aGroupIndex );
    }

    // ----------------------------------------------------------------------------------

    enum CellTopology
    Cut_Integration_Mesh::get_child_element_topology()
    {
        return xtk::determine_cell_topology( this->get_spatial_dim(), mXTKModel->ig_element_order(), CellShape::SIMPLEX );
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::set_child_mesh_subphase(
            moris_index                 aCMIndex,
            moris::Cell< moris_index >& aSubphasesGroups )
    {
        moris::Cell< std::shared_ptr< IG_Cell_Group > > tIgCellSubphases( aSubphasesGroups.size() );

        for ( uint i = 0; i < aSubphasesGroups.size(); i++ )
        {
            tIgCellSubphases( i ) = mSubPhaseCellGroups( aSubphasesGroups( i ) );
        }

        mChildMeshes( aCMIndex )->set_subphase_groups( tIgCellSubphases );
    }

    // ----------------------------------------------------------------------------------

    uint
    Cut_Integration_Mesh::get_num_subphases()
    {
        return mSubPhaseCellGroups.size();
    }

    // ----------------------------------------------------------------------------------

    uint
    Cut_Integration_Mesh::get_num_subphase_groups( moris_index aMeshIndexInList )
    {
        return mBsplineMeshInfos( aMeshIndexInList )->get_num_SPGs();
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< std::shared_ptr< Child_Mesh_Experimental > >&
    Cut_Integration_Mesh::get_owned_child_meshes()
    {
        return mOwnedChildMeshes;
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index >&
    Cut_Integration_Mesh::get_owned_subphase_indices()
    {
        return mOwnedSubphaseGroupsInds;
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index >&
    Cut_Integration_Mesh::get_not_owned_subphase_indices()
    {
        return mNotOwnedSubphaseGroupsInds;
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Cell*
    Cut_Integration_Mesh::get_subphase_parent_cell( const moris_index aSubPhaseIndex )
    {
        return mSubPhaseParentCell( aSubPhaseIndex );
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< IG_Cell_Group >
    Cut_Integration_Mesh::get_subphase_ig_cells( moris_index aSubPhaseIndex )
    {
        return mSubPhaseCellGroups( aSubPhaseIndex );
    }

    // ----------------------------------------------------------------------------------

    const moris::Cell< moris_index >&
    Cut_Integration_Mesh::get_ig_cells_in_SPG( moris_index aMeshIndexInList, moris_index aSubphaseGroupIndex )
    {
        return mBsplineMeshInfos( aMeshIndexInList )->mSubphaseGroups( aSubphaseGroupIndex )->get_ig_cell_indices_in_group();
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_subphase_id( moris_index aSubphaseIndex )
    {
        return mSubPhaseIds( aSubphaseIndex );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_subphase_index( moris_id aSubphaseId )
    {
        auto tIter = mGlobalToLocalSubphaseMap.find( aSubphaseId );

        MORIS_ASSERT( tIter != mGlobalToLocalSubphaseMap.end(),
                "Subphase id not in map: %i. par_rank:%i",
                aSubphaseId,
                par_rank() );

        return tIter->second;
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_subphase_bulk_phase( moris_index aSubphaseIndex )
    {
        return mSubPhaseBulkPhase( aSubphaseIndex );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_subphase_group_id(
            moris_index aSubphaseGroupIndex,
            moris_index aBsplineMeshListIndex )
    {
        return mBsplineMeshInfos( aBsplineMeshListIndex )->mSubphaseGroups( aSubphaseGroupIndex )->get_id();
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_subphase_group_index(
            moris_index aSubphaseGroupId,
            moris_index aBsplineMeshListIndex )
    {
        // get the index for the SPG ID
        return mBsplineMeshInfos( aBsplineMeshListIndex )->get_index_for_spg_id( aSubphaseGroupId );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_subphase_group_bulk_phase(
            moris_index aSubphaseGroupIndex,
            moris_index aBsplineMeshListIndex )
    {
        // get representative SP from SPG
        moris_index tSubphaseIndex = mBsplineMeshInfos( aBsplineMeshListIndex )->mSubphaseGroups( aSubphaseGroupIndex )->get_SP_indices_in_group()( 0 );

        // get and return Bulk phase index of SPG through SP
        return mSubPhaseBulkPhase( tSubphaseIndex );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_ig_cell_subphase_index( moris_index aIgCellIndex )
    {
        return mIntegrationCellToSubphaseIndex( aIgCellIndex );
    }

    // ----------------------------------------------------------------------------------

    bool
    Cut_Integration_Mesh::parent_cell_has_children( moris_index aParentCellIndex )
    {
        return (bool)mParentCellHasChildren( aParentCellIndex );
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::finalize_cut_mesh_construction()
    {
        Tracer tTracer( "XTK", "Cut_Integration_Mesh", "finalize_cut_mesh_construction", mXTKModel->mVerboseLevel, 1 );

        this->deduce_ig_cell_group_ownership();

        this->assign_controlled_ig_cell_ids();
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const &
    Cut_Integration_Mesh::get_parent_cell_subphases( moris_index aParentCellIndex )
    {
        return mParentCellToSubphase( aParentCellIndex );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_vertex_parent_index( moris_index const & aVertexIndex )
    {
        return mIgVertexParentEntityIndex( aVertexIndex );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_vertex_parent_rank( moris_index const & aVertexIndex )
    {
        return mIgVertexParentEntityRank( aVertexIndex );
    }

    // ----------------------------------------------------------------------------------

    uint
    Cut_Integration_Mesh::get_num_entities(
            enum EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        switch ( aEntityRank )
        {
            case EntityRank::NODE:
            {
                return mIntegrationVertices.size();
                break;
            }
            case EntityRank::ELEMENT:
            {
                return mIntegrationCells.size();
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }
                return 0;
        }
    }

    // ----------------------------------------------------------------------------------

    uint
    Cut_Integration_Mesh::get_num_base_ip_cells() const
    {
        return mBackgroundMesh->get_num_elems();
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Facet_Based_Connectivity >
    Cut_Integration_Mesh::get_face_connectivity()
    {
        return mIgCellFaceConnectivity;
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::set_face_connectivity( std::shared_ptr< Facet_Based_Connectivity > aFaceConnectivity )
    {
        mIgCellFaceConnectivity = aFaceConnectivity;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Facet_Based_Ancestry >
    Cut_Integration_Mesh::get_face_ancestry()
    {
        return mIgCellFaceAncestry;
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::set_face_ancestry( std::shared_ptr< Facet_Based_Ancestry > aFaceAncestry )
    {
        mIgCellFaceAncestry = aFaceAncestry;
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::set_interface_facets( moris::Cell< moris_index >& aInterfaces )
    {
        mInterfaceFacets = aInterfaces;
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const &
    Cut_Integration_Mesh::get_interface_facets()
    {
        return mInterfaceFacets;
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::set_bulk_phase_to_bulk_phase_dbl_side_interface( moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > >& aBpToBpDblSideInterfaces )
    {
        mBpToBpDblSideInterfaces = aBpToBpDblSideInterfaces;
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > > const &
    Cut_Integration_Mesh::get_bulk_phase_to_bulk_phase_dbl_side_interface()
    {
        return mBpToBpDblSideInterfaces;
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::set_background_facet_to_child_facet_connectivity( moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & aBgToChildFacet )
    {
        mBGFacetToChildFacet = aBgToChildFacet;
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const &
    Cut_Integration_Mesh::get_background_facet_to_child_facet_connectivity()
    {
        return mBGFacetToChildFacet;
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::set_subphase_neighborhood( std::shared_ptr< Subphase_Neighborhood_Connectivity > aSubphaseNeighborhood )
    {
        mSubphaseNeighborhood = aSubphaseNeighborhood;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Subphase_Neighborhood_Connectivity >
    Cut_Integration_Mesh::get_subphase_neighborhood()
    {
        return mSubphaseNeighborhood;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Subphase_Neighborhood_Connectivity >
    Cut_Integration_Mesh::get_subphase_group_neighborhood( moris_index aMeshIndex )
    {
        return mSubphaseGroupNeighborhood( aMeshIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< Bspline_Mesh_Info* >&
    Cut_Integration_Mesh::get_bspline_mesh_info()
    {
        return mBsplineMeshInfos;
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::setup_glob_to_loc_subphase_map()
    {
        mGlobalToLocalSubphaseMap.clear();
        for ( uint i = 0; i < this->get_num_subphases(); i++ )
        {
            moris_id tSubphaseId = this->get_subphase_id( (moris_id)i );
            if ( tSubphaseId == MORIS_INDEX_MAX )
            {
                std::cout << "IPCell Owner = " << this->get_subphase_parent_cell( (moris_index)i )->get_owner() << " on " << par_rank() << std::endl;
            }
            MORIS_ASSERT( tSubphaseId != MORIS_INDEX_MAX, "Subphase id set to max" );
            MORIS_ASSERT( mGlobalToLocalSubphaseMap.find( tSubphaseId ) == mGlobalToLocalSubphaseMap.end(), "Subphase id already in map" );
            mGlobalToLocalSubphaseMap[ tSubphaseId ] = i;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::construct_spg_id_to_index_map( moris_index aBsplineMeshListIndex )
    {
        // quickly get access to the B-spline mesh info for the current B-spline mesh
        Bspline_Mesh_Info* tBsplineMeshInfo = mBsplineMeshInfos( aBsplineMeshListIndex );

        // get access to the map stored in the B-spline mesh info
        std::unordered_map< moris_id, moris_index >& tSpgIdToIndexMap = tBsplineMeshInfo->mSpgIdToIndexMap;

        // clear any information that may already be in it
        tSpgIdToIndexMap.clear();

        // get the number of SPGs on the current B-spline mesh and processor
        uint tNumSpgs = tBsplineMeshInfo->get_num_SPGs();

        // loop over all SPGs and put their IDs in the map
        for ( uint iSPG = 0; iSPG < tNumSpgs; iSPG++ )
        {
            // get the ID of the currently treated SPG
            moris_id tSpgId = tBsplineMeshInfo->get_id_for_spg_index( iSPG );

            // check validity of SPG ID
            MORIS_ERROR( tSpgId != MORIS_ID_MAX,
                    "Cut_Integration_Mesh::construct_spg_id_to_index_map() - "
                    "Subphase Group ID not set. Should be set when this function is called." );

            MORIS_ASSERT( tSpgIdToIndexMap.find( tSpgId ) == tSpgIdToIndexMap.end(),
                    "Cut_Integration_Mesh::construct_spg_id_to_index_map() - "
                    "Subphase id already in map" );

            // fill map
            tSpgIdToIndexMap[ tSpgId ] = iSPG;
        }
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_cell_bulk_phase( moris_index aCellIndex )
    {
        return mIntegrationCellBulkPhase( aCellIndex );
    }

    // ----------------------------------------------------------------------------------

    Cell< moris_index >
    Cut_Integration_Mesh::register_side_set_names( moris::Cell< std::string > const & aSideSetNames )
    {
        uint tNumSetsToRegister = aSideSetNames.size();

        // block set ords
        Cell< moris_index > tSideSetOrds( tNumSetsToRegister );

        // iterate and add sets
        for ( uint i = 0; i < tNumSetsToRegister; i++ )
        {
            tSideSetOrds( i ) = mSideSetLabels.size();

            mSideSetLabels.push_back( aSideSetNames( i ) );
            mSideSetCellSides.push_back( nullptr );

            MORIS_ASSERT( mSideSideSetLabelToOrd.find( aSideSetNames( i ) ) == mSideSideSetLabelToOrd.end(),
                    "Duplicate side set in mesh" );

            mSideSideSetLabelToOrd[ aSideSetNames( i ) ] = tSideSetOrds( i );
        }
        return tSideSetOrds;
    }

    // ----------------------------------------------------------------------------------

    Cell< moris_index >
    Cut_Integration_Mesh::register_block_set_names(
            moris::Cell< std::string > const & aBlockSetNames,
            enum CellTopology                  aCellTopo )
    {
        uint tNumSetsToRegister = aBlockSetNames.size();

        // block set ords
        Cell< moris_index > tBlockOrds( tNumSetsToRegister );

        // iterate and add sets
        for ( uint i = 0; i < tNumSetsToRegister; i++ )
        {
            tBlockOrds( i ) = mBlockSetNames.size();

            mBlockSetNames.push_back( aBlockSetNames( i ) );
            mBlockSetCellGroup.push_back( nullptr );
            mBlockCellTopo.push_back( aCellTopo );

            MORIS_ASSERT( mBlockSetLabelToOrd.find( aBlockSetNames( i ) ) == mBlockSetLabelToOrd.end(),
                    "Duplicate block set in mesh" );

            mBlockSetLabelToOrd[ aBlockSetNames( i ) ] = tBlockOrds( i );
        }

        return tBlockOrds;
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::write_mesh(
            std::string aOutputPath,
            std::string aOutputFile )
    {
        Tracer tTracer( "XTK", "Cut Integration Mesh", "Write mesh", mXTKModel->mVerboseLevel, 0 );
        // get path to output XTK files to
        std::string tOutputPath = aOutputPath;
        std::string tOutputFile = aOutputFile;
        std::string tOutputBase = tOutputFile.substr( 0, tOutputFile.find( "." ) );
        std::string tOutputExt  = tOutputFile.substr( tOutputFile.find( "." ), tOutputFile.length() );

        MORIS_ASSERT( tOutputExt == ".exo" || tOutputExt == ".e", "Invalid file extension, needs to be .exo or .e" );

        // Write mesh
        moris::mtk::Writer_Exodus writer( this );

        writer.write_mesh(
                "", tOutputPath + tOutputFile, "", tOutputPath + "xtk_temp2.exo" );

        // Write the fields
        writer.set_time( 0.0 );
        writer.close_file();
    }

    // ----------------------------------------------------------------------------------

    Cell_Connectivity
    Cut_Integration_Mesh::get_background_cell_connectivity( moris_index aBGCellId ) const
    {
        if ( this->get_spatial_dim() == 3 )
        {
            return Cell_Connectivity(
                    mBackgroundMesh->get_entity_connected_to_entity_loc_inds( aBGCellId, EntityRank::ELEMENT, EntityRank::NODE ),
                    mBackgroundMesh->get_entity_connected_to_entity_loc_inds( aBGCellId, EntityRank::ELEMENT, EntityRank::EDGE ),
                    mBackgroundMesh->get_entity_connected_to_entity_loc_inds( aBGCellId, EntityRank::ELEMENT, EntityRank::FACE ) );
        }
        else if ( this->get_spatial_dim() == 2 )
        {
            return Cell_Connectivity(
                    mBackgroundMesh->get_entity_connected_to_entity_loc_inds( aBGCellId, EntityRank::ELEMENT, EntityRank::NODE ),
                    mBackgroundMesh->get_entity_connected_to_entity_loc_inds( aBGCellId, EntityRank::ELEMENT, EntityRank::EDGE ),
                    Matrix< IndexMat >( 0, 0 ) );
        }
        else
        {
            return Cell_Connectivity();
        }
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Cut_Integration_Mesh::get_coarsest_bspline_mesh_index_on_base_ip_cell( moris_index aBaseIpCellIndex ) const
    {
        // check that the input makes sense
        MORIS_ASSERT( (uint)aBaseIpCellIndex < mCoarsestBsplineMesh.size(),
                "Cut_Integration_Mesh::get_coarsest_bspline_mesh_index_on_base_ip_cell() - "
                "Map has not been constructed yet, or base IP cell index out of bounds." );

        // return the B-spline mesh index
        return mCoarsestBsplineMesh( aBaseIpCellIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const &
    Cut_Integration_Mesh::get_union_MSD_indices_for_base_IP_cell( const moris_index aBaseIpCellIndex ) const
    {
        MORIS_ASSERT( mUnionVoidMsdIndices.size() > 0,
                "Cut_Integration_Mesh::get_union_MSD_indices_for_base_IP_cell() - information not constructed yet" );
        MORIS_ASSERT( mUnionVoidMsdIndices.size() > (uint)aBaseIpCellIndex,
                "Cut_Integration_Mesh::get_union_MSD_indices_for_base_IP_cell() - IP cell index out of bounds" );
        return mUnionVoidMsdIndices( aBaseIpCellIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const &
    Cut_Integration_Mesh::get_bulk_phases_for_union_MSD_indices_for_base_IP_cell( const moris_index aBaseIpCellIndex ) const
    {
        MORIS_ASSERT( mUnionVoidMsdIndexBulkPhases.size() > 0,
                "Cut_Integration_Mesh::get_bulk_phases_for_union_MSD_indices_for_base_IP_cell() - information not constructed yet" );
        MORIS_ASSERT( mUnionVoidMsdIndexBulkPhases.size() > (uint)aBaseIpCellIndex,
                "Cut_Integration_Mesh::get_bulk_phases_for_union_MSD_indices_for_base_IP_cell() - IP cell index out of bounds" );
        return mUnionVoidMsdIndexBulkPhases( aBaseIpCellIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const &
    Cut_Integration_Mesh::get_material_SPG_indices_for_base_IP_cell(
            const moris_index aBsplineMeshListIndex,
            const moris_index aBaseIpCellIndex ) const
    {
        MORIS_ASSERT( mBsplineMeshInfos( aBsplineMeshListIndex )->mExtractionCellMaterialSpgs.size() > 0,
                "Cut_Integration_Mesh::get_material_SPG_indices_for_base_IP_cell() - information not constructed yet" );
        return mBsplineMeshInfos( aBsplineMeshListIndex )->mExtractionCellMaterialSpgs( aBaseIpCellIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const &
    Cut_Integration_Mesh::get_material_MSD_indices_for_base_IP_cell(
            const moris_index aBsplineMeshListIndex,
            const moris_index aBaseIpCellIndex ) const
    {
        MORIS_ASSERT( mBsplineMeshInfos( aBsplineMeshListIndex )->mExtractionCellMaterialMsdIndices.size() > 0,
                "Cut_Integration_Mesh::get_material_MSD_indices_for_base_IP_cell() - information not constructed yet" );
        return mBsplineMeshInfos( aBsplineMeshListIndex )->mExtractionCellMaterialMsdIndices( aBaseIpCellIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const &
    Cut_Integration_Mesh::get_void_SPG_indices_for_base_IP_cell(
            const moris_index aBsplineMeshListIndex,
            const moris_index aBaseIpCellIndex ) const
    {
        MORIS_ASSERT( mBsplineMeshInfos( aBsplineMeshListIndex )->mExtractionCellVoidSpgs.size() > 0,
                "Cut_Integration_Mesh::get_void_SPG_indices_for_base_IP_cell() - information not constructed yet" );
        return mBsplineMeshInfos( aBsplineMeshListIndex )->mExtractionCellVoidSpgs( aBaseIpCellIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const &
    Cut_Integration_Mesh::get_void_MSD_indices_for_base_IP_cell(
            const moris_index aBsplineMeshListIndex,
            const moris_index aBaseIpCellIndex ) const
    {
        MORIS_ASSERT( mBsplineMeshInfos( aBsplineMeshListIndex )->mExtractionCellVoidMsdIndices.size() > 0,
                "Cut_Integration_Mesh::get_void_MSD_indices_for_base_IP_cell() - information not constructed yet" );
        return mBsplineMeshInfos( aBsplineMeshListIndex )->mExtractionCellVoidMsdIndices( aBaseIpCellIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const &
    Cut_Integration_Mesh::get_free_void_MSD_indices_for_base_IP_cell(
            const moris_index aBsplineMeshListIndex,
            const moris_index aBaseIpCellIndex ) const
    {
        MORIS_ASSERT( mBsplineMeshInfos( aBsplineMeshListIndex )->mExtractionCellFreeVoidMsdIndices.size() > 0,
                "Cut_Integration_Mesh::get_free_void_MSD_indices_for_base_IP_cell() - information not constructed yet" );
        return mBsplineMeshInfos( aBsplineMeshListIndex )->mExtractionCellFreeVoidMsdIndices( aBaseIpCellIndex );
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::print_cells(
            bool        aOmitIndex,
            std::string aFile )
    {

        std::ostringstream tStringStream;
        // // Number of Controlled Cells:
        if ( aFile.empty() == true )
        {
            tStringStream << "Num IG Cells: " << this->mIntegrationCells.size() << std::endl;
            tStringStream << "Num Controlled IG Cells: " << this->mControlledIgCells.size() << std::endl;
        }
        // max num verts to cells
        uint tMaxVertsToCell = 0;
        for ( uint i = 0; i < this->get_num_entities( EntityRank::ELEMENT, 0 ); i++ )
        {
            mtk::Cell& tCell = this->get_mtk_cell( (moris_index)i );
            if ( tCell.get_number_of_vertices() > tMaxVertsToCell )
            {
                tMaxVertsToCell = tCell.get_number_of_vertices();
            }
        }

        tStringStream << "Cell_Id,";
        if ( !aOmitIndex ) { tStringStream << "Cell_Ind,"; }
        tStringStream << "Owner,";
        tStringStream << "PRank,";
        tStringStream << "Phase,";
        tStringStream << "Measure,";
        for ( uint iVH = 0; iVH < tMaxVertsToCell; iVH++ )
        {
            tStringStream << "Vert_" + std::to_string( iVH );

            if ( iVH != tMaxVertsToCell - 1 )
            {
                tStringStream << ",";
            }
        }
        tStringStream << "\n";

        for ( uint i = 0; i < this->get_num_entities( EntityRank::ELEMENT, 0 ); i++ )
        {
            mtk::Cell&                         tCell     = this->get_mtk_cell( (moris_index)i );
            moris::Cell< moris::mtk::Vertex* > tVertices = tCell.get_vertex_pointers();

            tStringStream << std::to_string( tCell.get_id() ) + ",";
            if ( !aOmitIndex ) { tStringStream << std::to_string( tCell.get_index() ) + ","; }
            tStringStream << std::to_string( tCell.get_owner() ) << ",";
            tStringStream << std::to_string( par_rank() ) << ",";
            tStringStream << std::to_string( this->get_cell_bulk_phase( i ) ) + ",";
            tStringStream << std::scientific << tCell.compute_cell_measure() << ",";

            for ( uint j = 0; j < tMaxVertsToCell; j++ )
            {
                if ( j < tVertices.size() )
                {
                    tStringStream << std::to_string( tVertices( j )->get_id() );
                }
                else
                {
                    tStringStream << std::to_string( MORIS_INDEX_MAX );
                }

                if ( j != tMaxVertsToCell - 1 )
                {
                    tStringStream << ",";
                }
            }
            tStringStream << "\n";
        }

        if ( aFile.empty() == false )
        {
            std::ofstream tOutputFile( aFile );
            tOutputFile << tStringStream.str() << std::endl;
            tOutputFile.close();
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::print_vertices(
            bool        aOmitIndex,
            std::string aFile )
    {
        std::ostringstream tStringStream;
        tStringStream.clear();
        tStringStream.str( "" );

        tStringStream << "Vert_Id,";
        if ( !aOmitIndex ) { tStringStream << "Vert Ind,"; }
        tStringStream << "Owner,";
        tStringStream << "Prank,";

        for ( uint iVH = 0; iVH < this->get_spatial_dim(); iVH++ )
        {
            tStringStream << "Coords_" + std::to_string( iVH );

            if ( iVH != this->get_spatial_dim() - 1 )
            {
                tStringStream << ",";
            }
        }

        tStringStream << std::endl;

        for ( uint i = 0; i < this->get_num_entities( EntityRank::NODE, 0 ); i++ )
        {
            mtk::Vertex& tVertex = this->get_mtk_vertex( (moris_index)i );
            tStringStream.precision( 16 );

            tStringStream << tVertex.get_id() << ",";
            if ( !aOmitIndex ) { tStringStream << tVertex.get_index() << ","; }
            tStringStream << tVertex.get_owner() << ",";
            tStringStream << par_rank() << ",";

            moris::Matrix< moris::DDRMat > tCoords = tVertex.get_coords();

            for ( uint iSp = 0; iSp < this->get_spatial_dim(); iSp++ )
            {
                tStringStream << std::scientific << tCoords( iSp );

                if ( iSp != this->get_spatial_dim() - 1 )
                {
                    tStringStream << ",";
                }
            }

            //
            tStringStream << std::endl;
        }
        if ( aFile.empty() == false )
        {
            std::ofstream tOutputFile( aFile );
            tOutputFile << tStringStream.str() << std::endl;
            tOutputFile.close();
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::print_groupings( std::string aFile )
    {
        std::ostringstream tStringStream;

        // headers
        tStringStream << "IP_Cell_Id,";
        tStringStream << "PRank,";
        tStringStream << "Group_Index,";

        // global max size of
        moris_index tLocalMaxIGCellGroupSize = 0;
        for ( uint iGroups = 0; iGroups < mIntegrationCellGroups.size(); iGroups++ )
        {
            if ( (moris_index)mIntegrationCellGroups( iGroups )->mIgCellGroup.size() > tLocalMaxIGCellGroupSize )
            {
                tLocalMaxIGCellGroupSize = mIntegrationCellGroups( iGroups )->mIgCellGroup.size();
            }
        }

        moris_index tGlbMaxIgCellGroupSize = moris::max_all( tLocalMaxIGCellGroupSize );

        for ( moris_index iCH = 0; iCH < tGlbMaxIgCellGroupSize; iCH++ )
        {
            tStringStream << "IG_Cell_ID_" + std::to_string( iCH );

            if ( iCH != tGlbMaxIgCellGroupSize - 1 )
            {
                tStringStream << ",";
            }
        }

        tStringStream << "\n";

        // iterate through the groups
        for ( uint iGroup = 0; iGroup < mIntegrationCellGroups.size(); iGroup++ )
        {
            tStringStream << mIntegrationCellGroupsParentCell( iGroup )->get_id() << ",";
            tStringStream << moris::par_rank() << ",";
            tStringStream << iGroup << ",";

            for ( moris_index iGC = 0; iGC < (moris_index)mIntegrationCellGroups( iGroup )->mIgCellGroup.size(); iGC++ )
            {
                tStringStream << mIntegrationCellGroups( iGroup )->mIgCellGroup( iGC )->get_id();

                if ( iGC != tGlbMaxIgCellGroupSize - 1 )
                {
                    tStringStream << ",";
                }
            }

            tStringStream << "\n";
        }

        if ( aFile.empty() == false )
        {
            std::ofstream tOutputFile( aFile );
            tOutputFile << tStringStream.str() << std::endl;
            tOutputFile.close();
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::setup_comm_map()
    {
        if ( par_size() > 1 )
        {
            // Build list of processors with whom current processor shares nodes
            std::unordered_map< moris_id, moris_id > tProcList;
            Cell< moris_index >                      tCellOfProcs;

            // Loop over all nodes in background mesh and get node's owner
            for ( uint i = 0; i < mBackgroundMesh->get_num_entities( EntityRank::NODE ); i++ )
            {
                moris_index tOwner = mBackgroundMesh->get_entity_owner( (moris_index)i, EntityRank::NODE );

                if ( tProcList.find( tOwner ) == tProcList.end() && tOwner != par_rank() )
                {
                    tCellOfProcs.push_back( tOwner );
                    tProcList[ tOwner ] = 1;
                }
            }

            // Initialize communication map
            mCommunicationMap.resize( 1, tCellOfProcs.size() );

            for ( uint i = 0; i < tCellOfProcs.size(); i++ )
            {
                mCommunicationMap( i ) = tCellOfProcs( i );
            }

            // Send communication maps to root processor (here 0)
            Cell< Matrix< IndexMat > > tGatheredMats;
            moris_index                tTag = 10009;
            if ( tCellOfProcs.size() == 0 )
            {
                Matrix< IndexMat > tDummy( 1, 1, MORIS_INDEX_MAX );
                all_gather_vector( tDummy, tGatheredMats, tTag, 0, 0 );
            }
            else
            {
                all_gather_vector( mCommunicationMap, tGatheredMats, tTag, 0, 0 );
            }

            // Initialize processor-to-processor communication table
            Cell< Matrix< IndexMat > > tReturnMats( par_size() );

            // Root processor (here 0) builds processor-to-processor communication table
            if ( par_rank() == 0 )
            {
                Cell< Cell< uint > > tProcToProc( par_size() );

                for ( uint i = 0; i < tGatheredMats.size(); i++ )
                {
                    for ( uint j = 0; j < tGatheredMats( i ).numel(); j++ )
                    {
                        if ( tGatheredMats( i )( j ) != MORIS_INDEX_MAX )
                        {
                            tProcToProc( tGatheredMats( i )( j ) ).push_back( (moris_index)i );
                        }
                    }
                }

                // convert to a matrix
                for ( uint i = 0; i < (uint)par_size(); i++ )
                {
                    tReturnMats( i ).resize( 1, tProcToProc( i ).size() );

                    for ( uint j = 0; j < tProcToProc( i ).size(); j++ )
                    {
                        tReturnMats( i )( j ) = tProcToProc( i )( j );
                    }
                    if ( tProcToProc( i ).size() == 0 )
                    {
                        tReturnMats( i ) = Matrix< IndexMat >( 1, 1, MORIS_INDEX_MAX );
                    }
                }

                // send processor-to-processor communication table back individual processors
                for ( uint i = 0; i < (uint)par_size(); i++ )
                {
                    nonblocking_send(
                            tReturnMats( i ),
                            tReturnMats( i ).n_rows(),
                            tReturnMats( i ).n_cols(),
                            i,
                            tTag );
                }
            }

            // receive processor-to-processor communication tables
            barrier();
            Matrix< IndexMat > tTempCommMap( 1, 1, 0 );
            receive( tTempCommMap, 1, 0, tTag );

            // add new processors to existing communication table
            for ( uint i = 0; i < tTempCommMap.numel(); i++ )
            {
                // Skip processors that do not share nodes
                if ( tTempCommMap( i ) != MORIS_INDEX_MAX )
                {
                    // Check if processor is already in communication table
                    if ( tProcList.find( tTempCommMap( i ) ) == tProcList.end()
                            && tTempCommMap( i ) != par_rank() )
                    {
                        moris_index tIndex = mCommunicationMap.numel();

                        mCommunicationMap.resize( 1, mCommunicationMap.numel() + 1 );

                        tProcList[ tTempCommMap( i ) ] = 1;

                        mCommunicationMap( tIndex ) = tTempCommMap( i );
                    }
                }
            }
            barrier();
        }

        // mark the index map as needing to be re-constructed
        mCommMapHasBeenConstructed = false;
    }

    // ----------------------------------------------------------------------------------
    void
    Cut_Integration_Mesh::deduce_ig_cell_group_ownership()
    {
        Tracer tTracer( "XTK", "Cut_Integration_Mesh", "deduce_ig_cell_group_ownership", mXTKModel->mVerboseLevel, 1 );

        mOwnedIntegrationCellGroupsInds.reserve( mIntegrationCellGroupsParentCell.size() );
        mNotOwnedIntegrationCellGroups.reserve( mIntegrationCellGroupsParentCell.size() );

        mOwnedChildMeshes.reserve( mChildMeshes.size() );
        mNotOwnedChildMeshes.reserve( mChildMeshes.size() );

        moris_index tParRank = moris::par_rank();

        // iterate through groups
        for ( uint iParentCell = 0; iParentCell < mIntegrationCellGroupsParentCell.size(); iParentCell++ )
        {
            if ( mIntegrationCellGroupsParentCell( iParentCell )->get_owner() == tParRank )
            {
                mOwnedIntegrationCellGroupsInds.push_back( (moris_index)iParentCell );
                mOwnedChildMeshes.push_back( mChildMeshes( iParentCell ) );
            }
            else
            {
                mNotOwnedIntegrationCellGroups.push_back( (moris_index)iParentCell );
                mNotOwnedChildMeshes.push_back( mChildMeshes( iParentCell ) );
            }
        }

        mOwnedIntegrationCellGroupsInds.shrink_to_fit();
        mNotOwnedIntegrationCellGroups.shrink_to_fit();
        mOwnedChildMeshes.shrink_to_fit();
        mNotOwnedChildMeshes.shrink_to_fit();
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::assign_controlled_ig_cell_ids()
    {
        // log this function when verbose output is requested
        Tracer tTracer( "XTK", "Cut Integration Mesh", "assign IG cell IDs", mXTKModel->mVerboseLevel, 1 );

        // get the communication table
        Matrix< IdMat > tCommTable     = this->get_communication_table();

        /* ---------------------------------------------------------------------------------------- */
        /* Step 1: Let each proc decide how many entity IDs it needs & communicate ID ranges */

        // get the number of owned IG cells
        uint tNumControlledCellsInCutMesh = mControlledIgCells.size();

        // reserve IDs for his proc
        moris_id tMyFirstId = (moris_id)get_processor_offset( tNumControlledCellsInCutMesh ) + mGlobalMaxCellId;

        // update the maximum ID allocated
        moris_id tMyMaxID = tMyFirstId + (moris_id)tNumControlledCellsInCutMesh;
        mGlobalMaxCellId  = moris::max_all( tMyMaxID );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 2: Assign IDs to owned entities */

        // assign IDs to owned IG cells
        this->assign_IDs_to_owned_IG_cells( tMyFirstId );

        /* ---------------------------------------------------------------------------------------- */
        /* The following steps are only necessary if code runs in parallel */

        if ( par_size() == 1 )    // serial
        {
            // check that all IG cells are owned in serial
            MORIS_ASSERT( mNotOwnedIntegrationCellGroups.size() == 0,
                    "Cut_Integration_Mesh::assign_controlled_ig_cell_ids() - "
                    "Code running in serial, but not all IG cell groups are owned by proc 0." );
        }
        else    // parallel
        {

            /* ---------------------------------------------------------------------------------------- */
            /* Step 3: Prepare requests for non-owned entities */

            // initialize lists of information that identifies IG cells (on other procs)
            Cell< Cell< moris_index > > tNotOwnedIgCellGroups;      // IG cell group index (local to current proc, just used for construction of arrays)
            Cell< Matrix< IdMat > >     tParentCellIds;             // IDs of the IG cells' parent cells
            Cell< Matrix< IndexMat > >  tNumIgCellsInParentCell;    // Number of IG cells in parent cell

            // fill the identifying information
            this->prepare_requests_for_not_owned_IG_cell_IDs( tNotOwnedIgCellGroups, tParentCellIds, tNumIgCellsInParentCell );

            /* ---------------------------------------------------------------------------------------- */
            /* Step 4: Send and Receive requests about non-owned entities to and from other procs */

            // initialize arrays for receiving
            Cell< Matrix< IdMat > >    tReceivedParentCellIds;
            Cell< Matrix< IndexMat > > tReceivedNumIgCellsInParentCell;

            // communicate information
            moris::communicate_mats( tCommTable, tParentCellIds, tReceivedParentCellIds );
            moris::communicate_mats( tCommTable, tNumIgCellsInParentCell, tReceivedNumIgCellsInParentCell );

            // clear memory not needed anymore
            tParentCellIds.clear();
            tNumIgCellsInParentCell.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 5: Find answers to the requests */

            // initialize lists of ID answers to other procs
            Cell< Matrix< IdMat > > tFirstIgCellIdsInCellGroups;

            // find the IG cell IDs requested by the other procs
            this->prepare_answers_for_owned_IG_cell_IDs( 
                    tFirstIgCellIdsInCellGroups,
                    tReceivedParentCellIds,
                    tReceivedNumIgCellsInParentCell );

            // clear memory from requests (the answers to which have been found)
            tReceivedParentCellIds.clear();
            tReceivedNumIgCellsInParentCell.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 6: Send and receive answers to and from other procs */

            // initialize arrays for receiving
            Cell< Matrix< IdMat > > tReceivedFirstIgCellIdsInCellGroups;

            // communicate answers
            moris::communicate_mats( tCommTable, tFirstIgCellIdsInCellGroups, tReceivedFirstIgCellIdsInCellGroups );

            // clear unused memory
            tFirstIgCellIdsInCellGroups.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 7: Use answers to assign IDs to non-owned entities */

            this->handle_requested_IG_cell_ID_answers( tNotOwnedIgCellGroups, tReceivedFirstIgCellIdsInCellGroups );

        }    // end if: parallel

    }    // end function: Cut_Integration_Mesh::assign_controlled_ig_cell_ids()

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::assign_IDs_to_owned_IG_cells( moris_id aFirstFreeId )
    {
        // make sure this function is not called before the owned IG cells have been determined
        MORIS_ASSERT( mOwnedIntegrationCellGroupsInds.size() > 0 ,
                "Cut_Integration_Mesh::assign_IDs_to_owned_IG_cells() - "
                "No IG cells are owned. This is likely because the owned IG cells have not been determined yet. " );

        // set child elements ids in the children meshes which the current proc owns and does not share
        for ( moris::size_t iCellGroup = 0; iCellGroup < mOwnedIntegrationCellGroupsInds.size(); iCellGroup++ )
        {
            // get the pointer to the current cell-group
            std::shared_ptr< IG_Cell_Group > tCellGroup = this->get_ig_cell_group( mOwnedIntegrationCellGroupsInds( iCellGroup ) );

            // iterate through the child cell elements in the group
            for ( uint iCellInGroup = 0; iCellInGroup < tCellGroup->mIgCellGroup.size(); iCellInGroup++ )
            {
                // get access to the current IG cell
                mtk::Cell const * tIgCell      = tCellGroup->mIgCellGroup( iCellInGroup );
                moris_index       tIgCellIndex = tIgCell->get_index();

                // get the current IG cell's index
                moris_index tIgCellControlledIndex = this->get_integration_cell_controlled_index( tIgCellIndex );

                // get the IG cell
                mControlledIgCells( tIgCellControlledIndex )->set_id( aFirstFreeId );

                // check that the ID
                MORIS_ASSERT(
                        aFirstFreeId == tIgCell->get_id(),
                        "Cut_Integration_Mesh::assign_controlled_ig_cell_ids() - "
                        "ID reported by IG cell different from ID just assigned to corresponding controlled IG cell." );
                MORIS_ASSERT(
                        mIntegrationCellIdToIndexMap.find( tIgCell->get_id() ) == mIntegrationCellIdToIndexMap.end(),
                        "Cut_Integration_Mesh::assign_controlled_ig_cell_ids() - "
                        "IG cell's ID already in the map, i.e. it has already been assigned before." );

                // populate ID to index map for IG cells
                mIntegrationCellIdToIndexMap[ aFirstFreeId ] = tIgCellIndex;

                // increment ID counter
                aFirstFreeId++;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::prepare_requests_for_not_owned_IG_cell_IDs( 
            Cell< Cell< moris_index > >& aNotOwnedIgCellGroups,  
            Cell< Matrix< IdMat > >&     aParentCellIds,
            Cell< Matrix< IndexMat > >&  aNumIgCellsInParentCell )
    {
        // get the communication table and map
        Matrix< IdMat > tCommTable     = this->get_communication_table();
        uint            tCommTableSize = tCommTable.numel();
        std::map< moris_id, moris_index > tProcIdToCommTableIndex = this->get_communication_map();

        // initialize lists of identifying information
        aNotOwnedIgCellGroups.resize( tCommTableSize );
        aParentCellIds.resize( tCommTableSize );
        aNumIgCellsInParentCell.resize( tCommTableSize );

        // get the number of IG cell group for which IG cells need to be communicated
        uint tNumNotOwnedIgCellGroup = mNotOwnedIntegrationCellGroups.size();

        // make sure this function is run in parallel only
        MORIS_ASSERT( tNumNotOwnedIgCellGroup > 0, 
                "Cut_Integration_Mesh::prepare_requests_for_not_owned_IG_cell_IDs() - "
                "No not owned IG cells assigned. Either the owned and not owned IG cells have not been established yet; "
                "or this function is called in serial (it should only be called in parallel)." );

        // sort non-owned parent cells (corresponding to IG cell groups) into lists associated with each of the processors communicated with
        for ( uint iNonOwnedCM = 0; iNonOwnedCM < tNumNotOwnedIgCellGroup; iNonOwnedCM++ )
        {
            // get the index of the current non-owned IG cell group
            moris_index tNonOwnedCellGroupIndex = mNotOwnedIntegrationCellGroups( iNonOwnedCM );

            // find the position of the current proc in the communication table
            moris_index tOwnerProc = mIntegrationCellGroupsParentCell( tNonOwnedCellGroupIndex )->get_owner();
            auto        tIter        = tProcIdToCommTableIndex.find( tOwnerProc );
            MORIS_ASSERT(
                    tIter != tProcIdToCommTableIndex.end(),
                    "Cut_Integration_Mesh::assign_controlled_ig_cell_ids() - "
                    "IG cell group owner (Proc #%i) not found in communication table of current proc #%i which is: %s",
                    tOwnerProc,
                    par_rank(),
                    ios::stringify_log( tCommTable ).c_str() );
            moris_index tProcDataIndex = tIter->second;

            // add the IG cell group to the list of child-meshes needing to be exchanged with this processor
            aNotOwnedIgCellGroups( tProcDataIndex ).push_back( tNonOwnedCellGroupIndex );
        }

        // populate the identifying information for each non-owned IG cell group (for each processor to communicate with)
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // number of IG cell group shared with this processor
            uint tNumIgCellGroupsOnProc = aNotOwnedIgCellGroups( iProc ).size();

            // resize matrix to accommodate all information
            aParentCellIds( iProc ).resize( tNumIgCellGroupsOnProc, 1 );
            aNumIgCellsInParentCell( iProc ).resize( tNumIgCellGroupsOnProc, 1 );

            // populate the identifying information for each non-owned IG cell group (for each IG cell group on each proc)
            for ( uint iIgCellGroup = 0; iIgCellGroup < tNumIgCellGroupsOnProc; iIgCellGroup++ )
            {
                // get the index of the IG cell group treated
                moris_index tIGCellGroupIndex = aNotOwnedIgCellGroups( iProc )( iIgCellGroup );

                // find and store the information for communication
                aParentCellIds( iProc )( iIgCellGroup )          = mIntegrationCellGroupsParentCell( tIGCellGroupIndex )->get_id();
                aNumIgCellsInParentCell( iProc )( iIgCellGroup ) = mIntegrationCellGroups( tIGCellGroupIndex )->mIgCellGroup.size();
            }
        }

        // size out unused memory
        aNotOwnedIgCellGroups.shrink_to_fit();
        aParentCellIds.shrink_to_fit();
        aNumIgCellsInParentCell.shrink_to_fit();
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::prepare_answers_for_owned_IG_cell_IDs( 
            Cell< Matrix< IdMat > >&     aFirstIgCellIdsInCellGroups,
            Cell< Matrix< IdMat > > const&     aReceivedParentCellIds,
            Cell< Matrix< IndexMat > > const&  aReceivedNumIgCellsInParentCell )
    {
        // initialize array of ID answers to other procs with correct size
        aFirstIgCellIdsInCellGroups.resize( aReceivedParentCellIds.size() );

        // answer requests from each proc
        for ( uint iProcInCommTable = 0; iProcInCommTable < aReceivedParentCellIds.size(); iProcInCommTable++ )
        {
            // get the number of CMs for which IG cell IDs need to be communicated for the current proc
            uint tNumIgCellGroupsCommunicatedWithProc = aReceivedParentCellIds( iProcInCommTable ).numel();

            // resize answer arrays
            aFirstIgCellIdsInCellGroups( iProcInCommTable ).set_size( tNumIgCellGroupsCommunicatedWithProc, 1 );

            // go through and answer for all CMs requested by the current proc
            for ( uint iIgCellGroup = 0; iIgCellGroup < tNumIgCellGroupsCommunicatedWithProc; iIgCellGroup++ )
            {
                // get the parent Cell ID for the current request
                moris_id tParentId = aReceivedParentCellIds( iProcInCommTable )( iIgCellGroup );

                // get this parent Cells index wrt to the executing proc
                moris_index tParentCellIndex = mBackgroundMesh->get_loc_entity_ind_from_entity_glb_id( tParentId, EntityRank::ELEMENT );

                // get the index of the attached IG cell group/Child mesh
                moris_index tIgCellGroupIndex = mParentCellCellGroupIndex( tParentCellIndex );

                // check the request
                MORIS_ASSERT( 
                        tIgCellGroupIndex != MORIS_INDEX_MAX,
                        "Cut_Integration_Mesh::prepare_answers_for_owned_IG_cell_IDs() - "
                        "Request is made for child element IDs on a parent cell not intersected" );
                MORIS_ASSERT( 
                        par_rank() == mIntegrationCellGroupsParentCell( tIgCellGroupIndex )->get_owner(),
                        "Cut_Integration_Mesh::prepare_answers_for_owned_IG_cell_IDs() - "
                        "Current proc does not own this entity that had info requested." );
                MORIS_ASSERT( 
                        mIntegrationCellGroups( tIgCellGroupIndex )->mIgCellGroup.size() == (uint)aReceivedNumIgCellsInParentCell( iProcInCommTable )( iIgCellGroup ),
                        "Cut_Integration_Mesh::prepare_answers_for_owned_IG_cell_IDs() - "
                        "Proc #%i: %i IG cells are in integration cell group %i, but the parent cell is marked to have %i IG cells (for communication).",
                        par_rank(),
                        mIntegrationCellGroups( tIgCellGroupIndex )->mIgCellGroup.size(),
                        tIgCellGroupIndex,
                        aReceivedNumIgCellsInParentCell( iProcInCommTable )( iIgCellGroup ) );

                // answer the request
                if ( mIntegrationCellGroups( tIgCellGroupIndex )->mIgCellGroup.size() > 0 )
                {
                    // get the first index in the child group
                    aFirstIgCellIdsInCellGroups( iProcInCommTable )( iIgCellGroup ) = mIntegrationCellGroups( tIgCellGroupIndex )->mIgCellGroup( 0 )->get_id();
                }
                else
                {
                    // return a default for non-decomposed background cells (as there are no 'controlled' IG cells that need an ID)
                    aFirstIgCellIdsInCellGroups( iProcInCommTable )( iIgCellGroup ) = MORIS_ID_MAX;
                }
            } // end for: each IG cell group communicated with the current processor
        } // end for: each proc communicated with
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::handle_requested_IG_cell_ID_answers( 
            Cell< Cell< moris_index > > const& aNotOwnedIgCellGroups, 
            Cell< Matrix< IdMat > > const&     aReceivedFirstIgCellIdsInCellGroups )
    {
        // get the communication table and map
        Matrix< IdMat > tCommTable     = this->get_communication_table();
        uint            tCommTableSize = tCommTable.numel();

        // answers received from each proc
        for ( uint iProcInCommTable = 0; iProcInCommTable < tCommTableSize; iProcInCommTable++ )
        {
            // get the number of CMs for which answers were communicated for the current proc
            uint tNumIgCellGroupsCommunicatedWithProc = aReceivedFirstIgCellIdsInCellGroups( iProcInCommTable ).numel();

            // check that all requests have been answered
            MORIS_ASSERT( tNumIgCellGroupsCommunicatedWithProc == aNotOwnedIgCellGroups( iProcInCommTable ).size(),
                    "Cut_Integration_Mesh::handle_requested_IG_cell_ID_answers() - "
                    "Proc #%i: Number of IG cell groups reportedly owned by proc %i is %i, but only answers for %i IG cell groups were received.",
                    par_rank(),
                    tCommTable( iProcInCommTable ),
                    aNotOwnedIgCellGroups( iProcInCommTable ).size(),
                    tNumIgCellGroupsCommunicatedWithProc );

            // go through and answers for all CMs requested by the current proc
            for ( uint iIgCellGroup = 0; iIgCellGroup < tNumIgCellGroupsCommunicatedWithProc; iIgCellGroup++ )
            {
                // get the ID that should be assigned to the first IG cell within this IG cell group
                moris_id tIgCellId = aReceivedFirstIgCellIdsInCellGroups( iProcInCommTable )( iIgCellGroup );

                // get the index of the current not-owned IG cell group
                moris_index tIgCellGroupIndexNotOwned = aNotOwnedIgCellGroups( iProcInCommTable )( iIgCellGroup );

                // if this is an empty cell group just skip it
                if ( tIgCellId == MORIS_ID_MAX )
                {
                    // check that this is indeed an empty IG cell group
                    MORIS_ASSERT(
                            mIntegrationCellGroups( tIgCellGroupIndexNotOwned )->mIgCellGroup.size() == 0,
                            "Cut_Integration_Mesh::handle_requested_IG_cell_ID_answers() - "
                            "MORIS_ID_MAX returned as ID for a IG cell group that is not empty." );

                    // skip to next IG cell group
                    continue;
                }

                // get access to the IG cell group
                std::shared_ptr< IG_Cell_Group > tCellGroup = this->get_ig_cell_group( tIgCellGroupIndexNotOwned );

                // subsequently assign indices to all
                for ( uint iCell = 0; iCell < tCellGroup->mIgCellGroup.size(); iCell++ )
                {
                    // get the index of the current IG cell in the group
                    moris_index tIgCellIndex = tCellGroup->mIgCellGroup( iCell )->get_index();

                    // get the controlled index of that IG cell
                    moris_index tControlledIndex = this->get_integration_cell_controlled_index( tIgCellIndex );

                    // access this IG cell and set its ID
                    mControlledIgCells( tControlledIndex )->set_id( tIgCellId );

                    // check that this ID has not already been assigned to another IG cell
                    MORIS_ASSERT( mIntegrationCellIdToIndexMap.find( tIgCellId ) == mIntegrationCellIdToIndexMap.end(),
                            "Cut_Integration_Mesh::handle_requested_IG_cell_ID_answers() - "
                            "Proc #%i: IG cell ID %i has already been assigned to another IG cell on this processor.",
                            par_rank(),
                            tIgCellId );

                    // populate the map relating the Proc-global IDs to the proc-local indices
                    mIntegrationCellIdToIndexMap[ tIgCellId ] = tIgCellIndex;

                    // increment the IG cell ID to get the ID for the next IG cell in the IG cell group
                    tIgCellId++;
                }

            }    // end for: loop over IG cell groups for each proc that were communicated

        }    // end for: loop over procs ID answers are received from

    } // end function: Cut_Integration_Mesh::handle_requested_IG_cell_ID_answers()

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::create_base_cell_blocks()
    {
        Tracer      tTracer( "XTK", "Integration_Mesh_Generator", "create_base_cell_blocks", mXTKModel->mVerboseLevel, 1 );
        std::string tBlockName = "bg_cells";

        if ( mBackgroundMesh->get_num_elems() > 0 )
        {
            moris::mtk::Cell const & tCell = mBackgroundMesh->get_mtk_cell( 0 );

            // MORIS_ERROR(  tCell.get_geometry_type() == mtk::Geometry_Type::HEX || tCell.get_geometry_type() == mtk::Geometry_Type::QUAD ,
            //     "Need to abstract by adding get cell topo to cell info class" );

            // decide on cell topology based on number of spatial dimensions
            enum CellTopology tCellTopo = tCell.get_cell_info()->get_cell_topology();

            Cell< moris_index > tBlockSetOrds = this->register_block_set_names( { tBlockName }, tCellTopo );

            this->mBlockSetCellGroup( tBlockSetOrds( 0 ) ) = std::make_shared< IG_Cell_Group >();

            this->mBlockSetCellGroup( tBlockSetOrds( 0 ) )->mIgCellGroup.append( mIntegrationCells );
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Cut_Integration_Mesh::trim_data()
    {
        // trim inner and outer cells
        shrink_to_fit_all( mControlledIgCells );
        shrink_to_fit_all( mIntegrationCells );
        shrink_to_fit_all( mIntegrationCellToCellGroupIndex );
        shrink_to_fit_all( mIntegrationCellBulkPhase );
    }

    // ----------------------------------------------------------------------------------

}    // namespace xtk
