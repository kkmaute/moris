/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Factory.cpp
 *
 */

#include "cl_VIS_Output_Manager.hpp"

#include <set>

#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Side_Cluster_Visualization.hpp"

#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Set_Communicator.hpp"
#include "cl_Tracer.hpp"

extern moris::Comm_Manager gMorisComm;

namespace moris
{
    namespace vis
    {
        //-----------------------------------------------------------------------------------------------------------

        VIS_Factory::VIS_Factory(
                std::shared_ptr< mtk::Mesh_Manager > aMesh,
                uint                                 aMeshPairIndex )
        {
            // get access to the integration and interpolation meshes
            aMesh->get_mesh_pair( aMeshPairIndex, mInterpolationMesh, mIntegrationMesh );
        }

        //-----------------------------------------------------------------------------------------------------------

        VIS_Factory::~VIS_Factory()
        {
            // make sure there is no memory leaking
            MORIS_ERROR( mVisMesh == nullptr,
                    "~VIS_Factory() - "
                    "The constructed VIS mesh has not been handed off to a new owner. Destructing this factory will cause memory to leak." );
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::initialize( moris::vis::Output_Data const & aOutputData )
        {
            // determine the type of mesh to be written
            switch ( aOutputData.mMeshType )
            {
                case ( vis::VIS_Mesh_Type::STANDARD ):
                    mOnlyPrimaryCells           = true;
                    mConstructDiscontinuousMesh = false;
                    break;

                case ( vis::VIS_Mesh_Type::STANDARD_WITH_OVERLAP ):
                    mOnlyPrimaryCells           = false;
                    mConstructDiscontinuousMesh = false;
                    break;

                case ( vis::VIS_Mesh_Type::FULL_DISCONTINUOUS ):
                    mOnlyPrimaryCells           = true;
                    mConstructDiscontinuousMesh = true;
                    break;

                case ( vis::VIS_Mesh_Type::FULL_DISCONTINUOUS_WITH_OVERLAP ):
                    mOnlyPrimaryCells           = false;
                    mConstructDiscontinuousMesh = true;
                    break;

                default:
                    MORIS_ERROR( false, "VIS_Factory::initialize() - Mesh type not specified or unknown." );
                    break;
            }

            // create a VIS mesh
            mVisMesh = new vis::Visualization_Mesh( mOnlyPrimaryCells );

            // Get set names
            mAllRequestedSetNames = aOutputData.mSetNames;

            // get the total number of requested sets
            uint tNumRequestedSets = mAllRequestedSetNames.size();

            // reserve memory for cells
            mRequestedBlockSetNames.reserve( tNumRequestedSets );
            mRequestedSideSetNames.reserve( tNumRequestedSets );
            mRequestedDoubleSideSetNames.reserve( tNumRequestedSets );
            mFemBlockSets.reserve( tNumRequestedSets );
            mFemSideSets.reserve( tNumRequestedSets );
            mFemDoubleSideSets.reserve( tNumRequestedSets );

            // figure out which sets are block, side, and dbl-side sets
            uint tNumBlockSets   = 0;
            uint tNumSideSets    = 0;
            uint tNumDblSideSets = 0;

            // count number of block sets
            for ( uint iSet = 0; iSet < tNumRequestedSets; iSet++ )
            {
                // get the current set's name
                std::string tSetName = mAllRequestedSetNames( iSet );

                // get the MTK set that corresponds to the current set name
                mtk::Set* tMtkMeshSet = mIntegrationMesh->get_set_by_name( tSetName );

                // sort sets by type
                switch ( tMtkMeshSet->get_set_type() )
                {
                    case mtk::SetType::BULK:
                    {
                        // count number of block sets
                        tNumBlockSets++;

                        // list the block set in the set lists
                        mRequestedBlockSetNames.push_back( tSetName );
                        mFemBlockSets.push_back( tMtkMeshSet );

                        break;
                    }
                    case mtk::SetType::SIDESET:
                    {
                        tNumSideSets++;
                        mRequestedSideSetNames.push_back( tSetName );
                        mFemSideSets.push_back( tMtkMeshSet );
                        break;
                    }
                    case mtk::SetType::DOUBLE_SIDED_SIDESET:
                    {
                        tNumDblSideSets++;
                        mRequestedDoubleSideSetNames.push_back( tSetName );
                        mFemDoubleSideSets.push_back( tMtkMeshSet );
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false,
                                "VIS_Factory::initialize() - "
                                "Mesh set with name '%s' is neither a block, side, nor dbl. side set. Unable to handle case.",
                                mAllRequestedSetNames( iSet ).c_str() );
                    }
                }    // end switch: type of mtk-mesh
            }        // end for: each requested set

            // copy set names into VIS mesh
            mVisMesh->mBlockSetNames      = mRequestedBlockSetNames;
            mVisMesh->mSideSetNames       = mRequestedSideSetNames;
            mVisMesh->mDoubleSideSetNames = mRequestedDoubleSideSetNames;
            mVisMesh->mAllSetNames.resize( tNumBlockSets + tNumSideSets + tNumDblSideSets );

            // assemble list of all set names and corresponding map
            moris_index tSetIndex = 0;
            for ( const std::string& iBlockSetName : mVisMesh->mBlockSetNames )
            {
                // make sure there are no repeated set names
                MORIS_ERROR( !mVisMesh->mSetNameToIndexMap.key_exists( iBlockSetName ),
                        "VIS::Visualization_Mesh::initialize() - "
                        "Set with name '%s' has been found twice in list of sets in VIS mesh. "
                        "Make sure that mesh sets are only listed once in the VIS parameter list.",
                        iBlockSetName.c_str() );

                mVisMesh->mAllSetNames( tSetIndex )           = iBlockSetName;
                mVisMesh->mSetNameToIndexMap[ iBlockSetName ] = tSetIndex;
                tSetIndex++;
            }
            for ( const std::string& iSideSetName : mVisMesh->mSideSetNames )
            {
                // make sure there are no repeated set names
                MORIS_ERROR( !mVisMesh->mSetNameToIndexMap.key_exists( iSideSetName ),
                        "VIS::Visualization_Mesh::initialize() - "
                        "Set with name '%s' has been found twice in list of sets in VIS mesh. "
                        "Make sure that mesh sets are only listed once in the VIS parameter list.",
                        iSideSetName.c_str() );

                mVisMesh->mAllSetNames( tSetIndex )          = iSideSetName;
                mVisMesh->mSetNameToIndexMap[ iSideSetName ] = tSetIndex;
                tSetIndex++;
            }
            for ( const std::string& iDblSideSetName : mVisMesh->mDoubleSideSetNames )
            {
                // make sure there are no repeated set names
                MORIS_ERROR( !mVisMesh->mSetNameToIndexMap.key_exists( iDblSideSetName ),
                        "VIS::Visualization_Mesh::initialize() - "
                        "Set with name '%s' has been found twice in list of sets in VIS mesh. "
                        "Make sure that mesh sets are only listed once in the VIS parameter list.",
                        iDblSideSetName.c_str() );

                mVisMesh->mAllSetNames( tSetIndex )             = iDblSideSetName;
                mVisMesh->mSetNameToIndexMap[ iDblSideSetName ] = tSetIndex;
                tSetIndex++;
            }

            // initialize lists in the VIS mesh
            mVisMesh->mListOfBlocks.resize( tNumBlockSets );
            mVisMesh->mListOfSideSets.resize( tNumSideSets );
            mVisMesh->mListOfDoubleSideSets.resize( tNumDblSideSets );
            mVisMesh->mListOfAllSets.resize( 0 );    // this will be sized correctly in the finalize step of the VIS mesh

            mVisMesh->mClustersOnBlockSets.resize( tNumBlockSets );
            mVisMesh->mClustersOnSideSets.resize( tNumSideSets );
            mVisMesh->mClustersOnDoubleSideSets.resize( tNumDblSideSets );
            mVisMesh->mLeaderSideClusters.resize( tNumDblSideSets );
            mVisMesh->mFollowerSideClusters.resize( tNumDblSideSets );

            // initialize  maps
            uint tNumFemCells = mIntegrationMesh->get_num_entities( mtk::EntityRank::ELEMENT ) + 1;

            mPrimaryFemCellIndexToVisCellIndex.resize( tNumFemCells, -1 );
            mPrimaryFemCellIndexToBlockIndex.resize( tNumFemCells, -1 );

            mBlockAndFemCellIndexToVisCellIndex.resize( tNumBlockSets );
            mBlockAndFemCellIndexToVisVertexIndices.resize( tNumBlockSets );

            for ( uint iBlockSet = 0; iBlockSet < tNumBlockSets; iBlockSet++ )
            {
                mBlockAndFemCellIndexToVisCellIndex( iBlockSet ).resize( tNumFemCells, -1 );
                mBlockAndFemCellIndexToVisVertexIndices( iBlockSet ).resize( tNumFemCells );
            }

        }    // end function: VIS_Factory::initialize()

        //-----------------------------------------------------------------------------------------------------------

        mtk::Mesh*
        VIS_Factory::create_visualization_mesh( moris::vis::Output_Data& aOutputData )
        {
            // time/log this operation
            Tracer tTracer( "VIS", "Factory", "Create VIS Mesh" );

            // initialize member data considering the requested output data
            this->initialize( aOutputData );

            // create the vertex objects
            this->create_visualization_vertices();

            // create the cells
            this->create_visualization_cells();

            // group cells into clusters
            this->create_visualization_clusters();

            // create side clusters
            this->create_visualization_side_clusters();

            // create the double sided side clusters
            this->create_visualization_double_side_clusters();

            // create block sets from the cell clusters
            this->create_visualization_blocks();

            // create side sets from the side clusters
            this->create_visualization_side_sets();

            // create the double sided side sets
            this->create_visualization_double_side_sets();

            // finalize data in the VIS mesh before handing over ownership
            mVisMesh->finalize();

            // return pointer to VIS mesh to new owning entity
            return this->hand_off_VIS_mesh();
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::create_visualization_vertices()
        {
            // depending on whether a standard or discontinuous output mesh has been requested, the logic for creating vertices changes
            if ( mConstructDiscontinuousMesh )
            {
                this->create_visualization_vertices_full_discontinuous();
            }
            else
            {
                this->create_visualization_vertices_standard();
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::create_visualization_vertices_standard()
        {
            // log/trace this operation
            Tracer tTracer( "VIS", "Factory", "Create Vertices (standard mesh)" );

            // NOTE: The "FEM vertices" are not unzipped.
            // NOTE: Each set needs to get its own vertices such that jumps along interfaces can be resolved in the output mesh.

            // sum up number of vertices on all block sets
            uint tNumVerticesInProcLocalVisMesh = 0;
            for ( uint iBlockSet = 0; iBlockSet < mRequestedBlockSetNames.size(); iBlockSet++ )
            {
                // sum up number of vertices on all block sets
                tNumVerticesInProcLocalVisMesh += mFemBlockSets( iBlockSet )->get_num_vertices_on_set( mOnlyPrimaryCells );
            }

            // reserve memory for this many vertices
            mVisMesh->mVertices.resize( tNumVerticesInProcLocalVisMesh );

            // get the max number of vertices that the VIS mesh could have
            uint tNumVerticesInFemMesh = mIntegrationMesh->get_num_entities( mtk::EntityRank::NODE ) + 1;

            // reserve enough memory in map to accommodate entries for all vertices
            mBlockAndFemCellIndexToVisVertexIndices.reserve( 10 * mRequestedBlockSetNames.size() * tNumVerticesInFemMesh );

            // update the id for each processor
            moris_id tNextFreeVisVertexID = get_processor_offset( tNumVerticesInProcLocalVisMesh );

            // initialize vertex index counter
            moris_index tVisVertexIndexCounter = 0;

            // loop over all requested block sets
            for ( uint iBlockSet = 0; iBlockSet < mRequestedBlockSetNames.size(); iBlockSet++ )
            {
                // get the number of vertices on the set, ignore if not block
                uint tNumMtkVerticesOnBlock = mFemBlockSets( iBlockSet )->get_num_vertices_on_set( mOnlyPrimaryCells );

                // skip set if it is empty
                if ( tNumMtkVerticesOnBlock == 0 )
                {
                    continue;
                }

                // ------------------
                // construct VIS vertex objects

                // initialize map relating the FEM/MTK vertex index to the VIS vertex index for the current block set
                Vector< moris_index > tFemVertexIndexToVisVertexIndex( tNumVerticesInFemMesh, -1 );

                // get vertex indices on block set
                moris::Matrix< DDSMat > tFemVertexIndicesInFemBlock = mFemBlockSets( iBlockSet )->get_ig_vertices_inds_on_block( mOnlyPrimaryCells );

                // Loop over all vertices on this set and create new vis vertices
                for ( uint iVertexOnSet = 0; iVertexOnSet < tNumMtkVerticesOnBlock; iVertexOnSet++ )
                {
                    // get the index of the current vertex on the current block
                    moris_index tFemVertexIndex = tFemVertexIndicesInFemBlock( iVertexOnSet );

                    // create vis vertex and with id and index
                    mVisMesh->mVertices( tVisVertexIndexCounter ) = new Vertex_Visualization(
                            tNextFreeVisVertexID,
                            tVisVertexIndexCounter,
                            &mIntegrationMesh->get_mtk_vertex( tFemVertexIndex ) );

                    // relate the FEM/MTK vertex index to the VIS vertex index for the current block set
                    tFemVertexIndexToVisVertexIndex( tFemVertexIndex ) = tVisVertexIndexCounter;

                    // update the index and ID counters for the next vertex
                    tNextFreeVisVertexID++;
                    tVisVertexIndexCounter++;
                }

                // ------------------
                // construct map relating the the VIS vertices with the FEM cells

                // get IG cell indices on set
                moris::Matrix< DDSMat > tFemCellIndicesInFemSet = mFemBlockSets( iBlockSet )->get_cell_inds_on_block( mOnlyPrimaryCells );

                // Loop over all cells on this set and create new vis cells
                for ( uint iCellOnBlockSet = 0; iCellOnBlockSet < tFemCellIndicesInFemSet.numel(); iCellOnBlockSet++ )
                {
                    // get the current cell's index in FEM/MTK
                    moris_index tFemCellIndex = tFemCellIndicesInFemSet( iCellOnBlockSet );

                    // get list of integration vertex indices for the current cell
                    Matrix< IndexMat > tFemVertexIndicesOnCell = mIntegrationMesh->get_mtk_cell( tFemCellIndex ).get_vertex_inds();
                    uint               tNumVertsOnCell         = tFemVertexIndicesOnCell.numel();

                    // initialize map
                    mBlockAndFemCellIndexToVisVertexIndices( iBlockSet )( tFemCellIndex ).resize( tNumVertsOnCell );

                    // relate each VIS vertex with the FEM cells on the blocks they were constructed from
                    for ( uint iVertexOnCell = 0; iVertexOnCell < tNumVertsOnCell; iVertexOnCell++ )
                    {
                        // get indices
                        moris_index tFemVertexIndex = tFemVertexIndicesOnCell( iVertexOnCell );
                        moris_index tVisVertexIndex = tFemVertexIndexToVisVertexIndex( tFemVertexIndex );

                        // fill map
                        mBlockAndFemCellIndexToVisVertexIndices( iBlockSet )( tFemCellIndex )( iVertexOnCell ) = tVisVertexIndex;
                    }

                }    // end for: each IG cell on current block set

            }    // end for: each block set

        }    // end function: VIS_Factory::create_visualization_vertices_standard()

        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::create_visualization_vertices_full_discontinuous()
        {
            // log/trace this operation
            Tracer tTracer( "VIS", "Factory", "Create Vertices (discontinuous mesh)" );

            // sum up number of vertices on all block sets
            uint tNumVerticesInProcLocalVisMesh = 0;
            for ( uint iBlockSet = 0; iBlockSet < mRequestedBlockSetNames.size(); iBlockSet++ )
            {
                // get the clusters in the current block
                Vector< mtk::Cluster const * > tClustersOnSet = mFemBlockSets( iBlockSet )->get_clusters_on_set();

                // count up the number of VIS vertices on this block by counting the number of vertices in each cluster
                // disregarding FEM vertices which may be shared between clusters
                for ( mtk::Cluster const * iCluster : tClustersOnSet )
                {
                    tNumVerticesInProcLocalVisMesh += iCluster->get_num_vertices_in_cluster();
                }
            }

            // reserve memory for this many vertices
            mVisMesh->mVertices.resize( tNumVerticesInProcLocalVisMesh );

            // get the max number of vertices that the VIS mesh could have
            uint tNumVerticesInFemMesh = mIntegrationMesh->get_num_entities( mtk::EntityRank::NODE ) + 1;

            // reserve enough memory in map to accommodate entries for all vertices
            mBlockAndFemCellIndexToVisVertexIndices.reserve( 10 * mRequestedBlockSetNames.size() * tNumVerticesInFemMesh );

            // update the id for each processor
            moris_id tNextFreeVisVertexID = get_processor_offset( tNumVerticesInProcLocalVisMesh );

            // initialize vertex index counter
            moris_index tVisVertexIndexCounter = 0;

            // loop over all requested block sets
            for ( uint iBlockSet = 0; iBlockSet < mRequestedBlockSetNames.size(); iBlockSet++ )
            {
                // get the number of vertices on the set, ignore if not block
                uint tNumMtkVerticesOnBlock = mFemBlockSets( iBlockSet )->get_num_vertices_on_set( mOnlyPrimaryCells );

                // skip set if it is empty
                if ( tNumMtkVerticesOnBlock == 0 )
                {
                    continue;
                }

                // ------------------
                // initialize map relating the the VIS vertices with the FEM cells

                // get IG cell indices on set
                moris::Matrix< DDSMat > tFemCellIndicesInFemSet = mFemBlockSets( iBlockSet )->get_cell_inds_on_block( mOnlyPrimaryCells );

                // Loop over all cells on this set and create new vis cells
                for ( uint iCellOnBlockSet = 0; iCellOnBlockSet < tFemCellIndicesInFemSet.numel(); iCellOnBlockSet++ )
                {
                    // get the current cell's index in FEM/MTK
                    moris_index tFemCellIndex = tFemCellIndicesInFemSet( iCellOnBlockSet );

                    // get list of integration vertex indices for the current cell
                    Matrix< IndexMat > tFemVertexIndicesOnCell = mIntegrationMesh->get_mtk_cell( tFemCellIndex ).get_vertex_inds();
                    uint               tNumVertsOnCell         = tFemVertexIndicesOnCell.numel();

                    // initialize map
                    mBlockAndFemCellIndexToVisVertexIndices( iBlockSet )( tFemCellIndex ).resize( tNumVertsOnCell, -1 );

                }    // end for: each IG cell on current block set

                // ------------------
                // construct VIS vertex objects

                // get the clusters in the current block
                Vector< mtk::Cluster const * > tClustersOnSet = mFemBlockSets( iBlockSet )->get_clusters_on_set();

                // create vertices for each cluster separately disregarding FEM vertices which may be shared between clusters
                for ( mtk::Cluster const * iCluster : tClustersOnSet )
                {
                    // initialize map relating FEM/MTK vertices to the VIS vertices to be constructed
                    map< moris_index, moris_index > tFemVertexIndexToVisVertexIndex;

                    // ------------------
                    // construct VIS vertex objects

                    // get the FEM/MTK vertices on the current cluster
                    Vector< mtk::Vertex const * > tFemVerticesInCluster = iCluster->get_vertices_in_cluster();

                    // create a VIS
                    for ( mtk::Vertex const * iFemVertexInCluster : tFemVerticesInCluster )
                    {
                        // get the index of the fem vertex on the current cluster
                        moris_index tFemVertexIndex = iFemVertexInCluster->get_index();

                        // create vis vertex and with id and index
                        mVisMesh->mVertices( tVisVertexIndexCounter ) = new Vertex_Visualization(
                                tNextFreeVisVertexID,
                                tVisVertexIndexCounter,
                                &mIntegrationMesh->get_mtk_vertex( tFemVertexIndex ) );

                        // record the VIS vertex constructed from the FEM vertex
                        tFemVertexIndexToVisVertexIndex[ tFemVertexIndex ] = tVisVertexIndexCounter;

                        // update the index and ID counters for the next vertex
                        tNextFreeVisVertexID++;
                        tVisVertexIndexCounter++;
                    }

                    // ------------------
                    // fill map relating the the VIS vertices with the FEM cells

                    // get the FEM/MTK cells on the current cluster
                    Vector< mtk::Cell const * > tFemCellsInCluster = iCluster->get_primary_cells_in_cluster();

                    // add or leave out void cells depending on the mesh type requested
                    if ( !mOnlyPrimaryCells )
                    {
                        tFemCellsInCluster.append( iCluster->get_void_cells_in_cluster() );
                    }

                    // go over each cell on the current cluster and relate the FEM/MTK vertices local to the element to the VIS vertices constructed from them
                    for ( auto iUsedFemCell : tFemCellsInCluster )
                    {
                        // get the index of the current FEM cell
                        moris_index tFemCellIndex = iUsedFemCell->get_index();

                        // get the indices of the FEM vertices attached to the cell
                        Matrix< IndexMat > tFemVertexIndicesOnCell = iUsedFemCell->get_vertex_inds();

                        // relate each VIS vertex back to the FEM vertex it was constructed from through the map
                        for ( uint iVertOnCell = 0; iVertOnCell < tFemVertexIndicesOnCell.numel(); iVertOnCell++ )
                        {
                            moris_index tFemVertexIndex                                                          = tFemVertexIndicesOnCell( iVertOnCell );
                            moris_index tVisVertexIndex                                                          = tFemVertexIndexToVisVertexIndex.find( tFemVertexIndex );
                            mBlockAndFemCellIndexToVisVertexIndices( iBlockSet )( tFemCellIndex )( iVertOnCell ) = tVisVertexIndex;
                        }
                    }

                }    // end for: each cluster on current block set

            }    // end for: each block set

        }    // end function: VIS_Factory::create_visualization_vertices_full_discontinuous()

        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::create_visualization_cells()
        {
            // time and log this function
            Tracer tTracer( "VIS", "Factory", "Create IG Cells" );

            // sum up number of cells on all sets
            uint tNumCellsInProcLocalVisMesh = 0;
            for ( uint iBlockSet = 0; iBlockSet < mRequestedBlockSetNames.size(); iBlockSet++ )
            {
                tNumCellsInProcLocalVisMesh += mFemBlockSets( iBlockSet )->get_num_cells_on_set( mOnlyPrimaryCells );
            }

            // reserve memory for this many cells
            mVisMesh->mCells.resize( tNumCellsInProcLocalVisMesh );

            // update the id for each processor
            moris_id tNextFreeCellID = get_processor_offset( tNumCellsInProcLocalVisMesh );

            // initialize cell index counter
            moris_index tCellIndexCounter = 0;

            // loop over all requested block sets
            for ( uint iBlockSet = 0; iBlockSet < mRequestedBlockSetNames.size(); iBlockSet++ )
            {
                // get number of cells on set
                uint tNumCellsOnSet = mFemBlockSets( iBlockSet )->get_num_cells_on_set( mOnlyPrimaryCells );

                // ignore the subsequent steps for empty sets and jump to next set
                if ( tNumCellsOnSet == 0 )
                {
                    continue;
                }

                // get cell indices on set
                moris::Matrix< DDSMat > tCellIndicesInFemSet = mFemBlockSets( iBlockSet )->get_cell_inds_on_block( mOnlyPrimaryCells );

                // Loop over all cells on this set and create new vis cells
                for ( uint iCellOnBlockSet = 0; iCellOnBlockSet < tNumCellsOnSet; iCellOnBlockSet++ )
                {
                    // get the current cells index in MTK
                    moris_index tFemCellIndex = tCellIndicesInFemSet( iCellOnBlockSet );

                    // get list of VIS vertex indices for the current cell
                    const Vector< moris_index > tVisVertexIndicesOnCell = mBlockAndFemCellIndexToVisVertexIndices( iBlockSet )( tFemCellIndex );
                    uint                      tNumVertsOnCell         = tVisVertexIndicesOnCell.size();

                    // Create List of VIS vertex pointers for this cell
                    Vector< mtk::Vertex* > tVisCellVertices( tNumVertsOnCell, nullptr );

                    // loop over integration vertices and get the corresponding vis vertices for this set
                    for ( uint iVertOnCell = 0; iVertOnCell < tNumVertsOnCell; iVertOnCell++ )
                    {
                        // get the index of the VIS vertex on the cell to be constructed
                        moris_index tVisVertIndex = tVisVertexIndicesOnCell( iVertOnCell );

                        // store pointer to the corresponding VIS vertices
                        tVisCellVertices( iVertOnCell ) = mVisMesh->mVertices( tVisVertIndex );
                    }

                    // create vis cells and renumber id and index
                    mVisMesh->mCells( tCellIndexCounter ) = new Cell_Visualization(
                            tNextFreeCellID,
                            tCellIndexCounter,
                            tVisCellVertices,
                            &mIntegrationMesh->get_mtk_cell( tFemCellIndex ) );

                    // relate the index of the new VIS cell back to the FEM cell it was constructed from for the current block set
                    mBlockAndFemCellIndexToVisCellIndex( iBlockSet )( tFemCellIndex ) = tCellIndexCounter;

                    // update the index and ID counters for the next vertex
                    tNextFreeCellID++;
                    tCellIndexCounter++;

                }    // end for: each cell on the current block set'


                // ------------------
                // fill the maps about the primary cells
                // this needs to be done in a separate step as for the case of an overlapping mesh the primary cells first need to be sorted out

                // get cell indices of primary cells on set
                Matrix< DDSMat > tPrimaryCellIndicesInFemSet = mFemBlockSets( iBlockSet )->get_cell_inds_on_block( true );

                // relate each primary cell in block back to the FEM/MTK cell it was constructed from
                for ( uint iPrimaryCellInBlock = 0; iPrimaryCellInBlock < tPrimaryCellIndicesInFemSet.numel(); iPrimaryCellInBlock++ )
                {
                    // get the current FEM/MTK cell's index
                    moris_index tFemCellIndex = tPrimaryCellIndicesInFemSet( iPrimaryCellInBlock );

                    // get the corresponding
                    moris_index tVisCellIndex = mBlockAndFemCellIndexToVisCellIndex( iBlockSet )( tFemCellIndex );

                    // there should never be more than one primary VIS cell be associated with a given MTK/FEM cell
                    MORIS_ASSERT(
                            mPrimaryFemCellIndexToVisCellIndex( tFemCellIndex ) == -1,
                            "VIS_Factory::create_visualization_cells() - "
                            "Assigning VIS cell with FEM/MTK cell #%i, but it already has a primary VIS cell associated with it.",
                            tFemCellIndex );

                    // fill the maps
                    mPrimaryFemCellIndexToVisCellIndex( tFemCellIndex ) = tVisCellIndex;
                    mPrimaryFemCellIndexToBlockIndex( tFemCellIndex )   = iBlockSet;
                }

            }    // end for: each block set

        }    // end function: VIS_Factory::create_visualization_cells()

        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::create_visualization_clusters()
        {
            // log and time this function
            Tracer tTracer( "VIS", "Factory", "Create Cell Clusters" );

            // loop over requested block sets
            for ( uint iBlockSet = 0; iBlockSet < mRequestedBlockSetNames.size(); iBlockSet++ )
            {
                // get number of clusters on set
                uint tNumClustersOnSet = mFemBlockSets( iBlockSet )->get_num_clusters_on_set();

                // ignore the subsequent steps for empty sets, and jump to next set
                if ( tNumClustersOnSet == 0 )
                {
                    continue;
                }

                // get list of old clusters on set
                Vector< mtk::Cluster const * > tClustersOnFemSet = mFemBlockSets( iBlockSet )->get_clusters_on_set();

                // resize list of clusters for this set
                mVisMesh->mClustersOnBlockSets( iBlockSet ).resize( tNumClustersOnSet, nullptr );

                // loop over clusters on set
                for ( uint iClusterOnSet = 0; iClusterOnSet < tNumClustersOnSet; iClusterOnSet++ )
                {
                    // create the VIS cluster object to feed the information to
                    Cell_Cluster_Visualization* tVisCellCluster = new Cell_Cluster_Visualization;

                    // mark as non-trivial if old cluster was not trivial
                    if ( !tClustersOnFemSet( iClusterOnSet )->is_trivial() )
                    {
                        tVisCellCluster->mark_as_nontrivial();
                    }

                    // ------------------------
                    // add cells to cluster

                    // get primary/void cells on old cluster
                    const Vector< moris::mtk::Cell const * >& tPrimaryFemCells = tClustersOnFemSet( iClusterOnSet )->get_primary_cells_in_cluster();
                    const Vector< moris::mtk::Cell const * >& tVoidFemCells    = tClustersOnFemSet( iClusterOnSet )->get_void_cells_in_cluster();

                    // resize primary/void cell list for VIS cluster
                    // (these will be filled with VIS cells that have been constructed from the corresponding mtk/fem cells)
                    Vector< moris::mtk::Cell const * > tClusterPrimaryVisCells( tPrimaryFemCells.size(), nullptr );
                    Vector< moris::mtk::Cell const * > tClusterVoidVisCells( tVoidFemCells.size(), nullptr );

                    // find VIS primary cells corresponding to mtk/fem cells
                    for ( uint iPrimaryCell = 0; iPrimaryCell < tPrimaryFemCells.size(); iPrimaryCell++ )
                    {
                        moris_index tPrimaryFemCellIndex        = tPrimaryFemCells( iPrimaryCell )->get_index();
                        moris_index tPrimaryVisCellIndex        = mBlockAndFemCellIndexToVisCellIndex( iBlockSet )( tPrimaryFemCellIndex );
                        tClusterPrimaryVisCells( iPrimaryCell ) = mVisMesh->mCells( tPrimaryVisCellIndex );
                    }

                    // find VIS void cells corresponding to mtk/fem cells
                    if ( !mOnlyPrimaryCells )
                    {
                        for ( uint iVoidCell = 0; iVoidCell < tVoidFemCells.size(); iVoidCell++ )
                        {
                            moris_index tVoidFemCellIndex     = tVoidFemCells( iVoidCell )->get_index();
                            moris_index tVoidVisCellIndex     = mBlockAndFemCellIndexToVisCellIndex( iBlockSet )( tVoidFemCellIndex );
                            tClusterVoidVisCells( iVoidCell ) = mVisMesh->mCells( tVoidVisCellIndex );
                        }
                    }

                    // add vis primary cells to vis cluster
                    tVisCellCluster->add_primary_integration_cell( tClusterPrimaryVisCells );

                    // add interpolation cell to vis cluster
                    tVisCellCluster->set_interpolation_cell( &tClustersOnFemSet( iClusterOnSet )->get_interpolation_cell() );

                    // add void cells to cluster if requested by mesh type
                    if ( tVoidFemCells.size() > 0 && !mOnlyPrimaryCells )
                    {
                        tVisCellCluster->add_void_integration_cell( tClusterVoidVisCells );
                    }

                    // ------------------------
                    // add vertices to cluster

                    // get vertices from the mtk/fem cluster
                    Vector< moris::mtk::Vertex const * > tFemVertices = tClustersOnFemSet( iClusterOnSet )->get_vertices_in_cluster();

                    // construct map relating fem vertex indices to their respective position in the list of vertices on the cluster
                    map< moris_index, moris_index > tFemVertexIndexToPosInFemClusterMap;
                    for ( uint iVertInFemCluster = 0; iVertInFemCluster < tFemVertices.size(); iVertInFemCluster++ )
                    {
                        moris_index tFemVertexIndex                            = tFemVertices( iVertInFemCluster )->get_index();
                        tFemVertexIndexToPosInFemClusterMap[ tFemVertexIndex ] = iVertInFemCluster;
                    }

                    // initialize map identifying which VIS vertices are in the VIS cluster to be constructed and at which position
                    map< moris_index, moris_index > tVisVertexIndexToPosInVisClusterMap;
                    uint                            tNextVertexInVisClusterPosition = 0;

                    // initialize a list of all VIS vertices used by the cluster
                    Vector< mtk::Vertex const * > tClusterVisVertices;
                    tClusterVisVertices.reserve( tFemVertices.size() );

                    // compile a list of IG cells used by the newly constructed VIS cluster
                    Vector< const moris::mtk::Cell* > tUsedClusterFemCells = tPrimaryFemCells;

                    // for overlapping meshes consider the void as well
                    if ( !mOnlyPrimaryCells )
                    {
                        tUsedClusterFemCells.append( tVoidFemCells );
                    }

                    // go over these cells and compile a list of vertices attached to them
                    for ( uint iCellInCluster = 0; iCellInCluster < tUsedClusterFemCells.size(); iCellInCluster++ )
                    {
                        // get the index of the currently considered FEM/MTK cell
                        moris_index tFemCellIndex = tUsedClusterFemCells( iCellInCluster )->get_index();

                        // access the list of VIS vertices attached to the cell
                        Vector< moris_index > tVisVertexIndicesOnCell = mBlockAndFemCellIndexToVisVertexIndices( iBlockSet )( tFemCellIndex );

                        MORIS_ASSERT( tVisVertexIndicesOnCell.size() > 0,
                                "VIS_Factory::create_visualization_clusters() - Map doesn't show any vertices connected to IG cell." );

                        // go over the VIS vertices on the cell and add them to the list of cluster cells if they haven't been already
                        for ( uint iVertOnCell = 0; iVertOnCell < tVisVertexIndicesOnCell.size(); iVertOnCell++ )
                        {
                            // get the index of the current VIS vertex
                            moris_index tVisVertexIndex = tVisVertexIndicesOnCell( iVertOnCell );

                            // check if this VIS vertex has already been added to the cluster, if not add it to the cluster
                            if ( !tVisVertexIndexToPosInVisClusterMap.key_exists( tVisVertexIndex ) )
                            {
                                // add VIS vertex to cluster
                                tVisVertexIndexToPosInVisClusterMap[ tVisVertexIndex ] = tNextVertexInVisClusterPosition;
                                tClusterVisVertices.push_back( mVisMesh->mVertices( tVisVertexIndex ) );
                                tNextVertexInVisClusterPosition++;
                            }
                        }
                    }    // end for: each used IG cell in the current cluster

                    // store pointer list of VIS vertices in the VIS cluster being constructed
                    tVisCellCluster->add_vertex_to_cluster( tClusterVisVertices );

                    // ------------------------
                    // add vertex coordinates

                    // get the coordinates for all vertices in the corresponding FEM cluster
                    Matrix< DDRMat > tLocalCoordsOfAllFemVertsOnCluster =
                            tClustersOnFemSet( iClusterOnSet )->get_vertices_local_coordinates_wrt_interp_cell();

                    // get the number of dimensions in parametric space
                    uint tNumParamDims = tLocalCoordsOfAllFemVertsOnCluster.n_cols();

                    // initialize matrix containing all local coords for the vertices in the VIS cluster
                    Matrix< DDRMat > tVisClusterVerticesLocalCoords( tNextVertexInVisClusterPosition, tNumParamDims, 0.0 );

                    // go through the vertices in the VIS cluster, find the corresponding FEM vertices, and copy over the local coordinates
                    for ( uint iVertInVisCluster = 0; iVertInVisCluster < tNextVertexInVisClusterPosition; iVertInVisCluster++ )
                    {
                        // get the corresponding FEM vertex's index
                        moris_index tFemVertexIndex = tClusterVisVertices( iVertInVisCluster )->get_base_vertex()->get_index();

                        // get the FEM vertex's position in the list of cluster vertices
                        moris_index tVertPosInFemCluster = tFemVertexIndexToPosInFemClusterMap.find( tFemVertexIndex );

                        // copy over coordinates
                        tVisClusterVerticesLocalCoords.set_row(
                                iVertInVisCluster,
                                tLocalCoordsOfAllFemVertsOnCluster.get_row( tVertPosInFemCluster ) );
                    }

                    // add local coordinates to vis cluster
                    tVisCellCluster->add_vertex_local_coordinates_wrt_interp_cell( tVisClusterVerticesLocalCoords );

                    // ------------------------
                    // store the constructed cluster in the VIS-mesh

                    // store away constructed cluster
                    mVisMesh->mClustersOnBlockSets( iBlockSet )( iClusterOnSet ) = tVisCellCluster;

                }    // end for: clusters on set

            }    // end for: each set

        }    // end function: VIS_Factory::create_visualization_clusters()


        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::create_visualization_side_clusters()
        {
            // time and log this function
            Tracer tTracer( "VIS", "Factory", "Create Side Clusters" );

            // loop over requested sets
            for ( uint iSideSet = 0; iSideSet < mRequestedSideSetNames.size(); iSideSet++ )
            {
                // get access to the current set
                mtk::Set* tFemSideSet = mFemSideSets( iSideSet );

                // get number of side clusters on set
                uint tNumSideClustersOnSet = tFemSideSet->get_num_clusters_on_set();

                // ignore the subsequent steps for empty sets, and jump to next set
                if ( tNumSideClustersOnSet == 0 )
                {
                    continue;
                }

                // get list of clusters on existing fem/mtk set
                Vector< mtk::Cluster const * > tSideClustersOnFemSet = tFemSideSet->get_clusters_on_set();

                // resize list of clusters for this set
                mVisMesh->mClustersOnSideSets( iSideSet ).resize( tNumSideClustersOnSet, nullptr );

                // loop over clusters in existing fem/mtk set
                for ( uint iSideClusterOnSet = 0; iSideClusterOnSet < tSideClustersOnFemSet.size(); iSideClusterOnSet++ )
                {
                    // create the VIS side cluster object. the below information will be feed into it
                    vis::Side_Cluster_Visualization* tVisSideCluster = new vis::Side_Cluster_Visualization();

                    // mark as non-trivial if old cluster was not trivial
                    if ( !tSideClustersOnFemSet( iSideClusterOnSet )->is_trivial() )
                    {
                        tVisSideCluster->mark_as_nontrivial();
                    }

                    // ------------------------
                    // IG cells that the side cluster is attached to

                    // get the IG cells the side elements are attached to
                    const Vector< const moris::mtk::Cell* >& tFemCellsInCluster =
                            tSideClustersOnFemSet( iSideClusterOnSet )->get_primary_cells_in_cluster();

                    MORIS_ASSERT( tFemCellsInCluster.size() > 0,
                            "VIS_Factory::create_visualization_side_clusters() - Empty side cluster in FEM mesh. This shouldn't happen." );

                    // initialize list of VIS IG cells the side set elements are attached to
                    Vector< moris::mtk::Cell const * > tSideClusterVisIgCells( tFemCellsInCluster.size(), nullptr );

                    // find the corresponding VIS cells constructed from the original mtk/fem cells and collect these in a list for constructing the VIS side cluster
                    for ( uint iIgCellOnCluster = 0; iIgCellOnCluster < tFemCellsInCluster.size(); iIgCellOnCluster++ )
                    {
                        // get the current Fem cell's index
                        moris_index tFemCellIndex = tFemCellsInCluster( iIgCellOnCluster )->get_index();

                        // get the corresponding VIS cell's index
                        moris_index tVisCellIndex = mPrimaryFemCellIndexToVisCellIndex( tFemCellIndex );

                        // make sure the requested VIS cell actually exists
                        MORIS_ERROR( tVisCellIndex > -1 && tVisCellIndex != MORIS_INDEX_MAX,
                                "VIS_Factory::create_visualization_side_clusters() - "
                                "No VIS cell index for the given FEM cell. "
                                "The current side cluster may be attached to a block set that is not part of the VIS mesh. "
                                "Make sure to include it." );

                        // store the corrsponding VIS cell
                        tSideClusterVisIgCells( iIgCellOnCluster ) = mVisMesh->mCells( tVisCellIndex );
                    }

                    // store list of cells on VIS cluster
                    tVisSideCluster->add_primary_integration_cell( tSideClusterVisIgCells );

                    // add interpolation cell to vis cluster
                    tVisSideCluster->set_interpolation_cell( &tSideClustersOnFemSet( iSideClusterOnSet )->get_interpolation_cell() );

                    // ------------------------
                    // side ordinals of the IG cells that the side cluster is attached to

                    // get the side ordinals of each IG cell that make up the side cluster
                    moris::Matrix< moris::IndexMat > tSideOrdinals = tSideClustersOnFemSet( iSideClusterOnSet )->get_cell_side_ordinals();

                    // store information in side cluster
                    tVisSideCluster->add_integration_cell_side_ordinals( tSideOrdinals );

                    // ------------------------
                    // add vertices to cluster

                    // get vertices from the mtk/fem cluster
                    Vector< moris::mtk::Vertex const * > tFemVertices = tSideClustersOnFemSet( iSideClusterOnSet )->get_vertices_in_cluster();

                    // get the number of vertices on the fem cluster for convenient access
                    uint tNumUsedVerticesInFemCluster = tFemVertices.size();

                    // construct map relating fem vertex indices to their respective position in the list of vertices on the cluster
                    map< moris_index, moris_index > tFemVertexIndexToPosInFemClusterMap;
                    for ( uint iVertInFemCluster = 0; iVertInFemCluster < tNumUsedVerticesInFemCluster; iVertInFemCluster++ )
                    {
                        moris_index tFemVertexIndex                            = tFemVertices( iVertInFemCluster )->get_index();
                        tFemVertexIndexToPosInFemClusterMap[ tFemVertexIndex ] = iVertInFemCluster;
                    }

                    // initialize list of VIS vertices constructed from the corresponding mtk/fem vertices
                    Vector< mtk::Vertex const * > tClusterVisVertices;
                    tClusterVisVertices.reserve( tNumUsedVerticesInFemCluster );

                    // initialize map identifying which VIS vertices are in the newly constructed VIS cluster and at which position
                    map< moris_index, moris_index > tVisVertexIndexToPosInVisClusterMap;
                    uint                            tNextVertexInVisClusterPosition = 0;

                    // go through IG cells the side cluster is attached to
                    for ( uint iIgCellOnCluster = 0; iIgCellOnCluster < tFemCellsInCluster.size(); iIgCellOnCluster++ )
                    {
                        // get the current Fem cell's index
                        moris_index tFemCellIndex = tFemCellsInCluster( iIgCellOnCluster )->get_index();

                        // get the block set index the side cluster is attached to
                        moris_index tPrimaryBlockSetIndex = mPrimaryFemCellIndexToBlockIndex( tFemCellIndex );

                        // get the VIS vertex indices
                        const Vector< moris_index > tVisVertexIndices = mBlockAndFemCellIndexToVisVertexIndices( tPrimaryBlockSetIndex )( tFemCellIndex );

                        // go through the vertices listed in the FEM cluster, find the corresponding VIS vertices, and store them in the above list
                        for ( uint iVertOnCell = 0; iVertOnCell < tVisVertexIndices.size(); iVertOnCell++ )
                        {
                            // get the VIS vertex's index
                            moris_index tVisVertexIndex = tVisVertexIndices( iVertOnCell );

                            // get the VIS vertex
                            mtk::Vertex const * tVisVertex = mVisMesh->mVertices( tVisVertexIndex );

                            // get the corresponding FEM/MTK vertex's index
                            moris_index tFemVertexIndex = tVisVertex->get_base_vertex()->get_index();

                            // don't add vis vertices to the cluster, if their corresponding FEM vertices are not in the cluster either
                            if ( !tFemVertexIndexToPosInFemClusterMap.key_exists( tFemVertexIndex ) )
                            {
                                continue;
                            }

                            // check whether this VIS vertex has already been identified as being part of the VIS cluster (if so, skip, nothing needs to be done)
                            if ( !tVisVertexIndexToPosInVisClusterMap.key_exists( tVisVertexIndex ) )
                            {
                                // list the vertex as being part of this cluster
                                tVisVertexIndexToPosInVisClusterMap[ tVisVertexIndex ] = tNextVertexInVisClusterPosition;

                                // list the vertex in the cluster
                                tClusterVisVertices.push_back( tVisVertex );

                                // update next available local vertex index in the cluster
                                tNextVertexInVisClusterPosition++;
                            }

                        }    // end for: each vertex in the side cluster

                    }    // end for: each IG cell the side cluster is attached to

                    // store pointer list of VIS vertices in the VIS cluster being constructed
                    tVisSideCluster->add_vertex_to_cluster( tClusterVisVertices );

                    // ------------------------
                    // add vertex coordinates

                    // get the coordinates for all vertices in the corresponding FEM cluster
                    Matrix< DDRMat > tLocalCoordsOfAllFemVertsOnCluster =
                            tSideClustersOnFemSet( iSideClusterOnSet )->get_vertices_local_coordinates_wrt_interp_cell();

                    // initialize matrix containing all local coords for the vertices in the VIS cluster
                    Matrix< DDRMat > tVisClusterVerticesLocalCoords( tNextVertexInVisClusterPosition, tLocalCoordsOfAllFemVertsOnCluster.n_cols(), 0.0 );

                    // go through the vertices in the VIS cluster, find the corresponding FEM vertices, and copy over the local coordinates
                    for ( uint iVertInVisCluster = 0; iVertInVisCluster < tNextVertexInVisClusterPosition; iVertInVisCluster++ )
                    {
                        // get the corresponding FEM vertex's index
                        moris_index tFemVertexIndex = tClusterVisVertices( iVertInVisCluster )->get_base_vertex()->get_index();

                        // get the FEM vertex's position in the list of cluster vertices
                        moris_index tVertInFemCluster = tFemVertexIndexToPosInFemClusterMap.find( tFemVertexIndex );

                        // copy over coordinates
                        tVisClusterVerticesLocalCoords.set_row(
                                iVertInVisCluster,
                                tLocalCoordsOfAllFemVertsOnCluster.get_row( tVertInFemCluster ) );
                    }

                    // add local coordinates to vis cluster
                    tVisSideCluster->add_vertex_local_coordinates_wrt_interp_cell( tVisClusterVerticesLocalCoords );

                    // ------------------------
                    // store away the constructed cluster

                    // store away constructed cluster
                    mVisMesh->mClustersOnSideSets( iSideSet )( iSideClusterOnSet ) = tVisSideCluster;

                }    // end for: each side cluster in fem/mtk set

            }    // end for: each set in FEM mesh

        }    // end function: VIS_Factory::create_visualization_side_clusters()

        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::create_visualization_double_side_clusters()
        {
            // time and log this function
            Tracer tTracer( "VIS", "Factory", "Create Double Sided Side Clusters" );

            // loop over requested sets
            for ( uint iDblSideSet = 0; iDblSideSet < mRequestedDoubleSideSetNames.size(); iDblSideSet++ )
            {
                // get access to the fem dbl side set
                mtk::Set* tFemDblSideSet = mFemDoubleSideSets( iDblSideSet );

                // get number of side clusters on set
                uint tNumDblSideClustersOnSet = tFemDblSideSet->get_num_clusters_on_set();

                // ignore the subsequent steps for empty sets, and jump to next set
                if ( tNumDblSideClustersOnSet == 0 )
                {
                    continue;
                }

                // get list of clusters on existing fem/mtk set
                Vector< mtk::Cluster const * > tDblSideClustersOnFemSet = tFemDblSideSet->get_clusters_on_set();

                // resize list of clusters for this set
                mVisMesh->mClustersOnDoubleSideSets( iDblSideSet ).resize( tNumDblSideClustersOnSet, nullptr );
                mVisMesh->mLeaderSideClusters( iDblSideSet ).resize( tNumDblSideClustersOnSet, nullptr );
                mVisMesh->mFollowerSideClusters( iDblSideSet ).resize( tNumDblSideClustersOnSet, nullptr );

                // loop over clusters in existing fem/mtk set
                for ( uint iDblSideClusterOnSet = 0; iDblSideClusterOnSet < tDblSideClustersOnFemSet.size(); iDblSideClusterOnSet++ )
                {
                    // get access to the current fem dbl side cluster and its leader and follower components
                    mtk::Cluster const * tFemDblSideCluster      = tDblSideClustersOnFemSet( iDblSideClusterOnSet );
                    mtk::Cluster const & tFemLeaderSideCluster   = tFemDblSideCluster->get_leader_side_cluster();
                    mtk::Cluster const & tFemFollowerSideCluster = tFemDblSideCluster->get_follower_side_cluster();

                    // create the VIS side cluster objects and populate these below
                    vis::Side_Cluster_Visualization* tVisLeaderSideCluster   = new vis::Side_Cluster_Visualization();
                    vis::Side_Cluster_Visualization* tVisFollowerSideCluster = new vis::Side_Cluster_Visualization();

                    // mark as non-trivial if old clusters were not trivial
                    if ( !tFemLeaderSideCluster.is_trivial() )
                    {
                        tVisLeaderSideCluster->mark_as_nontrivial();
                    }
                    if ( !tFemFollowerSideCluster.is_trivial() )
                    {
                        tVisFollowerSideCluster->mark_as_nontrivial();
                    }

                    // ------------------------
                    // IG cells that the side cluster are attached to

                    // get the IG cells the side elements are attached to
                    const Vector< const moris::mtk::Cell* >& tFemCellsInLeaderCluster   = tFemLeaderSideCluster.get_primary_cells_in_cluster();
                    const Vector< const moris::mtk::Cell* >& tFemCellsInFollowerCluster = tFemFollowerSideCluster.get_primary_cells_in_cluster();

                    MORIS_ASSERT(
                            tFemCellsInLeaderCluster.size() == tFemCellsInFollowerCluster.size(),
                            "VIS_Factory::create_visualization_dbl_side_clusters() - "
                            "Leader and follower side clusters obtained from double side cluster have different number of facets." );
                    MORIS_ASSERT(
                            tFemCellsInLeaderCluster.size() > 0,
                            "VIS_Factory::create_visualization_dbl_side_clusters() - "
                            "Empty leader side cluster in FEM mesh. This shouldn't happen." );
                    MORIS_ASSERT(
                            tFemCellsInFollowerCluster.size() > 0,
                            "VIS_Factory::create_visualization_dbl_side_clusters() - "
                            "Empty follower side cluster in FEM mesh. This shouldn't happen." );

                    // get the number of facets in the dbl sided side cluster
                    uint tNumFacetsInSideClusters = tFemCellsInLeaderCluster.size();

                    // initialize list of VIS IG cells the facets are attached to
                    Vector< mtk::Cell const * > tLeaderSideClusterVisIgCells( tNumFacetsInSideClusters, nullptr );
                    Vector< mtk::Cell const * > tFollowerSideClusterVisIgCells( tNumFacetsInSideClusters, nullptr );

                    // find the corresponding VIS cells constructed from the original mtk/fem cells and collect these in a list for constructing the VIS side cluster
                    for ( uint iFacet = 0; iFacet < tNumFacetsInSideClusters; iFacet++ )
                    {
                        // get the FEM IG cells' indices the facet is attached to
                        moris_index tFemLeaderCellIndex   = tFemCellsInLeaderCluster( iFacet )->get_index();
                        moris_index tFemFollowerCellIndex = tFemCellsInFollowerCluster( iFacet )->get_index();

                        // get the corresponding VIS cells' indices
                        moris_index tVisLeaderCellIndex   = mPrimaryFemCellIndexToVisCellIndex( tFemLeaderCellIndex );
                        moris_index tVisFollowerCellIndex = mPrimaryFemCellIndexToVisCellIndex( tFemFollowerCellIndex );

                        // make sure the requested VIS cells actually exists
                        MORIS_ERROR(
                                tVisLeaderCellIndex > -1 && tVisLeaderCellIndex != MORIS_INDEX_MAX && tVisFollowerCellIndex > -1 && tVisFollowerCellIndex != MORIS_INDEX_MAX,
                                "VIS_Factory::create_visualization_dbl_side_clusters() - "
                                "No VIS cell index for the given FEM cell. "
                                "The current dbl side cluster may be attached to a block set that is not part of the VIS mesh. "
                                "Make sure to include it." );

                        // store the corrsponding VIS cells
                        tLeaderSideClusterVisIgCells( iFacet )   = mVisMesh->mCells( tVisLeaderCellIndex );
                        tFollowerSideClusterVisIgCells( iFacet ) = mVisMesh->mCells( tVisFollowerCellIndex );

                    }    // end for: each facet (i.e. for each IG cell)

                    // store list of cells on VIS cluster
                    tVisLeaderSideCluster->add_primary_integration_cell( tLeaderSideClusterVisIgCells );
                    tVisFollowerSideCluster->add_primary_integration_cell( tFollowerSideClusterVisIgCells );

                    // add interpolation cell to vis cluster
                    tVisLeaderSideCluster->set_interpolation_cell( &tFemLeaderSideCluster.get_interpolation_cell() );
                    tVisFollowerSideCluster->set_interpolation_cell( &tFemFollowerSideCluster.get_interpolation_cell() );

                    // ------------------------
                    // side ordinals of the IG cells that the side cluster is attached to

                    // get the side ordinals of each IG cell that make up the side cluster
                    Matrix< IndexMat > tLeaderSideOrdinals   = tFemLeaderSideCluster.get_cell_side_ordinals();
                    Matrix< IndexMat > tFollowerSideOrdinals = tFemFollowerSideCluster.get_cell_side_ordinals();

                    // store information in side cluster
                    tVisLeaderSideCluster->add_integration_cell_side_ordinals( tLeaderSideOrdinals );
                    tVisFollowerSideCluster->add_integration_cell_side_ordinals( tFollowerSideOrdinals );

                    // ------------------------
                    // add vertices to cluster
                    // NOTE: the VIS dbl side clusters will only contain vertices which sit on the interface, as the rest is not needed for visualization purposes

                    // get vertices from the mtk/fem cluster
                    Vector< mtk::Vertex const * > tFemLeaderVertices   = tFemLeaderSideCluster.get_vertices_in_cluster();
                    Vector< mtk::Vertex const * > tFemFollowerVertices = tFemFollowerSideCluster.get_vertices_in_cluster();

                    // get the number of vertices on the fem cluster for convenient access
                    uint tNumFemVerticesInLeaderCluster   = tFemLeaderVertices.size();
                    uint tNumFemVerticesInFollowerCluster = tFemFollowerVertices.size();

                    // construct map relating fem vertex indices to their respective position in the list of vertices on the clusters
                    // this is needed to later retrieve the local coordinates for the vertices
                    map< moris_index, moris_index > tFemVertexIndexToPosInFemLeaderClusterMap;
                    map< moris_index, moris_index > tFemVertexIndexToPosInFemFollowerClusterMap;
                    for ( uint iVertInFemLeaderCluster = 0; iVertInFemLeaderCluster < tNumFemVerticesInLeaderCluster; iVertInFemLeaderCluster++ )
                    {
                        moris_index tFemLeaderVertexIndex                                  = tFemLeaderVertices( iVertInFemLeaderCluster )->get_index();
                        tFemVertexIndexToPosInFemLeaderClusterMap[ tFemLeaderVertexIndex ] = iVertInFemLeaderCluster;
                    }
                    for ( uint iVertInFemFollowerCluster = 0; iVertInFemFollowerCluster < tNumFemVerticesInFollowerCluster; iVertInFemFollowerCluster++ )
                    {
                        moris_index tFemFollowerVertexIndex                                    = tFemFollowerVertices( iVertInFemFollowerCluster )->get_index();
                        tFemVertexIndexToPosInFemFollowerClusterMap[ tFemFollowerVertexIndex ] = iVertInFemFollowerCluster;
                    }

                    // initialize sets that will carry a list of FEM vertex indices that correspond to the VIS vertices
                    // which are attached to the IG cells that carry the respective sides of the dbl side cluster
                    std::set< moris_index > tActiveFemVertexIndicesOnVisLeaderCluster;
                    std::set< moris_index > tActiveFemVertexIndicesOnVisFollowerCluster;

                    // go through the IG cells the leader and follower side clusters are attached to and collect the FEM vertices on these cells
                    for ( uint iFacet = 0; iFacet < tNumFacetsInSideClusters; iFacet++ )
                    {
                        // access the FEM cells attached to the current facet
                        const moris::mtk::Cell* tLeaderFemCell   = tFemCellsInLeaderCluster( iFacet );
                        const moris::mtk::Cell* tFollowerFemCell = tFemCellsInFollowerCluster( iFacet );

                        // get the current FEM cells' vertex indices
                        Matrix< IndexMat > tLeaderCellFemVertInds   = tLeaderFemCell->get_vertex_inds();
                        Matrix< IndexMat > tFollowerCellFemVertInds = tFollowerFemCell->get_vertex_inds();

                        // go through the FEM vertices in the leader cell, and store them in the list
                        for ( uint iVertOnLeaderCell = 0; iVertOnLeaderCell < tLeaderCellFemVertInds.numel(); iVertOnLeaderCell++ )
                        {
                            moris_index tFemLeaderVertexIndex = tLeaderCellFemVertInds( iVertOnLeaderCell );
                            tActiveFemVertexIndicesOnVisLeaderCluster.insert( tFemLeaderVertexIndex );
                            // NOTE: std::set's make sure that each index will only appear once in the set;
                            // NOTE: hence, inserting the same vertex index twice does nothing
                        }

                        // go through the FEM vertices in the follower cell, and store them in the list
                        for ( uint iVertOnFollowerCell = 0; iVertOnFollowerCell < tFollowerCellFemVertInds.numel(); iVertOnFollowerCell++ )
                        {
                            moris_index tFemFollowerVertexIndex = tFollowerCellFemVertInds( iVertOnFollowerCell );
                            tActiveFemVertexIndicesOnVisFollowerCluster.insert( tFemFollowerVertexIndex );
                        }

                    }    // end for: each facet on the dbl side cluster, collect the the used vertices on each side

                    // the intersection of the active fem vertices on each side are the shared vertices,
                    // i.e. the interface vertices which we are looking for
                    std::set< moris_index > tFemVerticesOnInterface;
                    std::set_intersection(
                            tActiveFemVertexIndicesOnVisLeaderCluster.begin(),
                            tActiveFemVertexIndicesOnVisLeaderCluster.end(),
                            tActiveFemVertexIndicesOnVisFollowerCluster.begin(),
                            tActiveFemVertexIndicesOnVisFollowerCluster.end(),
                            std::inserter( tFemVerticesOnInterface, tFemVerticesOnInterface.begin() ) );

                    /* count number interface nodes and give them a position within the list of FEM vertices on the interface
                     * NOTE: this position map is used for both the leader and follower side clusters to have matching ordering
                     * of the VIS vertices on the interface, such that the vertex pairing is trivial*/
                    moris_index                     tInterfaceVertexCounter = 0;
                    map< moris_index, moris_index > tFemVertexIndexToPosInVisClusterMap;
                    for ( const moris_index& iInterfaceFemVertexIndex : tFemVerticesOnInterface )
                    {
                        tFemVertexIndexToPosInVisClusterMap[ iInterfaceFemVertexIndex ] = tInterfaceVertexCounter;
                        tInterfaceVertexCounter++;
                    }

                    uint tNumInterfaceVertices = (uint)tInterfaceVertexCounter;

                    // initialize list of VIS vertices constructed from the mtk/fem vertices on the interface
                    Vector< mtk::Vertex const * > tLeaderClusterVisVertices( tNumInterfaceVertices, nullptr );
                    Vector< mtk::Vertex const * > tFollowerClusterVisVertices( tNumInterfaceVertices, nullptr );

                    // go through the IG cells the side cluster is attached to
                    for ( uint iFacet = 0; iFacet < tNumFacetsInSideClusters; iFacet++ )
                    {
                        // get the current Fem cells' indices
                        moris_index tFemLeaderCellIndex   = tFemCellsInLeaderCluster( iFacet )->get_index();
                        moris_index tFemFollowerCellIndex = tFemCellsInFollowerCluster( iFacet )->get_index();

                        // get the block set index the side cluster is attached to
                        moris_index tLeaderPrimaryBlockSetIndex   = mPrimaryFemCellIndexToBlockIndex( tFemLeaderCellIndex );
                        moris_index tFollowerPrimaryBlockSetIndex = mPrimaryFemCellIndexToBlockIndex( tFemFollowerCellIndex );

                        // get the VIS vertex indices
                        const Vector< moris_index > tVisLeaderVertexIndices =
                                mBlockAndFemCellIndexToVisVertexIndices( tLeaderPrimaryBlockSetIndex )( tFemLeaderCellIndex );
                        const Vector< moris_index > tVisFollowerVertexIndices =
                                mBlockAndFemCellIndexToVisVertexIndices( tFollowerPrimaryBlockSetIndex )( tFemFollowerCellIndex );

                        // go through the vertices listed in the leader FEM cluster, find the corresponding VIS vertices, and store them if they are interface vertices
                        for ( uint iVertOnLeaderCell = 0; iVertOnLeaderCell < tVisLeaderVertexIndices.size(); iVertOnLeaderCell++ )
                        {
                            // get the VIS vertex's index
                            moris_index tLeaderVisVertexIndex = tVisLeaderVertexIndices( iVertOnLeaderCell );

                            // get the VIS vertex
                            mtk::Vertex const * tLeaderVisVertex = mVisMesh->mVertices( tLeaderVisVertexIndex );

                            // get the corresponding FEM/MTK vertex's index
                            moris_index tFemLeaderVertexIndex = tLeaderVisVertex->get_base_vertex()->get_index();

                            // check if this vertex is part of the interface
                            bool tVertexIsOnInterface = tFemVertexIndexToPosInVisClusterMap.key_exists( tFemLeaderVertexIndex );

                            // if vertex does not sit on the interface, skip it
                            if ( !tVertexIsOnInterface )
                            {
                                continue;
                            }

                            // get the vertex's position within cluster
                            moris_index tVertexPosInCluster = tFemVertexIndexToPosInVisClusterMap.find( tFemLeaderVertexIndex );

                            // add vertex to the list of VIS vertices in the
                            tLeaderClusterVisVertices( tVertexPosInCluster ) = tLeaderVisVertex;

                        }    // end for: each active vertex on the leader side cluster

                        // go through the vertices listed in the follower FEM cluster, find the corresponding VIS vertices, and store them if they are interface vertices
                        for ( uint iVertOnFollowerCell = 0; iVertOnFollowerCell < tVisFollowerVertexIndices.size(); iVertOnFollowerCell++ )
                        {
                            // get the VIS vertex's index
                            moris_index tFollowerVisVertexIndex = tVisFollowerVertexIndices( iVertOnFollowerCell );

                            // get the VIS vertex
                            mtk::Vertex const * tFollowerVisVertex = mVisMesh->mVertices( tFollowerVisVertexIndex );

                            // get the corresponding FEM/MTK vertex's index
                            moris_index tFemFollowerVertexIndex = tFollowerVisVertex->get_base_vertex()->get_index();

                            // check if this vertex is part of the interface
                            bool tVertexIsOnInterface = tFemVertexIndexToPosInVisClusterMap.key_exists( tFemFollowerVertexIndex );

                            // if vertex does not sit on the interface, skip it
                            if ( !tVertexIsOnInterface )
                            {
                                continue;
                            }

                            // get the vertex's position within cluster
                            moris_index tVertexPosInCluster = tFemVertexIndexToPosInVisClusterMap.find( tFemFollowerVertexIndex );

                            // add vertex to the list of VIS vertices in the
                            tFollowerClusterVisVertices( tVertexPosInCluster ) = tFollowerVisVertex;

                        }    // end for: each active vertex on the follower side cluster

                    }    // end for: each facet in the dbl side cluster

                    // store pointer lists of VIS vertices in the VIS cluster being constructed
                    tVisLeaderSideCluster->add_vertex_to_cluster( tLeaderClusterVisVertices );
                    tVisFollowerSideCluster->add_vertex_to_cluster( tFollowerClusterVisVertices );

                    // ------------------------
                    // add vertex coordinates

                    // get the coordinates for all vertices in the corresponding FEM clusters
                    Matrix< DDRMat > tLocalCoordsOfAllFemVertsOnLeaderCluster   = tFemLeaderSideCluster.get_vertices_local_coordinates_wrt_interp_cell();
                    Matrix< DDRMat > tLocalCoordsOfAllFemVertsOnFollowerCluster = tFemFollowerSideCluster.get_vertices_local_coordinates_wrt_interp_cell();

                    // deduce number of spatial dimensions from coordinates
                    uint tNumSpatialDims = tLocalCoordsOfAllFemVertsOnLeaderCluster.n_cols();

                    // initialize matrices containing all local coords for the vertices in the VIS clusters
                    Matrix< DDRMat > tLeaderVisClusterVerticesLocalCoords( tNumInterfaceVertices, tNumSpatialDims, 0.0 );
                    Matrix< DDRMat > tFollowerVisClusterVerticesLocalCoords( tNumInterfaceVertices, tNumSpatialDims, 0.0 );

                    // go through the vertices in the VIS leader cluster, find the corresponding FEM vertices, and copy over the local coordinates
                    for ( uint iVertInVisLeaderCluster = 0; iVertInVisLeaderCluster < tNumInterfaceVertices; iVertInVisLeaderCluster++ )
                    {
                        // get the corresponding FEM vertex's index
                        moris_index tLeaderFemVertexIndex = tLeaderClusterVisVertices( iVertInVisLeaderCluster )->get_base_vertex()->get_index();

                        // get the FEM vertex's position in the list of cluster vertices
                        moris_index tVertInFemLeaderCluster = tFemVertexIndexToPosInFemLeaderClusterMap.find( tLeaderFemVertexIndex );

                        // copy over coordinates
                        tLeaderVisClusterVerticesLocalCoords.set_row(
                                iVertInVisLeaderCluster,
                                tLocalCoordsOfAllFemVertsOnLeaderCluster.get_row( tVertInFemLeaderCluster ) );
                    }

                    // go through the vertices in the VIS leader cluster, find the corresponding FEM vertices, and copy over the local coordinates
                    for ( uint iVertInVisFollowerCluster = 0; iVertInVisFollowerCluster < tNumInterfaceVertices; iVertInVisFollowerCluster++ )
                    {
                        // get the corresponding FEM vertex's index
                        moris_index tFollowerFemVertexIndex = tFollowerClusterVisVertices( iVertInVisFollowerCluster )->get_base_vertex()->get_index();

                        // get the FEM vertex's position in the list of cluster vertices
                        moris_index tVertInFemFollowerCluster = tFemVertexIndexToPosInFemFollowerClusterMap.find( tFollowerFemVertexIndex );

                        // copy over coordinates
                        tFollowerVisClusterVerticesLocalCoords.set_row(
                                iVertInVisFollowerCluster,
                                tLocalCoordsOfAllFemVertsOnFollowerCluster.get_row( tVertInFemFollowerCluster ) );
                    }

                    // add local coordinates to vis clusters
                    tVisLeaderSideCluster->add_vertex_local_coordinates_wrt_interp_cell( tLeaderVisClusterVerticesLocalCoords );
                    tVisFollowerSideCluster->add_vertex_local_coordinates_wrt_interp_cell( tFollowerVisClusterVerticesLocalCoords );

                    // ------------------------
                    // store away the constructed clusters

                    /* create the VIS dbl side cluster object. the below information will be feed into it
                     * NOTE: Since we use the same ordering of the vertices on the leader and follower side the vertex pairing is trivial
                     * and the list of follower cluster vertices directly corresponds to the leader cluster vertices */
                    mtk::Double_Side_Cluster* tVisDblSideCluster =
                            new mtk::Double_Side_Cluster( tVisLeaderSideCluster, tVisFollowerSideCluster, tFollowerClusterVisVertices );

                    // store away the dbl sided side cluster
                    mVisMesh->mClustersOnDoubleSideSets( iDblSideSet )( iDblSideClusterOnSet ) = tVisDblSideCluster;

                    // store away the constructed single side clusters
                    mVisMesh->mLeaderSideClusters( iDblSideSet )( iDblSideClusterOnSet )   = tVisLeaderSideCluster;
                    mVisMesh->mFollowerSideClusters( iDblSideSet )( iDblSideClusterOnSet ) = tVisFollowerSideCluster;


                }    // end for: each dbl side cluster in set

            }    // end for: each double sided side set requested for output

        }    // end function: VIS_Factory::create_visualization_double_side_clusters()

        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::create_visualization_blocks()
        {
            // log/trace this function
            Tracer tTracer( "VIS", "Factory", "Create Block Sets" );

            // Go over the block sets in the integration mesh. Store the block meta information for later writing the mesh.
            for ( uint iBlockSet = 0; iBlockSet < mVisMesh->mListOfBlocks.size(); iBlockSet++ )
            {
                // get access to the FEM set corresponding to the current block set
                mtk::Set* tFemBlockSet = mFemBlockSets( iBlockSet );

                // create block object to write mesh set information to
                mVisMesh->mListOfBlocks( iBlockSet ) = new moris::mtk::Block_Set(
                        tFemBlockSet->get_set_name(),
                        mVisMesh->mClustersOnBlockSets( iBlockSet ),
                        tFemBlockSet->get_set_colors(),
                        tFemBlockSet->get_spatial_dim() );

                // populate the block with meta information
                mVisMesh->mListOfBlocks( iBlockSet )->set_cell_topology( tFemBlockSet->get_cell_topology() );
                mVisMesh->mListOfBlocks( iBlockSet )->set_IG_cell_shape( tFemBlockSet->get_IG_cell_shape() );
            }

            // communicate block sets
            mtk::Set_Communicator tSetCommunicator( mVisMesh->mListOfBlocks );

        }    // end function: VIS_Factory::create_visualization_blocks()

        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::create_visualization_side_sets()
        {
            // log/trace this function
            Tracer tTracer( "VIS", "Factory", "Create Side Sets" );

            // Go over the sets in the integration mesh. Store the set meta information for later writing the mesh.
            for ( uint iSideSet = 0; iSideSet < mVisMesh->mListOfSideSets.size(); iSideSet++ )
            {
                // get access to the FEM set corresponding to the current block set
                mtk::Set* tFemSideSet = mFemSideSets( iSideSet );

                // create block object to write mesh set information to
                mVisMesh->mListOfSideSets( iSideSet ) = new moris::mtk::Side_Set(
                        tFemSideSet->get_set_name(),
                        mVisMesh->mClustersOnSideSets( iSideSet ),
                        tFemSideSet->get_set_colors(),
                        tFemSideSet->get_spatial_dim() );

                // populate the set with meta information
                mVisMesh->mListOfSideSets( iSideSet )->set_IG_cell_shape( tFemSideSet->get_IG_cell_shape() );
            }

            // communicate information across all sets
            mtk::Set_Communicator tSetCommunicator( mVisMesh->mListOfSideSets );

        }    // end function: VIS_Factory::create_visualization_side_sets()

        //-----------------------------------------------------------------------------------------------------------

        void
        VIS_Factory::create_visualization_double_side_sets()
        {
            // log/trace this function
            Tracer tTracer( "VIS", "Factory", "Create Double Sided Side Sets" );

            // Go over the sets in the integration mesh. Store the set meta information for later writing the mesh.
            for ( uint iDblSideSet = 0; iDblSideSet < mVisMesh->mListOfDoubleSideSets.size(); iDblSideSet++ )
            {
                // get access to the FEM set corresponding to the current block set
                mtk::Set* tFemDblSideSet = mFemDoubleSideSets( iDblSideSet );

                // create block object to write mesh set information to
                mVisMesh->mListOfDoubleSideSets( iDblSideSet ) = new moris::mtk::Double_Side_Set(
                        tFemDblSideSet->get_set_name(),
                        mVisMesh->mClustersOnDoubleSideSets( iDblSideSet ),
                        tFemDblSideSet->get_set_colors(),
                        tFemDblSideSet->get_spatial_dim() );

                // populate the set with meta information
                mVisMesh->mListOfDoubleSideSets( iDblSideSet )->set_IG_cell_shape( tFemDblSideSet->get_IG_cell_shape() );
            }

            // communicate information across all sets
            mtk::Set_Communicator tSetCommunicator( mVisMesh->mListOfDoubleSideSets );
        }

        //-----------------------------------------------------------------------------------------------------------

    }    // namespace vis
}    // namespace moris
