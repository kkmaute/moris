/*
 * cl_MTK_Mesh_STK_New.cpp
 *
 *  Created on: Sep 18, 2018
 *      Author: doble
 */
#include "Ioss_Region.h"                // for Region, NodeBlockContainer
#include <Teuchos_RCP.hpp>              // for RCP::RCP<T>, RCP::operator*, etc

#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/Selector.hpp>     // for Selector
#include <stk_mesh/base/FEMHelpers.hpp>   // for Selector
#include "stk_io/DatabasePurpose.hpp"     // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/CoordinateSystems.hpp" // for Cartesian
#include "stk_mesh/base/CreateFaces.hpp"  // for handling faces
#include "stk_mesh/base/CreateEdges.hpp"  // for handling faces
#include "stk_mesh/base/Bucket.hpp"       // for buckets
#include "stk_mesh/base/Field.hpp"    // for coordinates
#include "stk_mesh/base/GetEntities.hpp"    // for coordinates
#include "stk_mesh/base/FieldParallel.hpp"  // for handling parallel fields

#include "fn_assert.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_MTK_Mesh_STK.hpp"
namespace moris
{
namespace mtk
{
//##############################################
// Build mesh functionality
//##############################################
    // ----------------------------------------------------------------------------
    Mesh_STK::~Mesh_STK()
    {
        // Delete member variables that contain the mesh database
        delete mMeshReader;
        delete mMtkMeshBulkData;
        delete mMtkMeshMetaData;
    }
    // ----------------------------------------------------------------------------

    Mesh_STK::Mesh_STK(
            std::string    aFileName,
            MtkSetsInfo*   aSetsInfo )
    {
        // Call the function that handles the communication between stk and moris
        this->build_mesh( aFileName, aSetsInfo );
    }

    // ----------------------------------------------------------------------------

    void
    Mesh_STK::build_mesh(
            std::string    aFileName,
            MtkSetsInfo*   aSetsInfo )
    {
        //The 'generated:' syntax in fileName makes a hex mesh to be generated in memory.
        //If an Exodus file name is used instead, then the mesh is read from file. The code
        //below remains the same in either case.

        // Declare MPI communicator
        MPI_Comm aCommunicator = MPI_COMM_WORLD;
        stk::mesh::BulkData::AutomaticAuraOption aAutoAuraOption = stk::mesh::BulkData::AutomaticAuraOption::AUTO_AURA;

        // Generate MetaData and Bulk Data instances (later to be pointed to member variables)
        stk::mesh::MetaData * meshMeta = new stk::mesh::MetaData;
        stk::mesh::BulkData * meshBulk = new stk::mesh::BulkData( *meshMeta, aCommunicator, aAutoAuraOption );

        // Set member variables as pointers to meta_data and bulk_data
        mMtkMeshMetaData = ( meshMeta );
        mMtkMeshBulkData = ( meshBulk );

        // Use STK IO to populate a STK Mesh
        mMeshReader = new stk::io::StkMeshIoBroker( aCommunicator );

        // Create mesh database using the IO broker
        mMeshReader->set_bulk_data( *meshBulk );
        mMeshReader->add_mesh_database( aFileName, stk::io::READ_MESH );
        mMeshReader->create_input_mesh();

        // Include mesh fields and populate the database
        mMeshReader->add_all_mesh_fields_as_input_fields( stk::io::MeshField::CLOSEST );

        // Create nodesets and sidesets
        MORIS_ASSERT( aSetsInfo == NULL, "Sets other than the ones provided by the input file are not currently supported." );
        mSetNames.resize( 3 ); // Number of ranks for sets

        mMeshReader->populate_bulk_data();

        // Determine number of time increments on input database in region
        Teuchos::RCP<Ioss::Region> tIo_region = mMeshReader->get_input_io_region();
        int tTimestep_count                   = tIo_region->get_property( "state_count" ).get_int();

        // Loop over all time increments and get field data
        for ( int step=1; step <= tTimestep_count; step++ )
        {
            mMeshReader->read_defined_input_fields( step );
        }

        // Create mesh edge entities
        stk::mesh::create_edges( *mMtkMeshBulkData );
        // Create mesh face entities
        stk::mesh::create_faces( *mMtkMeshBulkData, true );

        // Create communication tables in parallel.
        // NOTE1: the information to be created in the function below duplicates communication-related data
        // already created in the STK mesh database. However, one can access that data on a 1 by 1 case scenario
        // which could be inefficient in certain operations. For that reason, communication lists that live in
        // the mesh class as member variables are used instead.
        // NOTE2: this information is already provided by the user in meshes generated from data, which implies
        // that this call is only required here (mesh from string or file).

        mEntityLocaltoGlobalMap = moris::Cell<moris::Matrix< IdMat >>((uint) EntityRank::END_ENUM, moris::Matrix< IndexMat >(1, 1, 0));
        mEntitySendList         = moris::Cell<moris::Cell<moris::Matrix< IndexMat >>>((uint) EntityRank::END_ENUM, moris::Cell<moris::Matrix< IndexMat >>(par_size(),moris::Matrix< IdMat >(1,1)));
        mEntityReceiveList      = moris::Cell<moris::Cell<moris::Matrix< IndexMat >>>((uint) EntityRank::END_ENUM, moris::Cell<moris::Matrix< IndexMat >>(par_size(),moris::Matrix< IdMat >(1,1)));
        create_communication_lists_and_local_to_global_map( EntityRank::NODE );
        create_communication_lists_and_local_to_global_map( EntityRank::EDGE );
        create_communication_lists_and_local_to_global_map( EntityRank::FACE );
        create_communication_lists_and_local_to_global_map( EntityRank::ELEMENT );

        set_up_vertices_and_cell();
    }

    // ----------------------------------------------------------------------------

    Mesh_STK::Mesh_STK(
            MtkMeshData   aMeshData )
    {
    }

    // ----------------------------------------------------------------------------


//##############################################
// General mesh information access
//##############################################
    uint
    Mesh_STK::get_num_entities(enum EntityRank aEntityRank) const
    {
        // Initialize
        stk::mesh::EntityRank requestedRank = this->get_stk_entity_rank(aEntityRank);

        // Get all entities from meta data
        stk::mesh::Selector allEntities = mMtkMeshMetaData->universal_part();

        uint tNumSTKEntities = (uint)stk::mesh::count_selected_entities(allEntities, mMtkMeshBulkData->buckets(requestedRank));

        return tNumSTKEntities;
    }


//##############################################
// Access Mesh Data by index Functions
//##############################################
    Matrix<IndexMat>
    Mesh_STK::get_entity_connected_to_entity_loc_inds(moris_index     aEntityIndex,
                                            enum EntityRank aInputEntityRank,
                                            enum EntityRank aOutputEntityRank) const
    {
        MORIS_ERROR(aInputEntityRank != aOutputEntityRank," Input and output entity rank cannot be the same (this is an invalid connectivity inside STK). Use get_element_to_element_loc_inds for element to element connectivity.");

        // Get Stk entity Id from local to global map
        moris_id tId = (moris_id)this->get_glb_entity_id_from_entity_loc_index(aEntityIndex, aInputEntityRank);
        stk::mesh::EntityId tStkEntityId = (stk::mesh::EntityId)tId ;

        // Call function that gets the connected entities
        stk::mesh::EntityRank tStkInputRank = this->get_stk_entity_rank(aInputEntityRank);
        stk::mesh::EntityRank tStkOutputRank = this->get_stk_entity_rank(aOutputEntityRank);
        stk::mesh::Entity tStkEntity = mMtkMeshBulkData->get_entity(tStkInputRank, tStkEntityId);

        std::vector<stk::mesh::Entity> tEntitiesConnected = this->entities_connected_to_entity_stk(&tStkEntity, tStkInputRank, tStkOutputRank);

        // Get the number of entities
        uint tNumOutputEntities = tEntitiesConnected.size();

        // Declare xtk::Mat that will contain the local entities
        Matrix< IndexMat > tLocalIndices(1, tNumOutputEntities);

        // Fill local ids to xtk::Mat
        for (uint i = 0; i < tNumOutputEntities; ++i)
        {
            tLocalIndices(0, i) = (moris_index) mMtkMeshBulkData->local_id(tEntitiesConnected[i]);
        }
        return tLocalIndices;
    }

    Matrix< IndexMat >
    Mesh_STK::get_element_connected_to_element_loc_inds(moris_index aElementIndex) const
    {
        // First get faces connected to element
        // Get Stk entity Id from local to global map
        moris_id tId = this->get_glb_entity_id_from_entity_loc_index(aElementIndex, EntityRank::ELEMENT);
        stk::mesh::EntityId tStkEntityId = (stk::mesh::EntityId) tId ;

        // Call function that gets the connected entities
        stk::mesh::EntityRank tStkInputRank  = stk::topology::ELEMENT_RANK;
        stk::mesh::EntityRank tStkOutputRank = stk::topology::FACE_RANK;
        stk::mesh::Entity tStkEntity = mMtkMeshBulkData->get_entity(tStkInputRank, tStkEntityId);

        std::vector<stk::mesh::Entity> tFacesInElem = this->entities_connected_to_entity_stk(&tStkEntity, tStkInputRank, tStkOutputRank);


        MORIS_ASSERT( ( tFacesInElem.size() != 0 ) || ( tFacesInElem.size() != 0 ),
                "No faces connected to element found. Maybe the CreateAllEdgesAndFaces flag is set to false. Check mesh struct." );

        // Then for each face get elements connected
        uint tCounter  = 0;
        uint tNumFaces = tFacesInElem.size();

        moris::Matrix< IndexMat > tElemsConnectedToElem(2, tNumFaces);

        for ( uint faceIt = 0; faceIt < tNumFaces; ++faceIt )
        {
            std::vector<stk::mesh::Entity> tDummyConnectivity = this->entities_connected_to_entity_stk( &tFacesInElem[faceIt],stk::topology::FACE_RANK, stk::topology::ELEMENT_RANK );

            // Faces in mesh boundaries do not have more than one element
            if ( tDummyConnectivity.size() > 0 )
            {
                if ( mMtkMeshBulkData->identifier(tDummyConnectivity[0]) !=
                     mMtkMeshBulkData->identifier(  tStkEntity ) )
                {
                    tElemsConnectedToElem( 0,tCounter ) = (moris_index) mMtkMeshBulkData->local_id(tDummyConnectivity[0]);
                    tElemsConnectedToElem( 1,tCounter ) = (moris_index) faceIt;
                    tCounter++;
                }
            }

            if ( tDummyConnectivity.size()  > 1 )
            {
                if ( mMtkMeshBulkData->identifier(tDummyConnectivity[1]) !=
                     mMtkMeshBulkData->identifier(  tStkEntity ) )
                {
                    tElemsConnectedToElem( 0,tCounter ) = (moris_index) mMtkMeshBulkData->local_id(tDummyConnectivity[1]);
                    tElemsConnectedToElem( 1,tCounter ) = (moris_index) faceIt;
                    tCounter++;
                }
            }

            MORIS_ASSERT( tDummyConnectivity.size()  <= 2,
                    "For some reason face has more than 2 elements connected to it... Check get_elements_connected_to_element." );
        }

        // Resize to include only ids added above and get rid of initialized extra zeros
        tElemsConnectedToElem.resize( 2,tCounter );

        return tElemsConnectedToElem;
    }

    // ----------------------------------------------------------------------------

//##############################################
// Local to global functions
//##############################################

    moris_id
    Mesh_STK::get_glb_entity_id_from_entity_loc_index(moris_index     aEntityIndex,
                                                      enum EntityRank aEntityRank) const
    {
       return mEntityLocaltoGlobalMap((uint)aEntityRank)(0,aEntityIndex);
    }

    // ----------------------------------------------------------------------------

    Matrix<IdMat>
    Mesh_STK::get_entity_connected_to_entity_glob_ids(
            moris_id     aEntityId,
            enum EntityRank aInputEntityRank,
            enum EntityRank aOutputEntityRank)
    {

        // Call function that gets the connected entities
        stk::mesh::EntityRank tStkInputRank = this->get_stk_entity_rank(aInputEntityRank);
        stk::mesh::EntityRank tStkOutputRank = this->get_stk_entity_rank(aOutputEntityRank);
        stk::mesh::Entity tStkEntity = mMtkMeshBulkData->get_entity(tStkInputRank, (stk::mesh::EntityId)aEntityId);

        std::vector<stk::mesh::Entity> tSTKEntitiesConnectedGlobIds = this->entities_connected_to_entity_stk(&tStkEntity, tStkInputRank, tStkOutputRank);

        // Get the number of entities
        uint tNumOutputEntities = tSTKEntitiesConnectedGlobIds.size();

        // Declare xtk::Mat that will contain the local entities
        Matrix< IndexMat > tEntitiesConnectedGlobIds(1, tNumOutputEntities);

        // Fill local ids to xtk::Mat
        for (uint i = 0; i < tNumOutputEntities; ++i)
        {
            tEntitiesConnectedGlobIds(0, i) = (moris_index) mMtkMeshBulkData->identifier(tSTKEntitiesConnectedGlobIds[i]);
        }

        return tEntitiesConnectedGlobIds;
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Mesh_STK::get_element_connected_to_element_glob_ids(moris_id aElementId) const
    {
        // First get faces connected to element
        // Get Stk entity Id from local to global map
        stk::mesh::EntityId tStkEntityId = (stk::mesh::EntityId) aElementId ;

        // Call function that gets the connected entities
        stk::mesh::EntityRank tStkInputRank  = stk::topology::ELEMENT_RANK;
        stk::mesh::EntityRank tStkOutputRank = stk::topology::FACE_RANK;
        stk::mesh::Entity tStkEntity = mMtkMeshBulkData->get_entity(tStkInputRank, tStkEntityId);

        std::vector<stk::mesh::Entity> tFacesInElem = this->entities_connected_to_entity_stk(&tStkEntity, tStkInputRank, tStkOutputRank);


        MORIS_ASSERT( ( tFacesInElem.size() != 0 ) || ( tFacesInElem.size() != 0 ),
                      "No faces connected to element found. Maybe the CreateAllEdgesAndFaces flag is set to false. Check mesh struct." );

        // Then for each face get elements connected
        uint tCounter  = 0;
        uint tNumFaces = tFacesInElem.size();

        moris::Matrix< IdMat > tElemsConnectedToElem(2, tNumFaces);

        for ( uint faceIt = 0; faceIt < tNumFaces; ++faceIt )
        {
            std::vector<stk::mesh::Entity> tDummyConnectivity = this->entities_connected_to_entity_stk( &tFacesInElem[faceIt],stk::topology::FACE_RANK, stk::topology::ELEMENT_RANK );

            // Faces in mesh boundaries do not have more than one element
            if ( tDummyConnectivity.size() > 0 )
            {
                if ( mMtkMeshBulkData->identifier(tDummyConnectivity[0]) !=
                        mMtkMeshBulkData->identifier(  tStkEntity ) )
                {
                    tElemsConnectedToElem( 0,tCounter ) = (moris_id) mMtkMeshBulkData->identifier(tDummyConnectivity[0]);
                    tElemsConnectedToElem( 1,tCounter ) = (moris_id) faceIt;
                    tCounter++;
                }
            }

            if ( tDummyConnectivity.size()  > 1 )
            {
                if ( mMtkMeshBulkData->identifier(tDummyConnectivity[1]) !=
                        mMtkMeshBulkData->identifier(  tStkEntity ) )
                {
                    tElemsConnectedToElem( 0,tCounter ) = (moris_id) mMtkMeshBulkData->identifier(tDummyConnectivity[1]);
                    tElemsConnectedToElem( 1,tCounter ) = (moris_id) faceIt;
                    tCounter++;
                }
            }

            MORIS_ASSERT( tDummyConnectivity.size()  <= 2,
                          "For some reason face has more than 2 elements connected to it... Check get_elements_connected_to_element." );
        }

        // Resize to include only ids added above and get rid of initialized extra zeros
        tElemsConnectedToElem.resize( 2,tCounter );

        return tElemsConnectedToElem;
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Mesh_STK::generate_unique_entity_ids( uint            aNumEntities,
                                          enum EntityRank aEntityRank) const
    {
        std::vector<stk::mesh::EntityId> tAvailableNodeIDs; // generate_new_ids requires a variable of this type
        mMtkMeshBulkData->generate_new_ids( get_stk_entity_rank(aEntityRank), aNumEntities, tAvailableNodeIDs );

        Matrix < IdMat > aAvailableNodeIDs ( aNumEntities, 1 );

        for ( uint i = 0; i < aNumEntities; i++ )
        {
            aAvailableNodeIDs( i, 0 ) = (moris_id) tAvailableNodeIDs[i];
        }

        return aAvailableNodeIDs;
    }

//##############################################
// Coordinate Field Functions
//##############################################

    Matrix< DDRMat >
    Mesh_STK::get_node_coordinate( moris_index aNodeIndex) const
    {
        // Get the coordinate field from stk
        stk::mesh::FieldBase const * coord = mMtkMeshMetaData->coordinate_field();

        // Get node id from provided index
        moris_id tId = get_glb_entity_id_from_entity_loc_index(aNodeIndex,EntityRank::NODE);

        stk::mesh::EntityId tNodeId = (stk::mesh::EntityId) tId ;

        // Declare node entity
        stk::mesh::Entity tNodeEntity = mMtkMeshBulkData->get_entity(stk::topology::NODE_RANK, tNodeId);

        // Get coordinates of node n
        double *fieldValue = static_cast<double *>(stk::mesh::field_data(*coord, tNodeEntity));

        Matrix < DDRMat > tNodeCoord (1,get_spatial_dim());


        for (uint dim = 0; dim < get_spatial_dim(); dim++)
        {
            tNodeCoord(0, dim) = fieldValue[dim];
        }

        return tNodeCoord;
    }

    //##############################################
    // Entity Ownership Functions
    //##############################################
    moris_id
    Mesh_STK::get_entity_owner( moris_index     aEntityIndex,
                                enum EntityRank aEntityRank ) const
    {

        // Convert index to ID
        moris_id tEntityId = get_glb_entity_id_from_entity_loc_index(aEntityIndex, aEntityRank);

        //Get entity Id
        stk::mesh::Entity tEntity = mMtkMeshBulkData->get_entity(get_stk_entity_rank(aEntityRank), tEntityId);

        //processor rank that owns entity
        moris_id tOwningProcessor = mMtkMeshBulkData->parallel_owner_rank(tEntity);

        return tOwningProcessor;
    }

    // ----------------------------------------------------------------------------

//##############################################
// moris::Cell and Vertex Pointer Functions
//##############################################

    mtk::Cell const &
    Mesh_STK::get_mtk_cell(moris_index aCellIndex)
    {
        return mMtkCells(aCellIndex);
    }

    // ----------------------------------------------------------------------------

    mtk::Vertex const &
    Mesh_STK::get_mtk_vertex(moris_index aVertexIndex)
    {
        return mMtkVertices(aVertexIndex);
    }

//##############################################
// Private functions to build mesh
//##############################################
    void
    Mesh_STK::create_communication_lists_and_local_to_global_map(enum EntityRank aEntityRank)
    {
        const int tParallelSize = mMtkMeshBulkData->parallel_size();
        const int tParallelRank = mMtkMeshBulkData->parallel_rank();

        // Declare vector of entity counts
        std::vector<uint> tEntityCounts;

        // Get all entities from meta data
        stk::mesh::Selector tSharedSelector = mMtkMeshMetaData->universal_part();

        // Count entities
        stk::mesh::count_entities( tSharedSelector, *mMtkMeshBulkData, tEntityCounts );

        uint tNumEntities = static_cast<uint>(tEntityCounts[ (uint)aEntityRank ]);

        // Resize comm lists to maximum possible
        for(int i = 0; i< tParallelSize; i++)
        {
            moris::Matrix< IndexMat > tSendMat(1,tNumEntities,(uint)0);
            moris::Matrix< IndexMat > tRecvMat(1,tNumEntities,(uint)0);

            mEntitySendList((uint)aEntityRank)(i) = tSendMat;
            mEntityReceiveList((uint)aEntityRank)(i) = tRecvMat;
        }

        moris::Matrix< IndexMat > tMapMat(1,tNumEntities,(uint)0);
        mEntityLocaltoGlobalMap((uint)aEntityRank)= tMapMat;

        stk::mesh::BucketVector const& shared_node_buckets =
                mMtkMeshBulkData->get_buckets( get_stk_entity_rank(aEntityRank) , tSharedSelector);

        uint tCurrentIndex = 0;
        // Initialize proc counter
        // moris::Cell #  = Proc rank
        moris::Cell<uint> tSendProcCounter(tParallelSize);
        moris::Cell<uint> tRecvProcCounter(tParallelSize);

        // Loop over shared nodes
        for(uint i = 0; i<shared_node_buckets.size();i++)
        {
            stk::mesh::Bucket& bucket = *shared_node_buckets[i];

            for(uint j = 0; j<bucket.size(); j++)
            {
                moris_id tEntityId = (moris_id) mMtkMeshBulkData->identifier(bucket[j]);
                int tOwnerProcRank = mMtkMeshBulkData->parallel_owner_rank(bucket[j]);

                // Set local to global map in mesh and STK
                mEntityLocaltoGlobalMap((uint)aEntityRank)(0, tCurrentIndex) = tEntityId;

                mMtkMeshBulkData->set_local_id( bucket[j], tCurrentIndex );

                std::vector<int> sharedProcs;

                // Get shared procs Ids
                mMtkMeshBulkData->comm_procs(mMtkMeshBulkData->entity_key(bucket[j]),sharedProcs);

                if(sharedProcs.size() != 0)
                {

                    // Sort (if current proc owns entity then add to send comm lists)
                    if(tOwnerProcRank == tParallelRank)
                    {
                        // loop over processors
                        for(uint p = 0; p<sharedProcs.size(); p++)
                        {
                            if(sharedProcs[p]!=mMtkMeshBulkData->parallel_rank())
                            {
                                uint tSharedProcRank = sharedProcs[p];
                                uint tSendCount = tSendProcCounter(tSharedProcRank);
                                mEntitySendList((uint)aEntityRank)(tSharedProcRank)(0,tSendCount) = tCurrentIndex;
                                tSendProcCounter(tSharedProcRank)++;
                            }
                        }
                    }
                    // (if current proc does not own entity then add to recv comm lists)
                    else if (tOwnerProcRank != tParallelRank)
                    {
                        for(uint p = 0; p<sharedProcs.size(); p++)
                        {
                            if(sharedProcs[p]!=tParallelRank)
                            {
                                uint tSharedProc = sharedProcs[p];
                                uint tRecvCount = tRecvProcCounter(tSharedProc);
                                mEntityReceiveList((uint)aEntityRank)(tSharedProc)(0,tRecvCount) = tCurrentIndex;
                                tRecvProcCounter(tSharedProc)++;
                            }
                        }
                    }
                }

                tCurrentIndex++;
            }
        }

        for(int pr = 0; pr<tParallelSize;pr++)
        {
            uint tRecvCount = tRecvProcCounter(pr);
            uint tSendCount = tSendProcCounter(pr);
            mEntitySendList((uint)aEntityRank)(pr).resize(1,tSendCount);
            mEntityReceiveList((uint)aEntityRank)(pr).resize(1,tRecvCount);
        }
    }

    void
    Mesh_STK::set_up_vertices_and_cell()
    {
        // Get information about the mesh
        uint tNumNodes        = this->get_num_entities(EntityRank::NODE);

        // Setup vertices
        mMtkVertices = moris::Cell<Vertex_STK>(tNumNodes);

        for( uint iVertInd = 0; iVertInd<tNumNodes; iVertInd++)
        {
            // pass global node ids, node index and a pointer to this mesh into the vertex
            mMtkVertices(iVertInd) = Vertex_STK(this->get_glb_entity_id_from_entity_loc_index(iVertInd,EntityRank::NODE),
                                                iVertInd,
                                                this);
        }

        // Setup Cells
        uint tNumElems        = this->get_num_entities(EntityRank::ELEMENT);
        // allocate member data
        mMtkCells = moris::Cell<mtk::Cell_STK>(tNumElems);
        Matrix< IndexMat > tElementToNode;

        for( moris_index iCellInd = 0; iCellInd<(moris_index)tNumElems; iCellInd++)
        {
            tElementToNode = get_entity_connected_to_entity_loc_inds(iCellInd, EntityRank::ELEMENT, EntityRank::NODE);

            // setup vertices of cells
            moris::Cell<Vertex*> tElementVertices(tElementToNode.numel());
            for(uint iNodes = 0; iNodes<tElementToNode.numel(); iNodes++)
            {
                tElementVertices(iNodes) = &mMtkVertices(tElementToNode(iNodes));
            }

            // Add cell to member data
            mMtkCells(iCellInd) = Cell_STK( CellType::HEX8,
                                            this->get_glb_entity_id_from_entity_loc_index(iCellInd,EntityRank::ELEMENT),
                                            iCellInd,
                                            tElementVertices,
                                            this);
        }
    }
// ----------------------------------------------------------------------------

//##############################################
// Private functions to access mesh information
//##############################################
    stk::mesh::EntityRank
    Mesh_STK::get_stk_entity_rank(enum EntityRank aMRSEntityRank) const
    {
        if (aMRSEntityRank == EntityRank::NODE)
        {
            return stk::topology::NODE_RANK;
        }
        else if (aMRSEntityRank == EntityRank::EDGE)
        {
            return stk::topology::EDGE_RANK;
        }

        else if (aMRSEntityRank == EntityRank::FACE)
        {
            return stk::topology::FACE_RANK;
        }
        else if (aMRSEntityRank == EntityRank::ELEMENT)
        {
            return stk::topology::ELEMENT_RANK;
        }
        else
        {
            return stk::topology::INVALID_RANK;
        }
    }

    // ----------------------------------------------------------------------------

//##############################################
// internal id functions
//##############################################
    std::vector<stk::mesh::Entity>
    Mesh_STK::entities_connected_to_entity_stk(
            stk::mesh::Entity*    const aInputEntity,
            stk::mesh::EntityRank const aInputEntityRank,
            stk::mesh::EntityRank const aOutputEntityRank) const
    {
        // Declare the object where we are going to store the shared faces and handlers
        std::vector<stk::mesh::Entity> tDesiredEntitiesConnectedToInputEntities;

        // Check if the connectivity exists (i.e., was already generated and is stored in mesh data)
        if (mMtkMeshBulkData->connectivity_map().valid(aInputEntityRank, aOutputEntityRank))
        {

            switch (aOutputEntityRank)
            {
                case stk::topology::NODE_RANK:

                    // Fill entities connected
                    if (mMtkMeshBulkData->num_nodes(aInputEntity[0]) > 0)
                    {
                        // Get pointers to the location of the connected nodes
                        stk::mesh::Entity const * tDesiredEntityStart = mMtkMeshBulkData->begin_nodes(aInputEntity[0]);
                        stk::mesh::Entity const * tDesiredEntityEnd = mMtkMeshBulkData->end_nodes(aInputEntity[0]);

                        // Store faces in output vector
                        tDesiredEntitiesConnectedToInputEntities.assign(tDesiredEntityStart, tDesiredEntityEnd);
                    }
                    break;

                case stk::topology::EDGE_RANK:

                    // Fill entities connected
                    if (mMtkMeshBulkData->num_edges(aInputEntity[0]) > 0)
                    {
                        // Get pointers to the location of the connected edges
                        stk::mesh::Entity const * tDesiredEntityStart = mMtkMeshBulkData->begin_edges(aInputEntity[0]);
                        stk::mesh::Entity const * tDesiredEntityEnd = mMtkMeshBulkData->end_edges(aInputEntity[0]);

                        // Store faces in output vector
                        tDesiredEntitiesConnectedToInputEntities.assign(tDesiredEntityStart, tDesiredEntityEnd);
                    }
                    break;

                case stk::topology::FACE_RANK:

                    // Fill entities connected
                    if (mMtkMeshBulkData->num_faces(aInputEntity[0]) > 0)
                    {
                        // Get pointers to the location of the connected faces
                        stk::mesh::Entity const * tDesiredEntityStart = mMtkMeshBulkData->begin_faces(aInputEntity[0]);
                        stk::mesh::Entity const * tDesiredEntityEnd = mMtkMeshBulkData->end_faces(aInputEntity[0]);

                        // Store faces in output vector
                        tDesiredEntitiesConnectedToInputEntities.assign(tDesiredEntityStart, tDesiredEntityEnd);
                    }
                    break;

                case stk::topology::ELEMENT_RANK:

                    // Fill entities connected
                    if (mMtkMeshBulkData->num_elements(aInputEntity[0]) > 0)
                    {
                        // Get pointers to the location of the connected elements
                        stk::mesh::Entity const * tDesiredEntityStart = mMtkMeshBulkData->begin_elements(aInputEntity[0]);
                        stk::mesh::Entity const * tDesiredEntityEnd = mMtkMeshBulkData->end_elements(aInputEntity[0]);

                        // Store faces in output vector
                        tDesiredEntitiesConnectedToInputEntities.assign(tDesiredEntityStart, tDesiredEntityEnd);
                    }
                    break;

                default:
                    std::cerr << " wrong topology in entities_connected_to_entity_stk ";
                    break;
            }
        }
        else
        {
            std::cerr << " STK already has valid connectivity maps. Check if you are trying to access invalid connectivity (e.g., edge to edge)";
        }
        return tDesiredEntitiesConnectedToInputEntities;
    }

    }
}


