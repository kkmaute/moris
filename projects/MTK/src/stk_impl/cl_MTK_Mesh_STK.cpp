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
#include "fn_isempty.hpp"
#include "fn_find.hpp"
#include "op_equal_equal.hpp"

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
            MtkMeshData   aMeshData ):
        mEntityLocaltoGlobalMap(4),
        mSetRankFlags( { false, false, false} )
    {
        // Flag for handling data generated mesh
        mDataGeneratedMesh = true;

        // Verify the minimum required data was provided and initialize uninitialized variables if necessary
        this->check_and_update_input_data( aMeshData );

        // Call the function that handles the communication between stk and moris
        this->build_mesh( aMeshData );
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

    moris_index
    Mesh_STK::get_loc_entity_ind_from_entity_glb_id(moris_id        aEntityId,
                                                    enum EntityRank aEntityRank) const
    {
        return mMtkMeshBulkData->local_id( mMtkMeshBulkData->get_entity( get_stk_entity_rank(aEntityRank), aEntityId ) );
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

    // ----------------------------------------------------------------------------
    void
    Mesh_STK::create_output_mesh(
            std::string  &aFileName )
    {
        if ( mDataGeneratedMesh )
            {
                // Generate data for mesh from mesh reader
                size_t outputFileIdx = mMeshReader->create_output_mesh( aFileName,stk::io::WRITE_RESULTS );
                mMeshReader->write_output_mesh( outputFileIdx );

                // Get fields and initialize fields information
                const std::vector<stk::mesh::FieldBase*>&fields = mMtkMeshMetaData->get_fields();

                // Add fields to output mesh
                std::string tFieldNoData = "dummyField";
                std::string tCoordField  = "coordinates";
                std::vector<stk::mesh::FieldBase *>::const_iterator fieldIterator = fields.begin();
                for ( ; fieldIterator != fields.end();++fieldIterator )
                {
                    // Get field name
                    std::string tIterFieldName = ( *fieldIterator )->name();

                    // Do not add dummy or coordinate fields to the output mesh
                    if ( ( tIterFieldName.compare( tFieldNoData ) != 0 ) && ( tIterFieldName.compare( tCoordField ) != 0 ) )
                    {
                        mMeshReader->add_field( outputFileIdx, *( stk::mesh::get_field_by_name( ( *fieldIterator )->name(), *mMtkMeshMetaData ) ) );
                    }
                }

                // Provisionally only handles static problems (hard-coded time)
                mMeshReader->begin_output_step( outputFileIdx, 0 );
                mMeshReader->write_defined_output_fields( outputFileIdx );
                mMeshReader->end_output_step( outputFileIdx );
            }
            else
            {
                // Generate data for mesh from mesh reader
                size_t fh = mMeshReader->create_output_mesh( aFileName, stk::io::WRITE_RESULTS );

                // write mesh with the information generated from the mesh reader
                mMeshReader->write_output_mesh( fh );
            }
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
            mMtkCells(iCellInd) = Cell_STK( CellTopology::HEX8,
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
    Matrix< IdMat >
    Mesh_STK::get_entities_owned_and_shared_by_current_proc(
            EntityRank   aEntityRank ) const
    {
        if ( aEntityRank == EntityRank::NODE )
        {
            return mEntityLocaltoGlobalMap(0);
        }
        else if ( aEntityRank == EntityRank::EDGE )
        {
            return mEntityLocaltoGlobalMap(1);
        }
        else if ( aEntityRank == EntityRank::FACE )
        {
            return mEntityLocaltoGlobalMap(2);
        }
        else if ( aEntityRank == EntityRank::ELEMENT )
        {
            return mEntityLocaltoGlobalMap(3);
        }
        else
        {
            MORIS_ASSERT( 0, "Invalid rank provided in get_entities_owned_and_shared_by_current_proc." );
        }

        Matrix< IdMat >  tDummyConn( 1, 1 );
        return tDummyConn;
    }

    // ----------------------------------------------------------------------------
    // Access entities in selected portion of the mesh
    moris::Matrix< DDUMat >
    Mesh_STK::get_entities_in_selector_interface(
            EntityRank            aRequestedEntityRank,
            stk::mesh::Selector   aSelectedEntities ) const
    {
        // Get selected entities
        std::vector<stk::mesh::Entity>tOutputEntityIDs;
        stk::mesh::get_selected_entities(
                aSelectedEntities, mMtkMeshBulkData->buckets( this->get_stk_entity_rank( aRequestedEntityRank ) ), tOutputEntityIDs );

        // Interface with STK, get the Ids and return them
        uint tNumEntity = tOutputEntityIDs.size();
        Matrix< DDUMat >  tOutputEntityIDMat( tNumEntity, 1, 0 );

        for ( uint i = 0; i<tNumEntity; i++)
        {
            tOutputEntityIDMat( i, 0 ) = ( uint ) mMtkMeshBulkData->identifier( tOutputEntityIDs[i] );
        }

        return tOutputEntityIDMat;
    }

    // ----------------------------------------------------------------------------

    // return processors sharing a particular entity
    moris::Matrix< DDUMat >
    Mesh_STK::get_procs_sharing_entity_by_id(
            uint              aEntityID,
            enum EntityRank   aEntityRank ) const
    {
        // Initialize returning mat with UINT_MAX in case no processors are sharing the given entity
        Matrix< DDUMat >  tSharedProcsMat( 1, 1, UINT_MAX );

        // Check if function is called in a serial run. If so, no parallel information is available.
        uint tProcSize = par_size();
        if ( tProcSize == 1 )
        {
            return tSharedProcsMat;
        }

        // Define entity
        const stk::mesh::Entity tEntity = mMtkMeshBulkData->get_entity( this->get_stk_entity_rank( aEntityRank ), aEntityID );

        // Intialize shared procs
        std::vector<int> tSharedProcsVec;

        // get shared processor IDs
        mMtkMeshBulkData->comm_shared_procs( mMtkMeshBulkData->entity_key( tEntity ), tSharedProcsVec );

        uint tNumSharedProcs = tSharedProcsVec.size();

        // Check if no processors are sharing the given entity
        if ( tNumSharedProcs != 0 )
        {
            tSharedProcsMat.resize( tNumSharedProcs, 1 );
        }

        // Transform from standard to moris call
        for ( uint iProc = 0; iProc < tNumSharedProcs; ++iProc )
        {
            tSharedProcsMat( iProc ) = tSharedProcsVec[iProc];
        }

        return tSharedProcsMat;
    }

    // ----------------------------------------------------------------------------

    moris::Cell < moris::Cell < uint > >
    Mesh_STK::get_shared_info_by_entity( uint aNumActiveSharedProcs, enum EntityRank  aEntityRank )
    {
        // Generate elements shared per processor list
            // -------------------------------------------
            Matrix< DDUMat >  tEntitiesShared = this->get_entities_glb_shared_current_proc( aEntityRank );
            uint tNumEntitiesShared = tEntitiesShared.length();
            moris_id tParallelRank  = par_rank();

            moris::Cell < moris::Cell < uint > > tTemporaryEntityMapSharingProcs( aNumActiveSharedProcs );

            // Loop over the number of nodes shared to get the shared processors
            for( uint iElemShared = 0; iElemShared < tNumEntitiesShared; ++iElemShared )
            {
                Matrix< DDUMat >  tProcsSharing = this->get_procs_sharing_entity_by_id( tEntitiesShared( iElemShared ), aEntityRank );
                for( uint iProc = 0; iProc < tProcsSharing.length(); ++iProc )
                {
                    if( tProcsSharing(iProc) != (uint)tParallelRank )
                    {
                        tTemporaryEntityMapSharingProcs( mProcsSharedToIndex[ tProcsSharing( iProc ) ] ).push_back( tEntitiesShared( iElemShared ) );
                    }
                }
            }

            return tTemporaryEntityMapSharingProcs;
    }
    // ----------------------------------------------------------------------------
    //##############################################
    // Build mesh from data functions internal
    //##############################################
    // Verifications of mesh essential information provided for meshes generated from data.
    void
    Mesh_STK::check_and_update_input_data(
            MtkMeshData&   aMeshData )
    {
        // Mesh main variables
        uint tNumNodes = aMeshData.NodeCoords[0].n_rows();
        uint tNumElems = aMeshData.ElemConn[0].n_rows();

        // Parallel initialization
        uint tProcSize = par_size();
        uint tProcRank = par_rank();
        // Access mesh data from struc and check input values for indispensable arguments
        MORIS_ASSERT( aMeshData.SpatialDim != NULL, "Number of spatial dimensions was not provided." );
        MORIS_ASSERT( aMeshData.ElemConn   != NULL, "Element connectivity was not provided.");
        MORIS_ASSERT( aMeshData.NodeCoords != NULL, "Node coordinates were not provided." );

        // Initialize number of dimensions
        mNumDims = aMeshData.SpatialDim[0];

        if ( ( aMeshData.EntProcOwner == NULL ) && ( tProcSize == 1 ) )
        {
            // Do nothing. This information is not needed in serial
        }
        else if ( aMeshData.EntProcOwner != NULL )
        {
            // Verify sizes
            MORIS_ASSERT( ( aMeshData.EntProcOwner[0].n_rows() == tNumNodes ) ||
                          ( aMeshData.EntProcOwner[0].n_rows() == tNumElems ),
                          "Number of rows for EntProcOwner should match number of nodes or elements." );
        }
        else
        {
            MORIS_ASSERT( 0, "Elements or nodes processor owner list must be provided in parallel." );
        }

        // If nodal map was not provided and the simulation is run with 1 processor,
        // just generate the map. The number of ids provided in the map should match
        // the number of coordinates given for each processor.
        if ( aMeshData.LocaltoGlobalNodeMap != NULL )
        {
            // Verify sizes
            MORIS_ASSERT( aMeshData.LocaltoGlobalNodeMap[0].n_rows() == tNumNodes, "Number of rows for LocaltoGlobalNodeMap should match number nodes." );

            mEntityLocaltoGlobalMap(0) = aMeshData.LocaltoGlobalNodeMap[0];
        }
        else if ( ( aMeshData.LocaltoGlobalNodeMap == NULL ) && ( tProcSize == 1 ) )
        {
            // Generate nodes maps
            mEntityLocaltoGlobalMap(0) = Matrix<IdMat>( tNumNodes, 1 );
            for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
            {
                mEntityLocaltoGlobalMap(0)( iNode ) = iNode + 1;
            }

            aMeshData.LocaltoGlobalNodeMap = &mEntityLocaltoGlobalMap(0);
        }
        else
        {
            // This set of input parameters only works in serial
            MORIS_ASSERT( 0, "You have to provide a local to global node map for coordinates if you want to create a mesh in parallel. " );
        }

        // If element map was not provided and the simulation, just generate the map.
        // WARNING: A threshold value that assumes a maximum number of elements
        //          per processors is used. This number needs to be changed if any
        //          processor has more elements than the threshold.
        if ( aMeshData.LocaltoGlobalElemMap != NULL )
        {
            // Verify sizes
            MORIS_ASSERT( aMeshData.LocaltoGlobalElemMap[0].n_rows() == tNumElems,
                         "Number of rows for LocaltoGlobalElemMap should match number of elements provided in connectivity map" );

            mEntityLocaltoGlobalMap(3) = aMeshData.LocaltoGlobalElemMap[0];
        }
        else
        {
            // Threshold value and number of elements in current processor
            uint tThresholdValue = 1000000;

            // Check that no one processor exceeds the threshold value.
            MORIS_ASSERT( tNumElems <= tThresholdValue, "Constructor generates a maximum of 1000000 elements per processor in current constructor. "
                                         "Change Threshold value or provide element local to global map." );

            // Generate IDs assuming there are no more than tThresholdValue elements per processor
            mEntityLocaltoGlobalMap(3) = Matrix<IdMat>( tNumElems, 1 );
            for ( uint iElem = 0; iElem < tNumElems; ++iElem )
            {
                mEntityLocaltoGlobalMap(3)( iElem ) = tThresholdValue*tProcRank + iElem + 1;
            }
            aMeshData.LocaltoGlobalElemMap = &mEntityLocaltoGlobalMap(3);
        }

        // No problem if field information was not provided (i.e., FieldsData, FieldsName, PartNames)
        // Only need to check if data given is consistent.
        if ( aMeshData.FieldsInfo != NULL )
        {
            this->check_and_update_fields_data( aMeshData );
        }
        mSetNames.resize( 3 ); // Number of ranks for set names
        // No problem if (block,node,side) sets information was not provided
        // Only need to check if data given is consistent.

        if ( aMeshData.SetsInfo != NULL )
        {
            this->check_and_update_sets_data( aMeshData );
        }
        else
        {
            // Create a block set tha contains the entire mesh by default
            mSetNames[2].resize( 1, "block_1" );
        }
    }

    // ----------------------------------------------------------------------------

    // Verification of fields data arrangement and consistency for meshes generated from data.
    void
    Mesh_STK::check_and_update_fields_data( MtkMeshData&   aMeshData )
    {
        mFieldInDataGiven = true;

        // Get the number of fields
        uint tNumFields = aMeshData.FieldsInfo[0].FieldsData[0].size();

        // Verify that all field ranks were given
        MORIS_ASSERT( aMeshData.FieldsInfo[0].FieldsRank.size() == tNumFields, "Number of field ranks should be the same as number of field data Mats." );

        // Check if set owner names were provided
        if ( aMeshData.FieldsInfo[0].SetsOwner != NULL )
        {
            MORIS_ASSERT( aMeshData.FieldsInfo[0].SetsOwner[0].size() == tNumFields ,
                    "Set owner container should have names for all fields declared. "
                    "If field is declared over universal part, provide empty string.");
        }

        // Loop over the number of fields
        for ( uint iField = 0; iField < tNumFields; ++iField )
        {
            // Verify that field sizes (number of columns) match the ones suppported
            Matrix< DDUMat >  tSupFieldComp = { {1, 2, 3, 4, 9} };
            uint tNumFieldComp = aMeshData.FieldsInfo[0].FieldsData[0]( iField ).n_cols();
            Matrix< DDBMat >  tDummy = ( tSupFieldComp == tNumFieldComp );
            Matrix< DDNIMat >  tCompFound = find ( tDummy );


            MORIS_ASSERT( !isempty( tCompFound ),
                    "Number of components (columns) for all FieldsData should "
                    "match one of the supported sizes {1, 2, 3, 4, 9}.");

            // Check if field names were provided
            if ( aMeshData.FieldsInfo[0].FieldsName( iField ).empty() )
            {
                aMeshData.FieldsInfo[0].FieldsName( iField ) = "genericFieldName_"+std::to_string( iField );
            }

            MORIS_ASSERT( aMeshData.FieldsInfo[0].FieldsRank( iField ) != EntityRank::INVALID, "Field rank was not provided.");
        }

        // Loop over the number of fields
        aMeshData.FieldsInfo[0].FieldsName.resize( mMaxNumFields );
        for ( uint iField = tNumFields; iField < mMaxNumFields; ++iField )
        {
            aMeshData.FieldsInfo[0].FieldsName( iField ) = "dummyField";
        }
    }
    // ----------------------------------------------------------------------------

    // Verifications of set data arrangement and consistency for meshes generated from data.
    void
    Mesh_STK::check_and_update_sets_data(
            MtkMeshData&   aMeshData )
    {
        ///////////////////////////
        // Checks for block sets //
        ///////////////////////////

        if ( aMeshData.SetsInfo[0].BlockSetsInfo != NULL )
        {
            uint tNumBlockSets  = aMeshData.SetsInfo[0].BlockSetsInfo[0].BSetInds[0].max() + 1;

            MORIS_ASSERT( aMeshData.SetsInfo[0].BlockSetsInfo[0].BSetInds[0].length() == aMeshData.ElemConn[0].n_rows(),
                          "Size of PartOwner vector should be the same as the number of elements." );

            // Communicate with other processors and see which one has the maximum
            uint tNumGlobalBlockSets = gather_value_and_bcast_max( tNumBlockSets );
            mSetRankFlags[2]         = true;
            mSetNames[2].resize( tNumGlobalBlockSets );

            // Loop over the number of block sets
            for ( uint iBSet = 0; iBSet < tNumGlobalBlockSets; ++iBSet )
            {
                // Check if set names were provided
                if ( aMeshData.SetsInfo[0].BlockSetsInfo[0].BSetNames( iBSet ).empty() )
                {
                    mSetNames[2][iBSet] = "BlockSet_" + std::to_string( iBSet );
                }
                else
                {
                    mSetNames[2][iBSet] = aMeshData.SetsInfo[0].BlockSetsInfo[0].BSetNames( iBSet );
                }
            }
        }
        else
        {
            // Create a block set tha contains the entire mesh by default
            mSetNames[2].resize( 1, "block_1" );
        }

        ///////////////////////////
        // Checks for side sets //
        ///////////////////////////

        if ( aMeshData.SetsInfo[0].SideSetsInfo != NULL )
        {
            uint tNumSideSets = aMeshData.SetsInfo[0].SideSetsInfo[0].ElemIdsAndSideOrds[0].size();
            mSetRankFlags[1]  = true;

            mSetNames[1].resize( tNumSideSets );

            // Loop over the number of block sets
            for ( uint iSSet = 0; iSSet < tNumSideSets; ++iSSet )
            {
                // Check if set names were provided
                if ( aMeshData.SetsInfo[0].SideSetsInfo[0].SSetNames( iSSet ).empty() )
                {
                    mSetNames[1][iSSet] = "SideSet_"+std::to_string( iSSet );
                }
                else
                {
                    mSetNames[1][iSSet] = aMeshData.SetsInfo[0].SideSetsInfo[0].SSetNames( iSSet );
                }

                // Check if side set specific info was provided
                std::string tTest = "Number of columns in side set should be equal to 2." ;
                MORIS_ASSERT( ( aMeshData.SetsInfo[0].SideSetsInfo[0].ElemIdsAndSideOrds[0]( iSSet ).n_cols() == 2 ) ||
                                    ( aMeshData.SetsInfo[0].SideSetsInfo[0].ElemIdsAndSideOrds[0]( iSSet ).n_cols() == 0 ) ,
                                "Number of columns in side set should be equal to 2. "
                                "The first column should have element Ids; and the second, side ordinals.");
            }
        }

        ///////////////////////////
        // Checks for node sets //
        ///////////////////////////

        if ( aMeshData.SetsInfo[0].NodeSetsInfo != NULL )
        {
            uint tNumNodeSets = aMeshData.SetsInfo[0].NodeSetsInfo[0].EntIds[0].size();
            mSetRankFlags[0]  = true;

            mSetNames[0].resize( tNumNodeSets );

            // Loop over the number of block sets
            for ( uint iNSet = 0; iNSet < tNumNodeSets; ++iNSet )
            {
                // Check if set names were provided
                if ( aMeshData.SetsInfo[0].NodeSetsInfo[0].NSetNames( iNSet ).empty() )
                {
                    mSetNames[0][iNSet] = "NodeSet_"+std::to_string( iNSet );
                }
                else
                {
                    mSetNames[0][iNSet] = aMeshData.SetsInfo[0].NodeSetsInfo[0].NSetNames( iNSet );
                }

                // Minimum check for node set
                MORIS_ASSERT( aMeshData.SetsInfo[0].SideSetsInfo[0].ElemIdsAndSideOrds[0]( iNSet ).n_rows() <= aMeshData.NodeCoords[0].n_rows(),
                              "Number of nodes in node set is greater than total number of nodes." );
            }
        }
    }

    // Main interface with STK that include calls to functions that provide specific implementation details.
    void
    Mesh_STK::build_mesh(
            MtkMeshData   aMeshData )
    {
        // A Mesh contains collections of entities, parts, fields, and field data. The STK Mesh API separates
        // these collections into 'MetaData' and 'BulkData'.
        //////////////////////////////////
        //   META DATA INITIALIZATION   //
        //////////////////////////////////

        // The MetaData component of a STK Mesh contains the definitions of its parts, the definitions of its
        // fields, and definitions of relationships among its parts and fields. For example, a subset relationship
        //  can be declared between two parts, and a field definition can be limited to specific parts.

        // Declare and initialize Stk mesh
        stk::mesh::MetaData * meshMeta = new stk::mesh::MetaData( mNumDims );

        // Set member variable as pointer to meta_data
        mMtkMeshMetaData = ( meshMeta );

        // Declare all additional parts provided by the user (default parts are always created by STK)
        this->declare_mesh_parts( aMeshData );

        // Declare all fields (including coordinates)
        this->declare_mesh_fields( aMeshData );

        // Commit MetaData before populating the BulkData
        mMtkMeshMetaData->commit();

        ////////////////////////////////
        // BULK DATA INITIALIZATION   //
        ////////////////////////////////

        // The BulkData component of a STK Mesh contains entities, entity ownership and ghosting
        // information, connectivity data, and field data. For efficiency, the BulkData API enables access to
        // data via buckets, in addition to via entity and rank.

        // Declare MPI communicator
        stk::ParallelMachine tPM = MPI_COMM_WORLD;

        // Declare aura
        stk::mesh::BulkData::AutomaticAuraOption aAutoAuraOption = stk::mesh::BulkData::AutomaticAuraOption::AUTO_AURA;

        // Create BulkData Object
        stk::mesh::BulkData * meshBulk = new stk::mesh::BulkData( *mMtkMeshMetaData, tPM, aAutoAuraOption );

        // Set member variable as pointer and bulk_data
        mMtkMeshBulkData = ( meshBulk );

        // Use STK IO to populate a STK Mesh
        MPI_Comm aCommunicator = MPI_COMM_WORLD;
        mMeshReader = new stk::io::StkMeshIoBroker( aCommunicator );

        // Create mesh database using the IO broker
        mMeshReader->set_bulk_data( *mMtkMeshBulkData );

        // Assign element to node connectivity
        this->populate_mesh_database( aMeshData );

        // Assign coordinates and any other field given by the user
        this->populate_mesh_fields( aMeshData );

        // Generate additional local to global maps (only for meshes generated from data).
        // Elemental and nodal information has been taken care of already in this case.
        if (aMeshData.CreateAllEdgesAndFaces )
        {
            this->create_additional_communication_lists_from_data();
        }
    }

    // ----------------------------------------------------------------------------

    // First declaration to structure the database before filling the data
    void
    Mesh_STK::declare_mesh_parts(
            MtkMeshData   aMeshData )
    {
        // Part is a general term for a subset of the entities in a mesh. STK Mesh automatically creates
        // four parts at startup: the universal part, the locally-owned part, the globally-shared part,
        // and the aura part. These parts are important to the basic understanding of ghosting. In addition,
        // Exodus parts, such as blocks, sidesets, and nodesets, are created if an Exodus file is read in.
        // Each entity in the mesh must be a member of one or more parts.
        uint tNumNodesPerElem = aMeshData.ElemConn[0].size( 1 );

        // Declare and initialize topology type. Also check if element type is supported
        stk::topology::topology_t tTopology = get_mesh_topology( mNumDims, tNumNodesPerElem );

        if ( aMeshData.SetsInfo != NULL ) // For all (block, node, side) sets
        {
            ////////////////////////
            // Declare block sets //
            ////////////////////////
            uint tNumBlockSets = mSetNames[2].size();

            for ( uint iSet = 0; iSet < tNumBlockSets; ++iSet )
            {
                // Declare part and add it to the IOBroker (needed for output).
                stk::mesh::Part& aSetPart = mMtkMeshMetaData->declare_part_with_topology( mSetNames[2][iSet], tTopology );
                // Add Part to the IOBroker (needed for output).
                stk::io::put_io_part_attribute( aSetPart );
            }

            ///////////////////////
            // Declare side sets //
            ///////////////////////
            if ( mSetRankFlags[1] )
            {
                uint tNumSideSets = aMeshData.SetsInfo[0].SideSetsInfo[0].ElemIdsAndSideOrds[0].size();

                for ( uint iSet = 0; iSet < tNumSideSets; ++iSet )
                {
                    // Declare part and add it to the IOBroker (needed for output).
                    stk::mesh::Part& aSetPart = mMtkMeshMetaData->declare_part( mSetNames[1][iSet], mMtkMeshMetaData->side_rank() );
                    // Add Part to the IOBroker (needed for output).
                    stk::io::put_io_part_attribute( aSetPart );
                }
            }

            ///////////////////////
            // Declare node sets //
            ///////////////////////
            if ( mSetRankFlags[0] )
            {
                uint tNumNodeSets = aMeshData.SetsInfo[0].NodeSetsInfo[0].EntIds[0].size();

                for ( uint iSet = 0; iSet < tNumNodeSets; ++iSet )
                {
                    // Declare part and add it to the IOBroker (needed for output).
                    stk::mesh::Part& aSetPart = mMtkMeshMetaData->declare_part( mSetNames[0][iSet], stk::topology::NODE_RANK );
                    // Add Part to the IOBroker (needed for output).
                    stk::io::put_io_part_attribute( aSetPart );
                }
            }
        }
        else
        {
            // Add default part if no block sets were provided
            stk::mesh::Part& tBlock = mMtkMeshMetaData->declare_part_with_topology( mSetNames[2][0], tTopology );
            // Add Part to the IOBroker (needed for output).
            stk::io::put_io_part_attribute( tBlock );
        }
    }

    // ----------------------------------------------------------------------------

    // Second declaration to structure the database before filling the data

    void
    Mesh_STK::declare_mesh_fields(
            MtkMeshData   aMeshData )
    {
        // Fields are data associated with mesh entities. Examples include coordinates, velocity,
        // displacement, and temperature. A field in STK Mesh can hold any data type (e.g., double or int)
        // and any number of scalars per entity (e.g., nodal velocity field has three doubles per node).
        // A field can be allocated (defined) on a whole mesh or on only a subset (part) of that mesh.
        // For example, a material property can be allocated on a specified element block.

        // Declare coordinates field
        Field3Comp* tCoord_field = &mMtkMeshMetaData->declare_field<Field3Comp>( stk::topology::NODE_RANK, "coordinates" );
        stk::mesh::put_field( *tCoord_field, mMtkMeshMetaData->universal_part() );

    //    Field3Comp* tCoord_field = mMtkMeshMetaData->declare_field< Field3Comp >( stk::topology::NODE_RANK, stk::io::CoordinateFieldName);
    //    stk::io::set_field_role( tCoord_field, Ioss::Field::MESH);
    //    mMtkMeshMetaData->set_coordinate_field( &tCoord_field );
    //    stk::mesh::put_field( *tCoord_field, mMtkMeshMetaData->universal_part() );

        // Declare all additional fields provided by the user
        if ( mFieldInDataGiven )
        {
            // WARNING: Currently hardcoded for 8 fields only
            MORIS_ASSERT( aMeshData.FieldsInfo[0].FieldsData[0].size() <= mMaxNumFields, "A maximum of 20 fields is currently supported");

            std::string tFieldNoData = "dummyField";

            if ( aMeshData.FieldsInfo[0].FieldsName( 0 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 0 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 1 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 1 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 2 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 2 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 3 ).compare( tFieldNoData ) != 0 )
            {
               this->internal_declare_mesh_field( aMeshData, 3 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 4 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 4 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 5 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 5 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 6 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 6 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 7 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 7 );
            }
            ///
            if ( aMeshData.FieldsInfo[0].FieldsName( 8 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 8 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 9 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 9 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 10 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 10 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 11 ).compare( tFieldNoData ) != 0 )
            {
               this->internal_declare_mesh_field( aMeshData, 11 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 12 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 12 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 13 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 13 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 14 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 14 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 15 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 15 );
            }
            ///
            if ( aMeshData.FieldsInfo[0].FieldsName( 16 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 16 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 17 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 17 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 18 ).compare( tFieldNoData ) != 0 )
            {
                this->internal_declare_mesh_field( aMeshData, 18 );
            }
            if ( aMeshData.FieldsInfo[0].FieldsName( 19 ).compare( tFieldNoData ) != 0 )
            {
               this->internal_declare_mesh_field( aMeshData, 19 );
            }
        }
    }


    // Declare size of a field (per entity) and throw an error if it is not supported
    void
    Mesh_STK::internal_declare_mesh_field(
            MtkMeshData   aMeshData,
            uint          iField )
    {
        // Get field variables
        uint tNumFieldComp     = aMeshData.FieldsInfo[0].FieldsData[0]( iField ).n_cols();
        std::string tFieldName = aMeshData.FieldsInfo[0].FieldsName( iField );
        EntityRank tFieldRank  = aMeshData.FieldsInfo[0].FieldsRank( iField );

        stk::mesh::EntityRank tStkFieldRank = this->get_stk_entity_rank( tFieldRank );
        stk::mesh::Selector aFieldPart      = mMtkMeshMetaData->universal_part();

        if ( aMeshData.FieldsInfo[0].SetsOwner != NULL )
        {
            if ( !aMeshData.FieldsInfo[0].SetsOwner[0]( iField ).empty() )
            {
                aFieldPart = *mMtkMeshMetaData->get_part( aMeshData.FieldsInfo[0].SetsOwner[0]( iField ) );
            }
        }

        switch ( tNumFieldComp )
        {
        case 1: // Scalar Field
        {
            // Declare fields
            mField1CompVec.push_back( & mMtkMeshMetaData->declare_field< Field1Comp >( tStkFieldRank, tFieldName ) );
            stk::mesh::put_field( *mField1CompVec.back(), aFieldPart, 1 );
            break;
        }
        case 2: // Vector Field with 2 components
        {
            // Declare fields
            mField2CompVec.push_back( & mMtkMeshMetaData->declare_field< Field2Comp >( tStkFieldRank, tFieldName ) );
            stk::mesh::put_field( *mField2CompVec.back(), aFieldPart );
            break;
        }
        case 3: // Vector Field with 3 components
        {
            // Declare fields
            mField3CompVec.push_back( & mMtkMeshMetaData->declare_field< Field3Comp >( tStkFieldRank, tFieldName ) );
            stk::mesh::put_field( *mField3CompVec.back(), aFieldPart );
            break;
        }
        case 4: // Tensor Field with 4 components
        {
            // Declare fields
            mField4CompVec.push_back( & mMtkMeshMetaData->declare_field< Field4Comp >( tStkFieldRank, tFieldName ) );
            stk::mesh::put_field( *mField4CompVec.back(), aFieldPart );
            break;
        }
        case 9: // Tensor Field with 9 components
        {
            // Declare fields
            mField9CompVec.push_back( & mMtkMeshMetaData->declare_field< Field9Comp >( tStkFieldRank, tFieldName ) );
            stk::mesh::put_field( *mField9CompVec.back(), aFieldPart );
            break;
        }
        default:
        {
            MORIS_ASSERT( 0, "Number of components (columns) for all FieldsData should match one of the supported sizes {1, 2, 3, 4, 9}." );
            break;
        }
        }
    }
    // ----------------------------------------------------------------------------

    // Add mesh information to database
    void
    Mesh_STK::populate_mesh_database(
            MtkMeshData   aMeshData )
    {
        ///////////////////////////////
        // Begin modification cycle  //
        ///////////////////////////////
        mMtkMeshBulkData->modification_begin();

        // Generate basic mesh information
        this->process_block_sets( aMeshData );

        // Declare node sets to mesh if they exist
        if ( mSetRankFlags[0] )
        {
            this->process_node_sets( aMeshData );
        }

        ///////////////////////////////
        // Close modification cycle  //
        ///////////////////////////////
        mMtkMeshBulkData->modification_end();

        // Declare node sets to mesh
        if ( mSetRankFlags[1] )
        {
            // If side sets were provided, generate only edges and/or faces of the
            // corresponding sets of the elements containing such entities. The elements
            // of the side sets entities will be moved to another part.
            this->process_side_sets( aMeshData );
        }
        else if ( aMeshData.CreateAllEdgesAndFaces )
        {    // If the user requires create all additional entities, use the functions below.
            // Note that this could potentially increase significantly memory usage and time for
            // creating the mesh for big amounts of data.

            stk::mesh::create_edges( *mMtkMeshBulkData );
            stk::mesh::create_faces( *mMtkMeshBulkData, true ); // Boolean to specify if want to connect faces to edges
        }
    }
// ----------------------------------------------------------------------------
    // Add all blocks information to database
    void
    Mesh_STK::process_block_sets(
            MtkMeshData   aMeshData )
    {
        // Get all sets provided by the user and go to the block set
        uint tNumElems     = aMeshData.ElemConn[0].size( 0 );
        uint tNumBlockSets = 1;

        Matrix< DDUMat >  aOwnerPartInds( tNumElems, 1, 0 );
        std::vector< stk::mesh::PartVector > aPartBlocks( 1 );

        // Update to number of blocks provided by the user
        if ( mSetRankFlags[2] )
        {
            tNumBlockSets  = aMeshData.SetsInfo[0].BlockSetsInfo[0].BSetInds[0].max() + 1;
            aOwnerPartInds = aMeshData.SetsInfo[0].BlockSetsInfo[0].BSetInds[0];
            aPartBlocks.resize( tNumBlockSets );
        }

        // Populate part blocks
        for ( uint iSet = 0; iSet < tNumBlockSets; ++iSet )
        {
            // Get block sets provided by user
            stk::mesh::Part* tBlock = mMtkMeshMetaData->get_part( mSetNames[2][iSet] );
            aPartBlocks[iSet]       = { tBlock };
        }

        // Declare MPI communicator
        stk::ParallelMachine tPM = MPI_COMM_WORLD;
        uint tParallelSize       = stk::parallel_machine_size( tPM );
        if ( tParallelSize == 1 )
        {
            // serial run
            this->populate_mesh_database_serial( aMeshData, aPartBlocks, aOwnerPartInds );
        }
        else
        {
            // Populating mesh is a bit more complicated in parallel because of entity sharing stuff
            this->populate_mesh_database_parallel( aMeshData, aPartBlocks, aOwnerPartInds );
        }
    }
// ----------------------------------------------------------------------------
    // Add all fields information to database
    void
    Mesh_STK::populate_mesh_fields(
            MtkMeshData   aMeshData )
    {
        // Get the coordinates field from Stk
        stk::mesh::FieldBase const* aCoord_field_i = mMtkMeshMetaData->coordinate_field();
        uint tNumNodes                             = aMeshData.NodeCoords[0].n_rows();

        // Loop over the number of nodes
        for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
        {
            // Get global Id of current node and create "node entity" for stk mesh
            uint aId                = aMeshData.LocaltoGlobalNodeMap[0]( iNode );
            stk::mesh::Entity aNode = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aId );

            // Store the coordinates of the current node
            if ( mMtkMeshBulkData->is_valid( aNode ) )
            {
                // Add coordinates information to the BulkData
                double* tCoord_data = static_cast <double*> ( stk::mesh::field_data ( *aCoord_field_i, aNode ) );
                for ( uint iDim = 0; iDim < mNumDims; ++iDim )
                {
                    tCoord_data[iDim] = aMeshData.NodeCoords[0]( iNode, iDim );
                }
            }
        }

        if ( mFieldInDataGiven )
        {
            // Get the number of fields
            uint tNumFields = aMeshData.FieldsInfo[0].FieldsData[0].size();
            std::vector< stk::mesh::FieldBase * > aFieldVector;

            // Loop over the number of fields
            for ( uint iField = 0; iField < tNumFields; ++iField )
            {
                // Get field variables
                uint tNumFieldComp     = aMeshData.FieldsInfo[0].FieldsData[0]( iField ).n_cols();
                std::string tFieldName = aMeshData.FieldsInfo[0].FieldsName( iField );
                EntityRank tFieldRank  = aMeshData.FieldsInfo[0].FieldsRank( iField ) ;
                uint tNumFieldEntities = aMeshData.FieldsInfo[0].FieldsData[0]( iField ).n_rows();

                stk::mesh::EntityRank tStkFieldRank  = this->get_stk_entity_rank( tFieldRank );

                // If set owner was provided, verify that the number of field values are the same as
                // the number of entities in the set (or if just one value was provided to populate the entire field).
                Matrix< IdMat >  tFieldIds;
                if ( aMeshData.FieldsInfo[0].SetsOwner != NULL )
                {
                    if ( !aMeshData.FieldsInfo[0].SetsOwner[0]( iField ).empty() )
                    {
                        tFieldIds = this->get_set_entity_ids( tStkFieldRank, aMeshData.FieldsInfo[0].SetsOwner[0]( iField ) );
                    }
                    else
                    {
                        tFieldIds = this->get_entities_owned_and_shared_by_current_proc( tFieldRank );
                    }

                    MORIS_ASSERT( ( tFieldIds.length() == tNumFieldEntities ) || ( tNumFieldEntities == 1 ),
                            "Field data should match the number of entities in set owner.");
                }
                else
                {
                    tFieldIds = this->get_entities_owned_and_shared_by_current_proc( tFieldRank );

                    MORIS_ASSERT( ( tFieldIds.length() == tNumFieldEntities ) || ( tNumFieldEntities == 1 ),
                            "Field data should match the number of entities in processor.");
                }

                // Update number of entities in field if only one value was provided
                if ( tNumFieldEntities == 1 )
                {
                    tNumFieldEntities = tFieldIds.length();
                }

                // Get field pointer
                stk::mesh::FieldBase * aFieldBase = mMtkMeshMetaData->get_field( tStkFieldRank, tFieldName );

                // Loop over field entities
                for ( uint iEntityInd = 0; iEntityInd < tNumFieldEntities; ++iEntityInd )
                {
                    // Get global Id of current entity based on rank (use map provided by the user).
                    // This is done only if the field is not contained in any set.
                    stk::mesh::Entity aEntity = mMtkMeshBulkData->get_entity( tStkFieldRank, tFieldIds( iEntityInd ) );

                    // Store the coordinates of the current entity
                    if ( mMtkMeshBulkData->is_valid( aEntity ) )
                    {
                        double* tFieldEntityData = static_cast <double*> ( stk::mesh::field_data ( *aFieldBase, aEntity ) );

                        // Add field data to the BulkData
                        for ( uint iComp = 0; iComp < tNumFieldComp; ++iComp )
                        {
                            if ( ( aMeshData.FieldsInfo[0].FieldsData[0]( iField ).n_rows() == 1 ) &&
                                 ( aMeshData.FieldsInfo[0].FieldsData[0]( iField ).n_cols() == 1 ) )
                            {
                                tFieldEntityData[iComp] = aMeshData.FieldsInfo[0].FieldsData[0]( iField )( 0, 0 );
                            }
                            else
                            {
                                tFieldEntityData[iComp] = aMeshData.FieldsInfo[0].FieldsData[0]( iField )( iEntityInd, iComp );
                            }
                        }
                    }
                }
                // Store field in vector
                aFieldVector.push_back( aFieldBase );
            }
            // do parallel stuff for syncing fields data
        }
    }
// ----------------------------------------------------------------------------
    // Provide element type (Hex8, Tri3, etc) and throw error if element is not supported yet.
    stk::topology::topology_t
    Mesh_STK::get_mesh_topology(
            uint   aModelDim,
            uint   aNumNodesInElem )
    {
        stk::topology::topology_t tTopology = stk::topology::INVALID_TOPOLOGY;

        // MTK supports the following 1D, 2D and 3D element topology temporarily

        if ( aModelDim == 1 ) // 1D
        {
            switch ( aNumNodesInElem )
            {
            case 2:
            {
                tTopology = stk::topology::LINE_2_1D;
                break;
            }
            case 3:
            {
                tTopology = stk::topology::LINE_3_1D;
                break;
            }
            default:
            {
                MORIS_ASSERT( 0, "MTK mesh build from data currently handles only LINE_2 for 1D elements.");
                break;
            }
            }
        }
        else if ( aModelDim == 2 ) // 2D
        {
            switch ( aNumNodesInElem )
            {
            case 3:
            {
                tTopology = stk::topology::TRI_3_2D;
                break;
            }
            case 4:
            {
                tTopology = stk::topology::QUAD_4_2D;
                break;
            }
            case 6:
            {
                tTopology = stk::topology::TRI_6_2D;
                break;
            }
            case 9:
            {
                tTopology = stk::topology::QUAD_9_2D;
                break;
            }
            default:
            {
                MORIS_ASSERT( 0, "MTK mesh build from data currently handles only TRI_3, TRI_6_2D, QUAD_4 and QUAD_9_2D for 2D elements.");
                break;
            }
            }
        }
        else if ( aModelDim == 3 ) // 3D
        {
            switch ( aNumNodesInElem )
            {
            case 4:
            {
                tTopology = stk::topology::TET_4;
                break;
            }
            case 8:
            {
                tTopology = stk::topology::HEX_8;
                break;
            }
            case 20:
            {
                tTopology = stk::topology::HEX_20;
                break;
            }
            case 27:
            {
                tTopology = stk::topology::HEX_27;
                break;
            }
            default:
            {
                MORIS_ASSERT( 0, "MTK mesh build from data currently handles only TET_4, HEX8, HEX_20 and HEX_27 for 3D elements.");
                break;
            }
            }
        }

        return tTopology;
    }
// ----------------------------------------------------------------------------
    // Function to create edges and faces communication lists in parallel for meshes generated from data
    void
    Mesh_STK::create_additional_communication_lists_from_data()
    {
        this->create_facets_communication_lists();
        this->create_owners_communication_lists();
        this->create_shared_communication_lists();
    }
// ----------------------------------------------------------------------------
    // Function to create edges and faces communication lists in parallel
    void
    Mesh_STK::create_facets_communication_lists()
    {
        // Create maps for edges if the problem is not 1D, and maps for faces if it is 3D.
        // NOTE: Not supporting 1D elements in 2D nor 3D space.
        if ( mNumDims > 1 )
        {
            uint tNumEdges = this->get_num_edges();

            // Access entities stored in mesh database
            Matrix< DDUMat >  tEdgeIds = this->get_entities_universal( EntityRank::EDGE );

            // resize member variable to its right size
            mEntityLocaltoGlobalMap(1).resize( tNumEdges, 1 );
            mEdgeMapToOwnerProc.resize( tNumEdges, 1 );

            // Populate internal member variable that contains the local index to
            // global id node communication information
            for ( uint iEdge = 0; iEdge < tNumEdges; ++iEdge )
            {
                // local to global and owner processor
                mEntityLocaltoGlobalMap(1)( iEdge ) = tEdgeIds( iEdge );
                mEdgeMapToOwnerProc( iEdge )   = this->get_entity_owner((moris_index) iEdge, EntityRank::EDGE );
            }
        }

        if ( mNumDims > 2 )
        {
            uint tNumFaces   = this->get_num_faces();

            // Access entities stored in mesh database
            Matrix< DDUMat >  tFaceIds = this->get_entities_universal( EntityRank::FACE );

            // resize member variable to its right size
            mEntityLocaltoGlobalMap(2).resize( tNumFaces, 1 );
            mFaceMapToOwnerProc.resize( tNumFaces, 1 );

            // Populate internal member variable that contains the local index to
            // global id element communication information
            for ( uint iFace = 0; iFace < tNumFaces; ++iFace )
            {
                // local to global and owner processor
                mEntityLocaltoGlobalMap(2)( iFace ) = tFaceIds( iFace );
                mFaceMapToOwnerProc( iFace )   = this->get_entity_owner( iFace, EntityRank::FACE );
            }
        }
    }
// ----------------------------------------------------------------------------
    // Function to create edges and faces communication lists in parallel
    void
    Mesh_STK::create_owners_communication_lists()
    {
        // Get basic mesh information
        uint tNumElems = this->get_num_elems();
        uint tNumNodes = this->get_num_nodes();

        // resize member variable to its right size
        mNodeMapToOwnerProc.resize( tNumNodes, 1 );

        // Populate internal member variable that contains the local index to
        // global id node communication information
        for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
        {
            // Owner processor
            mNodeMapToOwnerProc( iNode ) = this->get_entity_owner( iNode, EntityRank::NODE );
        }

        // resize member variable to its right size
        mElemMapToOwnerProc.resize( tNumElems, 1 );

        // Populate internal member variable that contains the local index to
        // global id element communication information
        for ( uint iElem = 0; iElem < tNumElems; ++iElem )
        {
            // Owner processor
            mElemMapToOwnerProc( iElem ) = this->get_entity_owner( iElem, EntityRank::ELEMENT );
        }
    }
// ----------------------------------------------------------------------------
    // Function to create edges and faces communication lists in parallel
    void
    Mesh_STK::create_shared_communication_lists()
    {
        // Get basic mesh information
        Matrix< DDUMat >  tNodesShared = this->get_entities_glb_shared_current_proc( EntityRank::NODE );
        uint tNumNodesShared      = tNodesShared.length();

        // Generate list of processors sharing information
        // -----------------------------------------------

        std::vector < uint > tActiveSharedProcs;

        // Loop over the number of nodes shared to get the shared processors
        for( uint iNodeShared = 0; iNodeShared < tNumNodesShared; ++iNodeShared )
        {
            Matrix< DDUMat >  tProcsSharing = this->get_procs_sharing_entity_by_id( tNodesShared( iNodeShared ), EntityRank::NODE );

            for( uint iProc = 0; iProc < tProcsSharing.length(); ++iProc )
            {
                tActiveSharedProcs.push_back( tProcsSharing( iProc ) );
            }
        }

        // Get processors shared excluding owner
        std::sort( tActiveSharedProcs.begin(), tActiveSharedProcs.end() );
        auto last = std::unique( tActiveSharedProcs.begin(), tActiveSharedProcs.end() );
        tActiveSharedProcs.erase( last, tActiveSharedProcs.end() );

        // remove current processor from list
        tActiveSharedProcs.erase( std::remove( tActiveSharedProcs.begin(), tActiveSharedProcs.end(), UINT_MAX ), tActiveSharedProcs.end() );

        // Populate sharing processors map
        uint aNumActiveSharedProcs = tActiveSharedProcs.size();
        for( uint iProcShared = 0; iProcShared < aNumActiveSharedProcs; ++iProcShared )
        {
            mProcsSharedToIndex.insert( std::pair< uint, uint > ( tActiveSharedProcs.at( iProcShared ), iProcShared ) );
        }

        // Generate nodes shared per processor list
        // ----------------------------------------
        mNodeMapToSharingProcs = this->get_shared_info_by_entity( aNumActiveSharedProcs, EntityRank::NODE );

        // Generate elements shared per processor list (because of aura)
        // -------------------------------------------
        mElemMapToSharingProcs = this->get_shared_info_by_entity( aNumActiveSharedProcs, EntityRank::ELEMENT );

        // Generate edges shared per processor list
        // ----------------------------------------
        mEdgeMapToSharingProcs = this->get_shared_info_by_entity( aNumActiveSharedProcs, EntityRank::EDGE );

        // Generate faces shared per processor list
        // ----------------------------------------
        mFaceMapToSharingProcs = this->get_shared_info_by_entity( aNumActiveSharedProcs, EntityRank::FACE );
    }
// ----------------------------------------------------------------------------
    // Add all node sets information to database
    void
    Mesh_STK::process_node_sets(
            MtkMeshData   aMeshData )
    {
        // Declare basic node set information
        uint tNumNodeSets = aMeshData.SetsInfo[0].NodeSetsInfo[0].EntIds[0].size();
        stk::mesh::EntityRank aStkSetRank  = stk::topology::NODE_RANK;

        for ( uint iSet = 0; iSet < tNumNodeSets; ++iSet )
        {
            // STK interface variables declaration
            stk::mesh::Part* aSetPart = mMtkMeshMetaData->get_part( mSetNames[0][iSet] );
            stk::mesh::PartVector aAddPart( 1, aSetPart );
            stk::mesh::EntityId aGlobalId;
            stk::mesh::Entity aEntity;

            // Populate node sets (change entity parts if nodes were declared already)
            uint tNumSetEntities = aMeshData.SetsInfo[0].NodeSetsInfo[0].EntIds[0]( iSet ).length();
            for ( uint iEntity = 0; iEntity < tNumSetEntities; ++iEntity )
            {
                // Declare new entity or add existing entity to declared part
                aGlobalId = aMeshData.SetsInfo[0].NodeSetsInfo[0].EntIds[0]( iSet )( iEntity, 0 );
                aEntity   = mMtkMeshBulkData->get_entity( aStkSetRank, aGlobalId );

                if ( !mMtkMeshBulkData->is_valid( aEntity ) )
                {
                    aEntity = mMtkMeshBulkData->declare_entity( aStkSetRank, aGlobalId, aAddPart );
                }
                else
                {
                    mMtkMeshBulkData->change_entity_parts( aEntity, aAddPart );
                }
            }// end of node sets declarations
        }
    }
// ----------------------------------------------------------------------------
    // Add all side sets information to database
    void
    Mesh_STK::process_side_sets(
            MtkMeshData   aMeshData )
    {
        // If the user requires create all additional entities, use the functions below.
        // Note that this could potentially increase significantly memory usage and time for
        // creating the mesh for big amounts of data.
        if ( aMeshData.CreateAllEdgesAndFaces )
        {
            stk::mesh::create_edges( *mMtkMeshBulkData );
            stk::mesh::create_faces( *mMtkMeshBulkData, true ); // Boolean to specify if want to connect faces to edges
        }

        // Get all sets provided by the user
        uint tNumSideSets = aMeshData.SetsInfo[0].SideSetsInfo[0].ElemIdsAndSideOrds[0].size();

        ///////////////////////////////
        // Begin modification cycle  //
        ///////////////////////////////
        mMtkMeshBulkData->modification_begin();

        for ( uint iSet = 0; iSet < tNumSideSets; ++iSet )
        {
            // STK interface variables declaration
            stk::mesh::Part* aSetPart = mMtkMeshMetaData->get_part( mSetNames[1][iSet] );
            stk::mesh::PartVector aAddPart( 1, aSetPart );
            stk::mesh::EntityId aGlobalElemId;
            stk::mesh::Entity aElemEntity;
            uint tRequestedSideOrd;

            uint tNumSetEntities = aMeshData.SetsInfo[0].SideSetsInfo[0].ElemIdsAndSideOrds[0]( iSet ).n_rows();
            for ( uint iEntity = 0; iEntity < tNumSetEntities; ++iEntity )
            {
                // First column contains element ids that will later be match with faces
                aGlobalElemId     = aMeshData.SetsInfo[0].SideSetsInfo[0].ElemIdsAndSideOrds[0]( iSet )( iEntity, 0 );
                aElemEntity       = mMtkMeshBulkData->get_entity( stk::topology::ELEMENT_RANK, aGlobalElemId );
                tRequestedSideOrd = aMeshData.SetsInfo[0].SideSetsInfo[0].ElemIdsAndSideOrds[0]( iSet )( iEntity, 1 );

                if ( !aMeshData.CreateAllEdgesAndFaces )
                {
                    // Create side entity
                    mMtkMeshBulkData->declare_element_side( aElemEntity, tRequestedSideOrd, aAddPart );
                }
                else
                {
                    // all faces and edges where created already. Only need to move entities to a
                    // different (previously declared) part assuming no new entities need to be created.
                    const stk::mesh::Entity * tSides                   = mMtkMeshBulkData->begin( aElemEntity, mMtkMeshMetaData->side_rank() );
                    const stk::mesh::ConnectivityOrdinal* tSideOrdinal = mMtkMeshBulkData->begin_ordinals( aElemEntity, mMtkMeshMetaData->side_rank() );
                    size_t tNumSides                                   = mMtkMeshBulkData->num_connectivity( aElemEntity, mMtkMeshMetaData->side_rank() );

                    for ( size_t sideI = 0 ; sideI < tNumSides ; ++sideI )
                    {
                        uint tCurrentSideOrd = static_cast<stk::mesh::ConnectivityOrdinal>( tSideOrdinal[sideI] );
                        if ( tRequestedSideOrd == tCurrentSideOrd )
                        {
                            // Move elements to side set part
                            mMtkMeshBulkData->change_entity_parts( tSides[sideI], aAddPart );
                            break;
                        }
                    }
                }
            }// end of side sets declarations
        }

        ///////////////////////////////
        // Close modification cycle  //
        ///////////////////////////////
        mMtkMeshBulkData->modification_end();
    }
// ----------------------------------------------------------------------------
    // Parallel specific implementation for blocks in database
    void
    Mesh_STK::populate_mesh_database_serial(
            MtkMeshData                            aMeshData,
            std::vector< stk::mesh::PartVector >   aElemParts,
            Matrix< DDUMat >                       aOwnerPartInds )
    {
        uint tNumElems        = aMeshData.ElemConn[0].size( 0 );
        uint aNumNodesPerElem = aMeshData.ElemConn[0].size( 1 );

        // Declare variables to access connectivity
        Matrix< IdMat >  tDummyMat( 1, aNumNodesPerElem );
        stk::mesh::EntityIdVector aCurrElemConn( aNumNodesPerElem );
        stk::mesh::EntityId aElemGlobalId;

        // Loop over the number of elements and interface between MORIS and Stk for connectivity
        for ( uint iElem = 0; iElem < tNumElems; ++iElem )
        {
            // Get row of nodes connected in moris variable and assign to STK variable
            tDummyMat.set_row( 0,aMeshData.ElemConn->get_row( iElem ));
            aCurrElemConn.assign( tDummyMat.data(), tDummyMat.data() + aNumNodesPerElem );

            // Declare element in function that also declares element-node relations internally
            aElemGlobalId = aMeshData.LocaltoGlobalElemMap[0]( iElem );
            stk::mesh::declare_element( *mMtkMeshBulkData, aElemParts[aOwnerPartInds( iElem )], aElemGlobalId, aCurrElemConn );
        }
    }
// ----------------------------------------------------------------------------
    // Parallel specific implementation for blocks in database
    void
    Mesh_STK::populate_mesh_database_parallel(
            MtkMeshData                          aMeshData,
            std::vector< stk::mesh::PartVector > aPartBlocks,
            Matrix< DDUMat >                     aOwnerPartInds )
    {
        // Mesh variables
        uint tNumNodes        = aMeshData.NodeCoords[0].n_rows();
        uint tNumElems        = aMeshData.ElemConn[0].size( 0 );
        uint aNumNodesPerElem = aMeshData.ElemConn[0].size( 1 );

        // Declare MPI communicator
        stk::ParallelMachine tPM = MPI_COMM_WORLD;
        uint tProcRank     = stk::parallel_machine_rank( tPM );

        // Check if the list provided is for nodes or elements shared.
        bool tNodesProcOwnerList = false;

        if ( aMeshData.EntProcOwner[0].n_rows() == tNumNodes )
        {
            tNodesProcOwnerList = true;
        }

        // Declare variables to access connectivity
        Matrix< IdMat >  tDummyMat( 1, aNumNodesPerElem );
        std::set<stk::mesh::EntityId> tNodesIDeclared;
        stk::mesh::EntityIdVector aCurrElemConn( aNumNodesPerElem );
        stk::mesh::EntityId aElemGlobalId;
        stk::mesh::Entity aNode;

        /////////////////////////////////
        // Using nodes proc owner list //
        /////////////////////////////////

        if ( tNodesProcOwnerList )
        {
            // Populate mesh with all elements since only locally owned were provided
            for ( uint iElem = 0; iElem < tNumElems; ++iElem )
            {
                // Get row of nodes connected in moris variable and assign to STK variable
                tDummyMat.get_row( 0 ) = aMeshData.ElemConn->get_row( iElem );
                aCurrElemConn.assign( tDummyMat.data(), tDummyMat.data() + aNumNodesPerElem );
                aElemGlobalId = aMeshData.LocaltoGlobalElemMap[0]( iElem );

                // declare element in function that also declares element-node relations internally
                stk::mesh::declare_element( *mMtkMeshBulkData, aPartBlocks[aOwnerPartInds( iElem )], aElemGlobalId, aCurrElemConn );
            }

            // Add directly nodes shared from shared list.
            for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
            {
                //Declare hanging/floating nodes, which do not relate to any element
                aNode   = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aMeshData.LocaltoGlobalNodeMap[0]( iNode, 0 ) );
                if ( !mMtkMeshBulkData->is_valid( aNode ) )
                {
                    aNode = mMtkMeshBulkData->declare_entity( stk::topology::NODE_RANK, aMeshData.LocaltoGlobalNodeMap[0]( iNode, 0 ) );
                }
                if ( aMeshData.EntProcOwner[0]( iNode ) != (moris_id)tProcRank )
                {
                    // Add node if it has not been assign to the mesh yet.
                    stk::mesh::Entity aNode = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aMeshData.LocaltoGlobalNodeMap[0]( iNode, 0 ) );
                    mMtkMeshBulkData->add_node_sharing( aNode, aMeshData.EntProcOwner[0]( iNode ) );
                }
            }
        }

        /////////////////////////////////
        // Using elems proc owner list //
        /////////////////////////////////

        else
        {
            stk::mesh::EntityId aNodeGlobalId;
            // Loop over the number of elements
            for ( uint iElem = 0; iElem < tNumElems; ++iElem )
            {
                // Add element to mesh in the corresponding part
                if ( aMeshData.EntProcOwner[0]( iElem ) == (moris_id)tProcRank )
                {
                    // Get row of nodes connected in moris variable and assign to STK variable
                    tDummyMat.get_row( 0 ) = aMeshData.ElemConn->get_row( iElem );
                    aCurrElemConn.assign( tDummyMat.data(), tDummyMat.data() + aNumNodesPerElem );

                    // Loop over the number of nodes in each element
                    for ( uint iNode = 0; iNode < aNumNodesPerElem; ++iNode )
                    {
                        // Get global Id of current node and create "node entity" for stk mesh
                        aNodeGlobalId = aMeshData.ElemConn[0]( iElem, iNode );
                        aNode = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aNodeGlobalId );

                        // Add node if it has not been assign to the mesh yet.
                        if ( !mMtkMeshBulkData->is_valid( aNode ) )
                        {
                            tNodesIDeclared.insert( aNodeGlobalId );
                        }
                    }

                    // declare element in function that also declares element-node relations internally
                    aElemGlobalId = aMeshData.LocaltoGlobalElemMap[0]( iElem );
                    stk::mesh::declare_element( *mMtkMeshBulkData, aPartBlocks[aOwnerPartInds( iElem )], aElemGlobalId, aCurrElemConn );
                }
            }

            // NOTE: This implementation requires the incorporation of the elements located at the boundaries of the processors
            // (i.e., aura elements) to the connectivity table and the elements processor owner list.

            // Loop over the number of elements
            for ( uint iElem = 0; iElem < tNumElems; ++iElem )
            {
                // Check if the element is not own by this processor
                if ( aMeshData.EntProcOwner[0]( iElem )!= (moris_id)tProcRank )
                {
                    // Loop over the number of nodes in each element
                    for ( uint iNode = 0; iNode < aNumNodesPerElem; ++iNode )
                    {
                        aNodeGlobalId = aMeshData.ElemConn[0]( iElem, iNode );

                        if ( tNodesIDeclared.find( aNodeGlobalId ) != tNodesIDeclared.end() )
                        {
                            // Add node if it has not been assign to the mesh yet.
                            aNode = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aNodeGlobalId );
                            mMtkMeshBulkData->add_node_sharing( aNode, aMeshData.EntProcOwner[0]( iElem ) );
                        }
                    }
                }
            }
        }
    }
// ----------------------------------------------------------------------------
    // Access set entity ids
    moris::Matrix< IdMat >
    Mesh_STK::get_set_entity_ids(
            stk::mesh::EntityRank   aEntityRank,
            std::string             aSetName ) const
    {
        // Get pointer to field defined by input name
        stk::mesh::Part* const tSetPart = mMtkMeshMetaData->get_part( aSetName );

        MORIS_ASSERT( tSetPart != NULL, "Set not found. Double check name provided." );

        // Access data through a selector
        stk::mesh::Selector tSetSelector( *tSetPart );
        stk::mesh::EntityVector aEntities;
        stk::mesh::get_selected_entities( tSetSelector, mMtkMeshBulkData->buckets( aEntityRank ), aEntities );

        // Get entity Ids
        uint tNumEntities = aEntities.size();
        Matrix< IdMat >  tOutputEntityIds ( tNumEntities, 1 );
        for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
        {
            tOutputEntityIds( iEntity ) = ( moris_id ) mMtkMeshBulkData->identifier( aEntities[iEntity] );
        }

        return tOutputEntityIds;
    }
// ----------------------------------------------------------------------------
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


