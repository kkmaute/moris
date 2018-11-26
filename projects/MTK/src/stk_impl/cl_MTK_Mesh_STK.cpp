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

        // Generate MetaData and Bulk Data instances (later to be pointed to member variables)
        stk::mesh::MetaData * meshMeta = new stk::mesh::MetaData;
        stk::mesh::BulkData * meshBulk = new stk::mesh::BulkData( *meshMeta, aCommunicator, this->get_aura_option() );

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


        stk::mesh::Field<real> & tNodeField = mMtkMeshMetaData->declare_field<stk::mesh::Field<real> >(stk::topology::NODE_RANK, "lsf1", 1);

        stk::mesh::put_field_on_entire_mesh(tNodeField);

        stk::io::set_field_role(tNodeField, Ioss::Field::TRANSIENT);

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
        // Initialize global to local map
        mEntityGlobaltoLocalMap = moris::Cell<std::unordered_map<moris_id,moris_index>>(4);
        setup_entity_global_to_local_map(EntityRank::NODE);
        setup_entity_global_to_local_map(EntityRank::EDGE);
        setup_entity_global_to_local_map(EntityRank::FACE);
        setup_entity_global_to_local_map(EntityRank::ELEMENT);

        set_up_vertices_and_cell();
    }

    // ----------------------------------------------------------------------------

    Mesh_STK::Mesh_STK(
            MtkMeshData &  aMeshData ):
        mEntityLocaltoGlobalMap(4),
        mEntityGlobaltoLocalMap(4),
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
            moris_id tId = mMtkMeshBulkData->identifier(tEntitiesConnected[i]);
            tLocalIndices(0, i) = get_loc_entity_ind_from_entity_glb_id(tId,aOutputEntityRank);
        }
        return tLocalIndices;
    }

    Matrix< IndexMat >
    Mesh_STK::get_elements_connected_to_element_and_face_ord_loc_inds(moris_index aElementIndex) const
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
                    tElemsConnectedToElem( 1,tCounter ) = faceIt;
                    tCounter++;
                }
            }

            if ( tDummyConnectivity.size()  > 1 )
            {
                if ( mMtkMeshBulkData->identifier(tDummyConnectivity[1]) !=
                     mMtkMeshBulkData->identifier(  tStkEntity ) )
                {
                    tElemsConnectedToElem( 0,tCounter ) = (moris_index) mMtkMeshBulkData->local_id(tDummyConnectivity[1]);
                    tElemsConnectedToElem( 1,tCounter ) = faceIt;
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

    Matrix< IndexMat >
    Mesh_STK::get_elements_connected_to_element_and_face_ind_loc_inds(moris_index aElementIndex) const
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
                     mMtkMeshBulkData->identifier( tStkEntity ) )
                {
                    tElemsConnectedToElem( 0,tCounter ) = (moris_index) mMtkMeshBulkData->local_id(tDummyConnectivity[0]);
                    tElemsConnectedToElem( 1,tCounter ) =  (moris_index) mMtkMeshBulkData->local_id(tFacesInElem[faceIt]);;
                    tCounter++;
                }
            }

            if ( tDummyConnectivity.size()  > 1 )
            {
                if ( mMtkMeshBulkData->identifier(tDummyConnectivity[1]) !=
                     mMtkMeshBulkData->identifier(  tStkEntity ) )
                {
                    tElemsConnectedToElem( 0,tCounter ) = (moris_index) mMtkMeshBulkData->local_id(tDummyConnectivity[1]);
                    tElemsConnectedToElem( 1,tCounter ) =  (moris_index) mMtkMeshBulkData->local_id(tFacesInElem[faceIt]);;
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
       return mEntityLocaltoGlobalMap((uint)aEntityRank)(aEntityIndex);
    }

    moris_index
    Mesh_STK::get_loc_entity_ind_from_entity_glb_id(moris_id        aEntityId,
                                                    enum EntityRank aEntityRank) const
    {

        auto tIter = mEntityGlobaltoLocalMap((uint)aEntityRank).find(aEntityId);
        MORIS_ERROR(tIter!=mEntityGlobaltoLocalMap((uint)aEntityRank).end(), "Provided Entity Id is not in the map, Has the map been initialized?");

        return tIter->second;
    }


    // ----------------------------------------------------------------------------

    Matrix<IdMat>
    Mesh_STK::get_entity_connected_to_entity_glob_ids(
            moris_id     aEntityId,
            enum EntityRank aInputEntityRank,
            enum EntityRank aOutputEntityRank) const
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

    Matrix< IndexMat >
    Mesh_STK::get_set_entity_loc_inds( enum EntityRank aSetEntityRank,
                                       std::string     aSetName) const
    {
        // Get pointer to field defined by input name
        stk::mesh::Part* const tSetPart = mMtkMeshMetaData->get_part( aSetName );

        MORIS_ASSERT( tSetPart != NULL, "Set not found. Double check name provided." );

        // Access data through a selector
        stk::mesh::Selector tSetSelector( *tSetPart );
        stk::mesh::EntityVector aEntities;
        stk::mesh::get_selected_entities( tSetSelector, mMtkMeshBulkData->buckets( this->get_stk_entity_rank( aSetEntityRank ) ), aEntities );

        // Get entity Ids
        uint tNumEntities = aEntities.size();
        Matrix< IndexMat >  tOutputEntityInds ( tNumEntities, 1 );
        for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
        {
            tOutputEntityInds( iEntity ) = this->get_loc_entity_ind_from_entity_glb_id(( moris_id ) mMtkMeshBulkData->identifier( aEntities[iEntity] ),aSetEntityRank);
        }

        return tOutputEntityInds;
    }

    // ----------------------------------------------------------------------------

    Matrix< DDRMat >
    Mesh_STK::get_entity_field_value_real_scalar(Matrix< IndexMat > const & aEntityIndices,
                                                 std::string        const & aFieldName,
                                                 enum EntityRank            aFieldEntityRank) const
    {
        //MORIS_ASSERT(aFieldEntityRank==EntityRank::NODE,"Only implemented for nodal scalar field");

        // Initialize Output
        size_t tNumEntities = aEntityIndices.n_cols();
        Matrix< DDRMat > tFieldValues(1,tNumEntities);

        // Get field by name and entity rank
        stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank(aFieldEntityRank);
        stk::mesh::Field<real> * tField   = mMtkMeshMetaData->get_field<stk::mesh::Field<real>>(tEntityRank,aFieldName);

        // make sure that field actually exists
        if( tField == NULL )
        {
            // select specifier for rank
            std::string tRank;

            switch( aFieldEntityRank )
            {
                case( moris::EntityRank::NODE ) :
                {
                    tRank = " node ";
                    break;
                }
                case( moris::EntityRank::ELEMENT ) :
                {
                    tRank = " element ";
                    break;
                }
                default :
                {
                    tRank = " ";
                    break;
                }
            }

            // assemble error message
            std::string tError = "Could not find" + tRank + "field " + aFieldName;

            // throw error
            MORIS_ERROR( tField != NULL, tError.c_str() );

        }
        // Loop over entities and access field value
        for (size_t i = 0; i < tNumEntities; i++ )
        {
            moris_id tId = get_glb_entity_id_from_entity_loc_index(aEntityIndices(0,i),aFieldEntityRank);
            stk::mesh::Entity tEntity = mMtkMeshBulkData->get_entity(tEntityRank, tId);

            // Store the coordinates of the current node
            real* tFieldData = stk::mesh::field_data ( *tField, tEntity );
            tFieldValues(0,i) = tFieldData[0];
        }

        return tFieldValues;
    }

    // ----------------------------------------------------------------------------

    void
    Mesh_STK::add_mesh_field_real_scalar_data_loc_inds(std::string      const & aFieldName,
                                                       enum EntityRank  const & aFieldEntityRank,
                                                       Matrix< DDRMat > const & aFieldData)
    {

        MORIS_ASSERT(aFieldEntityRank==EntityRank::NODE,"Only implemented for nodal scalar field");

        // Write Data to Field
        size_t tNumEntities = get_num_entities(aFieldEntityRank);

        // Get Field
        stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank(aFieldEntityRank);
        stk::mesh::Field<real> * tField = mMtkMeshMetaData->get_field<stk::mesh::Field<real>>(tEntityRank,aFieldName);
        for (size_t i = 0; i < tNumEntities; i++ )
        {
            // Get global Id of current node and create "node entity" for stk mesh
            //stk::mesh::EntityId nodeGlobalId = node_i;
            moris_id tId = get_glb_entity_id_from_entity_loc_index(i,aFieldEntityRank);
            stk::mesh::Entity tEntity = mMtkMeshBulkData->get_entity(tEntityRank, tId);

            // Store the coordinates of the current node
            real* tFieldData = stk::mesh::field_data ( *tField, tEntity );
            tFieldData[0] = aFieldData(i);
        }
    }


    // ----------------------------------------------------------------------------

//##############################################
// moris::Cell and Vertex Pointer Functions
//##############################################

    mtk::Cell &
    Mesh_STK::get_mtk_cell(moris_index aCellIndex)
    {
        return mMtkCells(aCellIndex);
    }

    // ----------------------------------------------------------------------------

    mtk::Vertex &
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
                mMeshReader->begin_output_step( outputFileIdx, mTimeStamp );
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
    moris::Matrix< IdMat >
    Mesh_STK::get_entities_in_selector_interface_glob_id(
            EntityRank            aRequestedEntityRank,
            stk::mesh::Selector   aSelectedEntities ) const
    {
        // Get selected entities
        std::vector<stk::mesh::Entity>tOutputEntityIDs;
        stk::mesh::get_selected_entities(
                aSelectedEntities, mMtkMeshBulkData->buckets( this->get_stk_entity_rank( aRequestedEntityRank ) ), tOutputEntityIDs );

        // Interface with STK, get the Ids and return them
        uint tNumEntity = tOutputEntityIDs.size();
        Matrix< IdMat >  tOutputEntityIDMat( tNumEntity, 1, 0 );

        for ( uint i = 0; i<tNumEntity; i++)
        {
            tOutputEntityIDMat( i, 0 ) = ( moris_id ) mMtkMeshBulkData->identifier( tOutputEntityIDs[i] );
        }

        return tOutputEntityIDMat;
    }

    // ----------------------------------------------------------------------------

    // return processors sharing a particular entity
    moris::Matrix< IdMat >
    Mesh_STK::get_procs_sharing_entity_by_id(
            moris_id              aEntityID,
            enum EntityRank   aEntityRank ) const
    {
        // Initialize returning mat with UINT_MAX in case no processors are sharing the given entity
        Matrix< IdMat >  tSharedProcsMat( 1, 1, INT_MAX );

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
            Matrix< IdMat >  tEntitiesShared = this->get_entities_glb_shared_current_proc( aEntityRank );
            uint tNumEntitiesShared = tEntitiesShared.length();
            moris_id tParallelRank  = par_rank();

            moris::Cell < moris::Cell < uint > > tTemporaryEntityMapSharingProcs( aNumActiveSharedProcs );

            // Loop over the number of nodes shared to get the shared processors
            for( uint iElemShared = 0; iElemShared < tNumEntitiesShared; ++iElemShared )
            {
                Matrix< IdMat >  tProcsSharing = this->get_procs_sharing_entity_by_id( tEntitiesShared( iElemShared ), aEntityRank );
                for( uint iProc = 0; iProc < tProcsSharing.length(); ++iProc )
                {
                    if( tProcsSharing(iProc) != (moris_id)tParallelRank )
                    {
                        tTemporaryEntityMapSharingProcs( mProcsSharedToIndex[ tProcsSharing( iProc ) ] ).push_back( tEntitiesShared( iElemShared ) );
                    }
                }
            }

            return tTemporaryEntityMapSharingProcs;
    }

    // ----------------------------------------------------------------------------

    void
    Mesh_STK::get_processors_whom_share_entity(moris_index       aEntityIndex,
                                               enum EntityRank   aEntityRank,
                                               Matrix< IdMat > & aProcsWhomShareEntity) const
    {
        // Convert index to ID
        stk::mesh::EntityId tEntityId = { (stk::mesh::EntityId)this->get_glb_entity_id_from_entity_loc_index(aEntityIndex, aEntityRank) };

        //Get entity
        stk::mesh::Entity tEntity = mMtkMeshBulkData->get_entity(get_stk_entity_rank(aEntityRank), tEntityId);

        // Intialize shared procs
        std::vector<int> tSharedProcs;

        // get shared processor IDs
        mMtkMeshBulkData->comm_procs(mMtkMeshBulkData->entity_key(tEntity), tSharedProcs);

        if (tSharedProcs.size() == 0)
        {
            tSharedProcs.push_back(mMtkMeshBulkData->parallel_owner_rank(tEntity));
        }

        // Initialize output
        aProcsWhomShareEntity.resize(1, tSharedProcs.size());

        // Cell to vector conversion
        for (uint i = 0; i < tSharedProcs.size(); i++)
        {
            aProcsWhomShareEntity(0, i) = tSharedProcs[i];
        }
    }

    // ----------------------------------------------------------------------------


    uint
    Mesh_STK::get_num_of_entities_shared_with_processor(moris_id        aProcessorRank,
                                                        enum EntityRank aEntityRank,
                                                        bool aSendFlag) const
    {
        uint tNumEntitiesInCommList = 0;
        if(aSendFlag)
        {
            tNumEntitiesInCommList = mEntitySendList((uint)aEntityRank)(aProcessorRank).n_cols();
        }
        else
        {
            tNumEntitiesInCommList = mEntityReceiveList((uint)aEntityRank)(aProcessorRank).n_cols();
        }

        return tNumEntitiesInCommList;
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

        // Parallel initialization
        uint tProcSize = par_size();
        // Access mesh data from struc and check input values for indispensable arguments
        MORIS_ASSERT( aMeshData.SpatialDim != NULL, "Number of spatial dimensions was not provided." );
        MORIS_ASSERT( aMeshData.ElemConn(0)!= NULL, "Element connectivity was not provided.");
        MORIS_ASSERT( aMeshData.NodeCoords != NULL, "Node coordinates were not provided." );

        // Initialize number of dimensions
        mNumDims = aMeshData.SpatialDim[0];

        if ( ( aMeshData.NodeProcOwner == NULL ) && ( tProcSize == 1 ) )
        {
            // Do nothing. This information is not needed in serial
        }
        else if ( aMeshData.NodeProcOwner != NULL )
        {
            // Verify sizes
            MORIS_ASSERT( ( aMeshData.NodeProcOwner[0].numel() == tNumNodes ) ||
                          ( aMeshData.NodeProcOwner[0].numel() == aMeshData.get_num_elements() ),
                          "Number of rows for EntProcOwner should match number of nodes or elements." );
        }
//        else
//        {
//            MORIS_ASSERT( 0, "Elements or nodes processor owner list must be provided in parallel." );
//        }

        // If nodal map was not provided and the simulation is run with 1 processor,
        // just generate the map. The number of ids provided in the map should match
        // the number of coordinates given for each processor.
        if ( aMeshData.LocaltoGlobalNodeMap != NULL )
        {
            // Verify sizes
            MORIS_ASSERT( aMeshData.LocaltoGlobalNodeMap[0].numel() == tNumNodes, "Number of rows for LocaltoGlobalNodeMap should match number nodes." );

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
        MORIS_ASSERT(aMeshData.LocaltoGlobalElemMap(0) != NULL, " No Local to global element map provided");
        // Verify sizes
        MORIS_ASSERT( aMeshData.size_local_to_global_elem_map() == aMeshData.get_num_elements(),
                      "Number of rows for LocaltoGlobalElemMap should match number of elements provided in connectivity map" );

        if(aMeshData.LocaltoGlobalElemMap.size() == 0)
        {
            mEntityLocaltoGlobalMap(3) =(* aMeshData.LocaltoGlobalElemMap(0) );
        }

        else
        {
            mEntityLocaltoGlobalMap(3) = aMeshData.collapse_element_map();
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
//        mFieldInDataGiven = true;
//
//        // Get the number of fields
//        uint tNumFields = aMeshData.FieldsInfo->get_num_fields();
//
//        // Verify that all field ranks were given
//        MORIS_ASSERT( aMeshData.FieldsInfo[0].FieldsRank.size() == tNumFields, "Number of field ranks should be the same as number of field data Mats." );
//
//        // Check if set owner names were provided
//        if ( aMeshData.FieldsInfo[0].SetsOwner != NULL )
//        {
//            MORIS_ASSERT( aMeshData.FieldsInfo[0].SetsOwner[0].size() == tNumFields ,
//                    "Set owner container should have names for all fields declared. "
//                    "If field is declared over universal part, provide empty string.");
//        }
//
//        // Loop over the number of fields
//        for ( uint iField = 0; iField < tNumFields; ++iField )
//        {
//            // Verify that field sizes (number of columns) match the ones suppported
//            Matrix< DDUMat >  tSupFieldComp = { {1, 2, 3, 4, 9} };
//            uint tNumFieldComp = aMeshData.FieldsInfo[0].FieldsData[0]( iField ).n_cols();
//            Matrix< DDBMat >  tDummy = ( tSupFieldComp == tNumFieldComp );
//            Matrix< DDNIMat >  tCompFound = find ( tDummy );
//
//
//            MORIS_ASSERT( !isempty( tCompFound ),
//                    "Number of components (columns) for all FieldsData should "
//                    "match one of the supported sizes {1, 2, 3, 4, 9}.");
//
//            // Check if field names were provided
//            if ( aMeshData.FieldsInfo[0].FieldsName( iField ).empty() )
//            {
//                aMeshData.FieldsInfo[0].FieldsName( iField ) = "genericFieldName_"+std::to_string( iField );
//            }
//
//            MORIS_ASSERT( aMeshData.FieldsInfo[0].FieldsRank( iField ) != EntityRank::INVALID, "Field rank was not provided.");
//        }
//
//        // Loop over the number of fields
//        aMeshData.FieldsInfo[0].FieldsName.resize( mMaxNumFields );
//        for ( uint iField = tNumFields; iField < mMaxNumFields; ++iField )
//        {
//            aMeshData.FieldsInfo[0].FieldsName( iField ) = "dummyField";
//        }
    }
    // ----------------------------------------------------------------------------

    // Verifications of set data arrangement and consistency for meshes generated from data.
    void
    Mesh_STK::check_and_update_sets_data(
            MtkMeshData &   aMeshData )
    {
        ///////////////////////////
        // Checks for block sets //
        ///////////////////////////
        if(aMeshData.has_mesh_sets())
        {
            if ( aMeshData.SetsInfo->get_num_block_sets() != 0 )
            {
                uint tNumBlockSets  = aMeshData.SetsInfo->get_num_block_sets();

                // Communicate with other processors and see which one has the maximum

                // TODO: It is not clear why this is needed, not only would
                // TODO: I need the proc with the maximum block sets but also the names?
                uint tNumGlobalBlockSets = gather_value_and_bcast_max( tNumBlockSets );
                mSetRankFlags[2]         = true;
                mSetNames[2].resize( tNumGlobalBlockSets );

                // Loop over the number of block sets
                for ( uint iBSet = 0; iBSet < tNumBlockSets; ++iBSet )
                {
                    // Get the block set
                    MtkBlockSetInfo* tBlockSet = aMeshData.SetsInfo->get_block_set(iBSet);

                    // Check if set names were provided
                    if ( !tBlockSet->blockset_has_name() )
                    {
                        mSetNames[2][iBSet]      = "BlockSet_" + std::to_string( iBSet );
                        tBlockSet->mBlockSetName = "BlockSet_" + std::to_string( iBSet );
                    }
                    else
                    {
                        mSetNames[2][iBSet] = tBlockSet->mBlockSetName;
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
            uint tNumSideSets = aMeshData.SetsInfo->get_num_side_sets();

            if ( tNumSideSets != 0 )
            {
                mSetRankFlags[1]  = true;

                mSetNames[1].resize( tNumSideSets );

                // Loop over the number of block sets
                for ( uint iSSet = 0; iSSet < tNumSideSets; ++iSSet )
                {
                    // Get the side set
                    MtkSideSetInfo* tSideSet = aMeshData.SetsInfo->get_side_set(iSSet);

                    // Check if set names were provided
                    if ( !tSideSet->sideset_has_name() )
                    {
                        mSetNames[1][iSSet] = "SideSet_"+std::to_string( iSSet );
                    }
                    else
                    {
                        mSetNames[1][iSSet] = tSideSet->mSideSetName;
                    }

                    // Check if side set specific info was provided
                    std::string tTest = "Number of columns in side set should be equal to 2." ;
                    MORIS_ASSERT( ( tSideSet->mElemIdsAndSideOrds->n_cols() == 2 ) ||
                                  ( tSideSet->mElemIdsAndSideOrds->n_cols() == 0 ) ,
                                  "Number of columns in side set should be equal to 2. "
                                  "The first column should have element Ids; and the second, side ordinals.");
                }
            }

            ///////////////////////////
            // Checks for node sets //
            ///////////////////////////
            uint tNumNodeSets = aMeshData.SetsInfo->get_num_node_sets();

            if ( tNumNodeSets != 0 )
            {
                mSetRankFlags[0]  = true;

                mSetNames[0].resize( tNumNodeSets );

                // Loop over the number of block sets
                for ( uint iNSet = 0; iNSet < tNumNodeSets; ++iNSet )
                {
                    // Get node set
                    MtkNodeSetInfo* tNodeSet = aMeshData.SetsInfo->get_node_set(iNSet);


                    // Check if set names were provided
                    if ( !tNodeSet->nodeset_has_name() )
                    {
                        mSetNames[0][iNSet] = "NodeSet_"+std::to_string( iNSet );
                    }
                    else
                    {
                        mSetNames[0][iNSet] = tNodeSet->mNodeSetName;
                    }

                }
            }
        }
    }

    // Main interface with STK that include calls to functions that provide specific implementation details.
    void
    Mesh_STK::build_mesh(
            MtkMeshData &   aMeshData )
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

        // Create BulkData Object
        stk::mesh::BulkData * meshBulk = new stk::mesh::BulkData( *mMtkMeshMetaData, tPM, this->get_aura_option() );

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
            setup_entity_global_to_local_map(EntityRank::FACE);
            setup_entity_global_to_local_map(EntityRank::EDGE);
        }


        // set timestamp
        mTimeStamp = aMeshData.TimeStamp;

        // set auto aura option
        mAutoAuraOption = aMeshData.AutoAuraOptionInSTK;
    }

    // ----------------------------------------------------------------------------

    // First declaration to structure the database before filling the data
    void
    Mesh_STK::declare_mesh_parts(
            MtkMeshData &  aMeshData )
    {
        // Part is a general term for a subset of the entities in a mesh. STK Mesh automatically creates
        // four parts at startup: the universal part, the locally-owned part, the globally-shared part,
        // and the aura part. In addition, Exodus parts, such as blocks, sidesets, and nodesets, are
        // created if an Exodus file is read in. Each entity in the mesh must be a member of one or more parts.

        uint tNumElementTypes = aMeshData.ElemConn.size();

        // Loop over the different element types and declare a part
        for(uint iET = 0; iET < tNumElementTypes; iET++)
        {
            uint tNumNodesPerElem = aMeshData.ElemConn(iET)->n_cols();
            // Declare and initialize topology type. Also check if element type is supported
            stk::topology::topology_t tTopology = get_mesh_topology( mNumDims, tNumNodesPerElem );

            // Add default part if no block sets were provided
            stk::mesh::Part& tBlock = mMtkMeshMetaData->declare_part_with_topology( "noblock_"+std::to_string(iET), tTopology );
            // Add Part to the IOBroker (needed for output).
            stk::io::put_io_part_attribute( tBlock );
        }



        if ( aMeshData.SetsInfo != NULL ) // For all (block, node, side) sets
        {
            ////////////////////////
            // Declare block sets //
            ////////////////////////
            uint tNumBlockSets = aMeshData.SetsInfo->get_num_block_sets();
            for ( uint iSet = 0; iSet < tNumBlockSets; ++iSet )
            {

                MtkBlockSetInfo* tBlockSet = aMeshData.SetsInfo->get_block_set(iSet);
                // Declare part and add it to the IOBroker (needed for output).
                stk::mesh::Part& aSetPart = mMtkMeshMetaData->declare_part_with_topology( tBlockSet->mBlockSetName,
                                                                                          get_stk_topo(tBlockSet->mBlockSetTopo) );
                // Add Part to the IOBroker (needed for output)
                stk::io::put_io_part_attribute( aSetPart );
            }

            ///////////////////////
            // Declare side sets //
            ///////////////////////
            if ( mSetRankFlags[1] )
            {
                uint tNumSideSets = aMeshData.SetsInfo->get_num_side_sets();

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
                uint tNumNodeSets = aMeshData.SetsInfo->get_num_node_sets();

                for ( uint iSet = 0; iSet < tNumNodeSets; ++iSet )
                {
                    // Declare part and add it to the IOBroker (needed for output).
                    stk::mesh::Part& aSetPart = mMtkMeshMetaData->declare_part( mSetNames[0][iSet], stk::topology::NODE_RANK );
                    // Add Part to the IOBroker (needed for output).
                    stk::io::put_io_part_attribute( aSetPart );
                }
            }
        }

    }

    // ----------------------------------------------------------------------------

    // Second declaration to structure the database before filling the data

    void
    Mesh_STK::declare_mesh_fields(
            MtkMeshData &  aMeshData )
    {
        // Fields are data associated with mesh entities. Examples include coordinates, velocity,
        // displacement, and temperature. A field in STK Mesh can hold any data type (e.g., double or int)
        // and any number of scalars per entity (e.g., nodal velocity field has three doubles per node).
        // A field can be allocated (defined) on a whole mesh or on only a subset (part) of that mesh.
        // For example, a material property can be allocated on a specified element block.

        // Declare coordinates field
        Field3CompReal* tCoord_field = &mMtkMeshMetaData->declare_field<Field3CompReal>( stk::topology::NODE_RANK, "coordinates" );
        stk::mesh::put_field( *tCoord_field, mMtkMeshMetaData->universal_part() );


        // Declare all additional fields provided by the user
        if ( aMeshData.FieldsInfo != NULL)
        {
            // Iterate over real scalar fields amd declare them
            uint tNumRealScalarFields = aMeshData.FieldsInfo->get_num_real_scalar_fields();
            for(uint iF = 0; iF<tNumRealScalarFields; iF++)
            {
                Scalar_Field_Info<DDRMat>* tRealScalarField     = (aMeshData.FieldsInfo->mRealScalarFields)(iF);
                enum EntityRank            tFieldEntityRank     = tRealScalarField->get_field_entity_rank();
                std::string                tFieldName           = tRealScalarField->get_field_name();
                stk::mesh::Selector        tFieldPart;
                if(!tRealScalarField->field_has_part_name())
                {
                    tFieldPart = mMtkMeshMetaData->universal_part();
                }
                else
                {
                    tFieldPart = *mMtkMeshMetaData->get_part( tRealScalarField->get_part_name() );
                }

                stk::mesh::Field<real> &   tSTKRealScalarField  =
                        mMtkMeshMetaData->declare_field<stk::mesh::Field<real> >(this->get_stk_entity_rank(tFieldEntityRank),tFieldName , 1);

                stk::mesh::put_field( tSTKRealScalarField, tFieldPart );

                stk::io::set_field_role(tSTKRealScalarField, Ioss::Field::TRANSIENT);
            }

            // iterate over real matrix fields and declare them

            // Iterate over real scalar fields amd declare them
            uint tNumRealMatrixFields = aMeshData.FieldsInfo->get_num_real_matrix_fields();
            for(uint iF = 0; iF<tNumRealMatrixFields; iF++)
            {
                Matrix_Field_Info<DDRMat>* tRealMatrixField = (aMeshData.FieldsInfo->mRealMatrixFields)(iF);
                enum EntityRank            tFieldEntityRank = tRealMatrixField->get_field_entity_rank();
                std::string                tFieldName       = tRealMatrixField->get_field_name();
                const uint                 tNumRows         = tRealMatrixField->get_num_rows();
                const uint                 tNumCols         = tRealMatrixField->get_num_cols();

                stk::mesh::Selector        tFieldPart;
                if(!tRealMatrixField->field_has_part_name())
                {
                    tFieldPart = mMtkMeshMetaData->universal_part();
                }
                else
                {
                    tFieldPart = *mMtkMeshMetaData->get_part( tRealMatrixField->get_part_name() );
                }

                internal_declare_mesh_real_matrix_fields(tFieldName,tFieldEntityRank,tFieldPart,tNumRows,tNumCols);
            }

        }
    }


    void
    Mesh_STK::internal_declare_mesh_real_matrix_fields(
                std::string         aFieldName,
                enum EntityRank     aFieldRank,
                stk::mesh::Selector aFieldPart,
                uint                aNumRows,
                uint                aNumCols)

    {

        uint tNumFieldComp = aNumCols*aNumRows;
        switch ( tNumFieldComp )
        {
            case 1: // Scalar Field
            {
                // Declare fields
                mField1CompVecsReal.push_back( & mMtkMeshMetaData->declare_field< Field1CompReal >( this->get_stk_entity_rank(aFieldRank), aFieldName ) );
                stk::mesh::put_field( *mField1CompVecsReal.back(), aFieldPart, 1 );
                stk::io::set_field_role(*mField1CompVecsReal.back(), Ioss::Field::TRANSIENT);
                break;
            }
            case 2: // Vector Field with 2 components
            {
                // Declare fields
                mField2CompVecsReal.push_back( & mMtkMeshMetaData->declare_field< Field2CompReal >( this->get_stk_entity_rank(aFieldRank), aFieldName ) );
                stk::mesh::put_field( *mField2CompVecsReal.back(), aFieldPart );
                stk::io::set_field_role(*mField2CompVecsReal.back(), Ioss::Field::TRANSIENT);
                break;
            }
            case 3: // Vector Field with 3 components
            {
                // Declare fields
                mField3CompVecsReal.push_back( & mMtkMeshMetaData->declare_field< Field3CompReal >( this->get_stk_entity_rank(aFieldRank), aFieldName ) );
                stk::mesh::put_field( *mField3CompVecsReal.back(), aFieldPart );
                stk::io::set_field_role(*mField3CompVecsReal.back(), Ioss::Field::TRANSIENT);
                break;
            }
            case 4: // Tensor Field with 4 components
            {
                // Declare fields
                mField4CompVecsReal.push_back( & mMtkMeshMetaData->declare_field< Field4CompReal >( this->get_stk_entity_rank(aFieldRank), aFieldName ) );
                stk::mesh::put_field( *mField4CompVecsReal.back(), aFieldPart );
                stk::io::set_field_role(*mField4CompVecsReal.back(), Ioss::Field::TRANSIENT);
                break;
            }
            case 9: // Tensor Field with 9 components
            {
                // Declare fields
                mField9CompVecsReal.push_back( & mMtkMeshMetaData->declare_field< Field9CompReal >( this->get_stk_entity_rank(aFieldRank), aFieldName ) );
                stk::mesh::put_field( *mField9CompVecsReal.back(), aFieldPart );
                stk::io::set_field_role(*mField9CompVecsReal.back(), Ioss::Field::TRANSIENT);
                break;
            }
            default:
            {
                MORIS_ASSERT( 0, "Number of components (columns) for all FieldsData should match one of the supported sizes {1, 2, 3, 4, 9}." );
                break;
            }
        }
    }

    // Declare size of a field (per entity) and throw an error if it is not supported
    void
    Mesh_STK::internal_declare_mesh_field(
            MtkMeshData &  aMeshData,
            uint          iField )
    {
//        // Get field variables
//        uint tNumFieldComp     = aMeshData.FieldsInfo[0].FieldsData[0]( iField ).n_cols();
//        std::string tFieldName = aMeshData.FieldsInfo[0].FieldsName( iField );
//        EntityRank tFieldRank  = aMeshData.FieldsInfo[0].FieldsRank( iField );
//
//        stk::mesh::EntityRank tStkFieldRank = this->get_stk_entity_rank( tFieldRank );
//        stk::mesh::Selector aFieldPart      = mMtkMeshMetaData->universal_part();
//
//        if ( aMeshData.FieldsInfo[0].SetsOwner != NULL )
//        {
//            if ( !aMeshData.FieldsInfo[0].SetsOwner[0]( iField ).empty() )
//            {
//                aFieldPart = *mMtkMeshMetaData->get_part( aMeshData.FieldsInfo[0].SetsOwner[0]( iField ) );
//            }
//        }
//
//        switch ( tNumFieldComp )
//        {
//        case 1: // Scalar Field
//        {
//            // Declare fields
//            mField1CompVec.push_back( & mMtkMeshMetaData->declare_field< Field1Comp >( tStkFieldRank, tFieldName ) );
//            stk::mesh::put_field( *mField1CompVec.back(), aFieldPart, 1 );
//            break;
//        }
//        case 2: // Vector Field with 2 components
//        {
//            // Declare fields
//            mField2CompVec.push_back( & mMtkMeshMetaData->declare_field< Field2Comp >( tStkFieldRank, tFieldName ) );
//            stk::mesh::put_field( *mField2CompVec.back(), aFieldPart );
//            break;
//        }
//        case 3: // Vector Field with 3 components
//        {
//            // Declare fields
//            mField3CompVec.push_back( & mMtkMeshMetaData->declare_field< Field3Comp >( tStkFieldRank, tFieldName ) );
//            stk::mesh::put_field( *mField3CompVec.back(), aFieldPart );
//            break;
//        }
//        case 4: // Tensor Field with 4 components
//        {
//            // Declare fields
//            mField4CompVec.push_back( & mMtkMeshMetaData->declare_field< Field4Comp >( tStkFieldRank, tFieldName ) );
//            stk::mesh::put_field( *mField4CompVec.back(), aFieldPart );
//            break;
//        }
//        case 9: // Tensor Field with 9 components
//        {
//            // Declare fields
//            mField9CompVec.push_back( & mMtkMeshMetaData->declare_field< Field9Comp >( tStkFieldRank, tFieldName ) );
//            stk::mesh::put_field( *mField9CompVec.back(), aFieldPart );
//            break;
//        }
//        default:
//        {
//            MORIS_ASSERT( 0, "Number of components (columns) for all FieldsData should match one of the supported sizes {1, 2, 3, 4, 9}." );
//            break;
//        }
//        }
    }
    // ----------------------------------------------------------------------------

    // Add mesh information to database
    void
    Mesh_STK::populate_mesh_database(
            MtkMeshData &  aMeshData )
    {
        ///////////////////////////////
        // Begin modification cycle  //
        ///////////////////////////////
        mMtkMeshBulkData->modification_begin();

        // Setup global to local element map
        setup_cell_global_to_local_map(aMeshData);

        // Setup global to local node map
        setup_vertex_global_to_local_map(aMeshData);

        // Generate basic mesh information
        this->process_block_sets( aMeshData );

        // Declare node sets to mesh if they exist
        if ( mSetRankFlags[0] )
        {
            this->process_node_sets( aMeshData );
        }

        // Add node sharing
        this->process_node_sharing_information(aMeshData);

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

//------------------------------------------------------------------------------

    void
    Mesh_STK::setup_cell_global_to_local_map(
                    MtkMeshData &   aMeshData )
    {
        uint tNumElemTypes = aMeshData.LocaltoGlobalElemMap.size();

        // Loop over element types
        moris_index tCount = 0;
        for(uint iET = 0; iET<tNumElemTypes; iET++)
        {
            uint tNumElements = aMeshData.LocaltoGlobalElemMap(iET)->numel();

            for(uint iElem = 0; iElem<tNumElements; iElem++)
            {
                uint tElemId = (*aMeshData.LocaltoGlobalElemMap(iET))(iElem);
                if(mEntityGlobaltoLocalMap(3).find(tElemId) == mEntityGlobaltoLocalMap(3).end())
                {
                    mEntityGlobaltoLocalMap(3)[tElemId] = tCount;
                    tCount++;
                }
                else
                {
                    MORIS_ERROR(0,"Element Id already in map, does your LocaltoGlobalElemMap have the same id twice?");
                }
            }
        }
    }

// ----------------------------------------------------------------------------

    void
    Mesh_STK::setup_vertex_global_to_local_map(
                    MtkMeshData &   aMeshData )
    {
        uint tNumNodes = aMeshData.LocaltoGlobalNodeMap->numel();

        // Loop over nodes
        moris_index tCount = 0;
        for(uint iNode = 0; iNode<tNumNodes; iNode++)
        {
            uint tNodeId = (*aMeshData.LocaltoGlobalNodeMap)(iNode);

            if(mEntityGlobaltoLocalMap(0).find(tNodeId) == mEntityGlobaltoLocalMap(0).end())
            {
                mEntityGlobaltoLocalMap(0)[tNodeId] = tCount;
                tCount++;
            }
            else
            {
                MORIS_ERROR(0,"Node Id already in map, does your LocaltoGlobalNodeMap have the same id twice?");
            }
        }

    }


    void
    Mesh_STK::setup_entity_global_to_local_map(enum EntityRank aEntityRank)
    {


        moris::Matrix<IdMat> tEntityIds = get_entities_universal_glob_id(aEntityRank);

        uint tNumEntities = tEntityIds.numel();
        moris_id tCount = 0;
        for(uint i = 0; i<tNumEntities; i++)
        {

            if(mEntityGlobaltoLocalMap((uint)aEntityRank).find(tEntityIds(i)) == mEntityGlobaltoLocalMap((uint)aEntityRank).end())
            {
                mEntityGlobaltoLocalMap((uint)aEntityRank)[tEntityIds(i)] = tCount;
                tCount++;
            }
            else
            {
                MORIS_ERROR(0,"Node Id already in map, does your LocaltoGlobalNodeMap have the same id twice?");
            }

        }


    }

// ----------------------------------------------------------------------------
    // Add all blocks information to database
    void
    Mesh_STK::process_block_sets(
            MtkMeshData &  aMeshData )
    {

        // Get all blocks an element belongs to
        Matrix< IndexMat > tElementToBlock = process_cell_block_membership(aMeshData);

        uint tNumBlocks = 0;
        if(aMeshData.has_mesh_sets())
        {
            tNumBlocks = aMeshData.SetsInfo->get_num_block_sets();
        }

        // Get part vector ordered by the order found in aMeshData block sets
        stk::mesh::PartVector tBlockParts = get_block_set_part_vector(aMeshData);

        // iterate over element types
        for(uint iET = 0; iET<aMeshData.ElemConn.size(); iET++)
        {
            // iterate over elements of this type
            for(uint iElem = 0; iElem<aMeshData.LocaltoGlobalElemMap(iET)->numel(); iElem++)
            {
                moris_id    tCellId    = (*aMeshData.LocaltoGlobalElemMap(iET))(iElem);
                moris_index tCellIndex = get_loc_entity_ind_from_entity_glb_id(tCellId,EntityRank::ELEMENT);

                // get the part vector associated with this cell
                stk::mesh::PartVector tCellParts;
                for(uint iBl = 0; iBl <tNumBlocks; iBl++)
                {
                    if(tElementToBlock(tCellIndex,iBl) != std::numeric_limits<moris_index>::max())
                    {
                        uint tBlockIndex = tElementToBlock(tCellIndex,iBl);
                        stk::mesh::Part* tPart = tBlockParts[tBlockIndex];
                        tCellParts.push_back(tPart);
                    }
                    else
                    {
                        break;
                    }
                }

                if( tCellParts.size() == 0 )
                {
                    stk::mesh::Part* tBlock = mMtkMeshMetaData->get_part("noblock_"+std::to_string(iET) );
                    tCellParts.push_back(tBlock);
                }

                // Declare element
                stk::mesh::Entity tElement = mMtkMeshBulkData->declare_entity(stk::topology::ELEM_RANK,tCellId,tCellParts);

                for (uint node_i = 0; node_i < aMeshData.ElemConn(iET)->n_cols(); ++node_i )
                {
                    moris_id tNodeGlobalId = (stk::mesh::EntityId)(*aMeshData.ElemConn(iET))(iElem,node_i);

                    if(tNodeGlobalId != 0)
                    {
                        stk::mesh::Entity tNode = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, tNodeGlobalId );

                        if(!mMtkMeshBulkData->is_valid( tNode ))
                        {
                            tNode = mMtkMeshBulkData->declare_entity(stk::topology::NODE_RANK, tNodeGlobalId,mMtkMeshMetaData->universal_part());
                        }

                        mMtkMeshBulkData->declare_relation(tElement, tNode, node_i);
                    }
                }
            }
        }
    }



//            // Get all sets provided by the user and go to the block set
//            uint tNumElems     = aMeshData.ElemConn(iET)->numel();
//            uint tNumBlockSets = 1;
//
//            Matrix< IndexMat >  aOwnerPartInds( tNumElems, 1, 0 );
//            std::vector< stk::mesh::PartVector > aPartBlocks( 1 );
//
////            // Update to number of blocks provided by the user
////            if ( mSetRankFlags[2] )
////            {
////                tNumBlockSets  = aMeshData.SetsInfo->get_num_block_sets();
////                aOwnerPartInds = aMeshData.SetsInfo[0].BlockSetsInfo[0].BSetInds[0];
////                aPartBlocks.resize( tNumBlockSets );
////            }
////
////
////            // Populate part blocks
////            for ( uint iSet = 0; iSet < tNumBlockSets; ++iSet )
////            {
////                // Get block sets provided by user
////                stk::mesh::Part* tBlock = mMtkMeshMetaData->get_part( mSetNames[2][iSet] + std::to_string(iET) );
////                aPartBlocks[iSet]       = { tBlock };
////            }
////
////            // Declare MPI communicator
////            stk::ParallelMachine tPM = MPI_COMM_WORLD;
////            uint tParallelSize       = stk::parallel_machine_size( tPM );
////            if ( tParallelSize == 1 )
////            {
////                // serial run
////                this->populate_mesh_database_serial( iET, aMeshData, aPartBlocks, aOwnerPartInds );
////            }
////            else
////            {
////                // Populating mesh is a bit more complicated in parallel because of entity sharing stuff
////                this->populate_mesh_database_parallel( aMeshData, aPartBlocks, aOwnerPartInds );
////            }
////        }

    // ----------------------------------------------------------------------------


    Matrix< IndexMat >
    Mesh_STK::process_cell_block_membership(
                                             MtkMeshData  & aMeshData)
    {

        uint tNumBlocks = 0;
        Matrix< IndexMat > tElementToBlock;
        if(aMeshData.has_mesh_sets())
        {
            tNumBlocks      = aMeshData.SetsInfo->get_num_block_sets();
            tElementToBlock = Matrix< IndexMat >(aMeshData.get_num_elements(),
                                                 tNumBlocks,
                                                 std::numeric_limits<moris_index>::max());
        }

        else
        {
            tElementToBlock = Matrix< IndexMat >(aMeshData.get_num_elements(),
                                                 1,
                                                 std::numeric_limits<moris_index>::max());
        }


        // counter
        Matrix< IndexMat > tElementToBlockCounter(aMeshData.get_num_elements(), 1, 0);

        // Iterate through blocks
        for(uint iBlock = 0; iBlock<tNumBlocks; iBlock++)
        {
            // get pointer to block set
            MtkBlockSetInfo* tBlockSet = aMeshData.SetsInfo->get_block_set(iBlock);

            // Iterate through elements in block
            uint tNumElemsInBlock = tBlockSet->mCellIdsInSet->numel();
            for(uint iElem = 0; iElem<tNumElemsInBlock; iElem++)
            {
                // Cell index from cell id
                moris_index tElemIndex = get_loc_entity_ind_from_entity_glb_id((*tBlockSet->mCellIdsInSet)(iElem), EntityRank::ELEMENT);

                // number of blocks this element belongs to
                uint tElemCount = tElementToBlockCounter(tElemIndex);

                // Add block index to element to block matrix
                tElementToBlock(tElemIndex,tElemCount) = iBlock;

                tElementToBlockCounter(tElemIndex)++;
            }
        }

        return tElementToBlock;
    }

// ----------------------------------------------------------------------------

    stk::mesh::PartVector
    Mesh_STK::get_block_set_part_vector(MtkMeshData &  aMeshData)
    {
        if(aMeshData.has_mesh_sets())
        {
            uint tNumBlocks = aMeshData.SetsInfo->get_num_block_sets();
            stk::mesh::PartVector tBlockSetParts(tNumBlocks);

            // Iterate over block sets and get block part from the block set name
            for(uint iBlock = 0; iBlock<tNumBlocks; iBlock++)
            {
                std::string tBlockName = aMeshData.SetsInfo->get_block_set(iBlock)->mBlockSetName;
                stk::mesh::Part* tBlock = mMtkMeshMetaData->get_part( tBlockName );
                tBlockSetParts[iBlock] = { tBlock };
            }


            return tBlockSetParts;
        }

        else
        {
            return stk::mesh::PartVector(0);
        }
    }

    void
    Mesh_STK::process_node_sharing_information(MtkMeshData& aMeshData)
    {
        if(aMeshData.has_node_sharing_info())
        {
            moris_id tParRank = par_rank();
            uint tNumNodes = aMeshData.get_num_nodes();
            uint tMaxNumProcsShared = aMeshData.NodeProcsShared->n_cols();

            for(uint iNode = 0; iNode<tNumNodes; iNode++ )
            {
                stk::mesh::Entity aNode = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, (*aMeshData.LocaltoGlobalNodeMap)(iNode) );

                for(uint iShare = 0; iShare<tMaxNumProcsShared; iShare++)
                {
                    moris_id tShareProcRank = (*aMeshData.NodeProcsShared)( iNode,iShare );
                    MORIS_ERROR(tParRank != tShareProcRank,"Cannot share a node with self. This causes an issue in STK");
                    if(tShareProcRank != MORIS_ID_MAX )
                    {
                        mMtkMeshBulkData->add_node_sharing( aNode, tShareProcRank );
                    }

                }
            }
        }
    }

// ----------------------------------------------------------------------------
    // Add all fields information to database
    void
    Mesh_STK::populate_mesh_fields(
            MtkMeshData &  aMeshData )
    {
        std::cout<<"populate_mesh_fields start"<<std::endl;
        // Get the coordinates field from Stk
        stk::mesh::FieldBase const* aCoord_field_i = mMtkMeshMetaData->coordinate_field();
        uint tNumNodes                             = aMeshData.NodeCoords->n_rows();

        // Loop over the number of nodes
        for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
        {
            // Get global Id of current node and create "node entity" for stk mesh
            uint aId                = (*aMeshData.LocaltoGlobalNodeMap)( iNode );
            stk::mesh::Entity aNode = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aId );

            // Store the coordinates of the current node
            if ( mMtkMeshBulkData->is_valid( aNode ) )
            {
                // Add coordinates information to the BulkData
                double* tCoord_data = static_cast <double*> ( stk::mesh::field_data ( *aCoord_field_i, aNode ) );
                for ( uint iDim = 0; iDim < mNumDims; ++iDim )
                {
                    tCoord_data[iDim] = (*aMeshData.NodeCoords)( iNode, iDim );
                }
            }
        }

        if ( aMeshData.FieldsInfo!=NULL )
        {
            // Get the number of fields
            uint tNumRealScalarFields = aMeshData.FieldsInfo->get_num_real_scalar_fields();
            for(uint iF = 0; iF<tNumRealScalarFields; iF++)
            {
                Scalar_Field_Info<DDRMat>* tRealScalarField = (aMeshData.FieldsInfo->mRealScalarFields)(iF);

                if(tRealScalarField->field_has_data())
                {
                    populate_field_data_scalar_field(tRealScalarField);
                }
            }

            // Get the number of fields : fixme: causes segfault
            /*uint tNumSintScalarFields = aMeshData.FieldsInfo->get_num_sint_scalar_fields();
            for(uint iF = 0; iF<tNumSintScalarFields; iF++)
            {
                Scalar_Field_Info<DDSMat>* tSintScalarField = (aMeshData.FieldsInfo->mSintScalarFields)(iF);
                if(tSintScalarField->field_has_data())
                {
                    populate_field_data_scalar_field( tSintScalarField );
                }
            }*/

            // Get the number of fields
            uint tNumRealMatrixFields = aMeshData.FieldsInfo->get_num_real_matrix_fields();
            for(uint iF = 0; iF<tNumRealMatrixFields; iF++)
            {
                Matrix_Field_Info<DDRMat>* tRealMatrixField = (aMeshData.FieldsInfo->mRealMatrixFields)(iF);
                if(tRealMatrixField->field_has_data())
                {
                    populate_field_data_matrix_field(tRealMatrixField);
                }
            }
        }
        std::cout<<"populate_mesh_fields end"<<std::endl;
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
                std::cout<<"Quad 4s"<<std::endl;
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
                std::cout<<"Quad 9s"<<std::endl;
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
                std::cout<<"Hex 8s"<<std::endl;
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
                std::cout<<"Hex 27s"<<std::endl;
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


//    stk::topology::topology_t
//    Mesh_STK::get_stk_entity_rank(enum EntityRank aMTKEntityRank)
//    {
//
//        stk::topology::topology_t tTopology = stk::topology::INVALID_TOPOLOGY;
//
//        switch ( aMTKEntityRank )
//        {
//            case EntityRank::NODE:
//            {
//                tTopology = stk::topology::NODE_RANK;
//                break;
//            }
//            case EntityRank::EDGE:
//            {
//                tTopology = stk::topology::EDGE_RANK;
//                break;
//            }
//            case EntityRank::FACE:
//            {
//                tTopology = stk::topology::FACE_RANK;
//                break;
//            }
//            case EntityRank::ELEMENT:
//            {
//                tTopology = stk::topology::ELEMENT_RANK;
//                break;
//            }
//            default:
//            {
//                MORIS_ASSERT( 0, "Invalid entity rank.");
//                break;
//            }
//        }
//    }


    // Provide element type (Hex8, Tri3, etc) and throw error if element is not supported yet.
    stk::topology::topology_t
    Mesh_STK::get_stk_topo(enum CellTopology aMTKCellTopo )
    {
        stk::topology::topology_t tTopology = stk::topology::INVALID_TOPOLOGY;

        switch ( aMTKCellTopo )
        {
            case CellTopology::TET4:
            {
                tTopology = stk::topology::TET_4;
                break;
            }
            case CellTopology::HEX8:
            {
                tTopology = stk::topology::HEX_8;
                break;
            }
            default:
            {
                MORIS_ASSERT( 0, "MTK mesh build from data currently handles only TET_4, HEX8,  and  for 3D elements.");
                break;
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
            Matrix< IdMat >  tEdgeIds = this->get_entities_universal_glob_id( EntityRank::EDGE );

            // resize member variable to its right size
            mEntityLocaltoGlobalMap(1) = Matrix<IdMat>( tNumEdges, 1 );

            // Populate internal member variable that contains the local index to
            // global id node communication information
            for ( uint iEdge = 0; iEdge < tNumEdges; ++iEdge )
            {
                // local to global and owner processor
                mEntityLocaltoGlobalMap(1)( iEdge ) = tEdgeIds( iEdge );
            }
        }

        if ( mNumDims > 2 )
        {
            uint tNumFaces   = this->get_num_faces();

            // Access entities stored in mesh database
            Matrix< IdMat >  tFaceIds = this->get_entities_universal_glob_id( EntityRank::FACE );

            // resize member variable to its right size
            mEntityLocaltoGlobalMap(2)= Matrix<IdMat>( tNumFaces, 1 );

            // Populate internal member variable that contains the local index to
            // global id element communication information
            for ( uint iFace = 0; iFace < tNumFaces; ++iFace )
            {
                // local to global and owner processor
                mEntityLocaltoGlobalMap(2)( iFace ) = tFaceIds( iFace );
            }
        }
    }
// ----------------------------------------------------------------------------
    // Function to create edges and faces communication lists in parallel
    void
    Mesh_STK::create_owners_communication_lists()
    {
    }
// ----------------------------------------------------------------------------
    // Function to create edges and faces communication lists in parallel
    void
    Mesh_STK::create_shared_communication_lists()
    {
        // Get basic mesh information
        Matrix< IdMat >  tNodesShared = this->get_entities_glb_shared_current_proc( EntityRank::NODE );
        uint tNumNodesShared      = tNodesShared.length();

        // Generate list of processors sharing information
        // -----------------------------------------------

        std::vector < uint > tActiveSharedProcs;

        // Loop over the number of nodes shared to get the shared processors
        for( uint iNodeShared = 0; iNodeShared < tNumNodesShared; ++iNodeShared )
        {
            Matrix< IdMat >  tProcsSharing = this->get_procs_sharing_entity_by_id( tNodesShared( iNodeShared ), EntityRank::NODE );

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
        uint tNumNodeSets = aMeshData.SetsInfo->get_num_node_sets();
        stk::mesh::EntityRank aStkSetRank  = stk::topology::NODE_RANK;

        for ( uint iSet = 0; iSet < tNumNodeSets; ++iSet )
        {
            // STK interface variables declaration
            stk::mesh::Part* aSetPart = mMtkMeshMetaData->get_part( mSetNames[0][iSet] );
            stk::mesh::PartVector aAddPart( 1, aSetPart );
            stk::mesh::EntityId aGlobalId;
            stk::mesh::Entity aEntity;

            // Get node set and size of node set
            MtkNodeSetInfo* tNodeSet = aMeshData.SetsInfo->get_node_set(iSet);
            uint tNumNodesInSet = tNodeSet->mNodeIds->numel();

            // Populate node sets (change entity parts if nodes were declared already)
            for ( uint iEntity = 0; iEntity < tNumNodesInSet; ++iEntity )
            {
                // Declare new entity or add existing entity to declared part
                aGlobalId = (*tNodeSet->mNodeIds)( iEntity );
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
        uint tNumSideSets = aMeshData.SetsInfo->get_num_side_sets();

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

            // Get side set and size of node set
            MtkSideSetInfo* tSideSet = aMeshData.SetsInfo->get_side_set(iSet);
            uint tNumSidesInSet = tSideSet->mElemIdsAndSideOrds->n_rows();

            for ( uint iEntity = 0; iEntity < tNumSidesInSet; ++iEntity )
            {
                // First column contains element ids that will later be match with faces
                aGlobalElemId     = (*tSideSet->mElemIdsAndSideOrds)( iEntity, 0 );
                aElemEntity       = mMtkMeshBulkData->get_entity( stk::topology::ELEMENT_RANK, aGlobalElemId );
                tRequestedSideOrd = (*tSideSet->mElemIdsAndSideOrds)( iEntity, 1 );

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
            moris::uint  aElementTypeInd,
            MtkMeshData                            aMeshData,
            std::vector< stk::mesh::PartVector >   aElemParts,
            Matrix< IdMat >                       aOwnerPartInds )
    {

            uint tNumElems        = aMeshData.ElemConn(aElementTypeInd)->n_rows();
            uint aNumNodesPerElem = aMeshData.ElemConn(aElementTypeInd)->n_cols();

            // Declare variables to access connectivity
            Matrix< IdMat >  tDummyMat( 1, aNumNodesPerElem );
            stk::mesh::EntityIdVector aCurrElemConn( aNumNodesPerElem );
            stk::mesh::EntityId aElemGlobalId;

            // Loop over the number of elements and interface between MORIS and Stk for connectivity
            for ( uint iElem = 0; iElem < tNumElems; ++iElem )
            {
                // Get row of nodes connected in moris variable and assign to STK variable
                tDummyMat.set_row( 0,aMeshData.ElemConn(aElementTypeInd)->get_row( iElem ));
                aCurrElemConn.assign( tDummyMat.data(), tDummyMat.data() + aNumNodesPerElem );

                // Declare element in function that also declares element-node relations internally
                aElemGlobalId = (*aMeshData.LocaltoGlobalElemMap(aElementTypeInd))( iElem );
                stk::mesh::declare_element( *mMtkMeshBulkData, aElemParts[aOwnerPartInds( iElem )], aElemGlobalId, aCurrElemConn );
            }
    }
// ----------------------------------------------------------------------------
    // Parallel specific implementation for blocks in database
    void
    Mesh_STK::populate_mesh_database_parallel(
            MtkMeshData                          aMeshData,
            std::vector< stk::mesh::PartVector > aPartBlocks,
            Matrix< IdMat >                     aOwnerPartInds )
    {
        for(uint iET = 0; iET<aMeshData.ElemConn.size(); iET++)
        {
            // Mesh variables
            uint tNumNodes        = aMeshData.NodeCoords[0].n_rows(  );
            uint tNumElems        = aMeshData.ElemConn(iET)->n_rows(  );
            uint aNumNodesPerElem = aMeshData.ElemConn(iET)->n_cols(  );

            // Declare MPI communicator
            stk::ParallelMachine tPM = MPI_COMM_WORLD;
            uint tProcRank     = stk::parallel_machine_rank( tPM );

            // Check if the list provided is for nodes or elements shared.
            bool tNodesProcOwnerList = false;

            if ( aMeshData.NodeProcOwner[0].n_rows() == tNumNodes )
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
                    tDummyMat.get_row( 0 ) = aMeshData.ElemConn(iET)->get_row( iElem );
                    aCurrElemConn.assign( tDummyMat.data(), tDummyMat.data() + aNumNodesPerElem );
                    aElemGlobalId = (*aMeshData.LocaltoGlobalElemMap(iET))( iElem );

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
                    if ( aMeshData.NodeProcOwner[0]( iNode ) != (moris_id)tProcRank )
                    {
                        // Add node if it has not been assign to the mesh yet.
                        stk::mesh::Entity aNode = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aMeshData.LocaltoGlobalNodeMap[0]( iNode, 0 ) );
                        mMtkMeshBulkData->add_node_sharing( aNode, aMeshData.NodeProcOwner[0]( iNode ) );
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
                    if ( aMeshData.NodeProcOwner[0]( iElem ) == (moris_id)tProcRank )
                    {
                        // Get row of nodes connected in moris variable and assign to STK variable
                        tDummyMat.get_row( 0 ) = aMeshData.ElemConn(iET)->get_row( iElem );
                        aCurrElemConn.assign( tDummyMat.data(), tDummyMat.data() + aNumNodesPerElem );

                        // Loop over the number of nodes in each element
                        for ( uint iNode = 0; iNode < aNumNodesPerElem; ++iNode )
                        {
                            // Get global Id of current node and create "node entity" for stk mesh
                            aNodeGlobalId = (*aMeshData.ElemConn(iET))( iElem, iNode );
                            aNode = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aNodeGlobalId );

                            // Add node if it has not been assign to the mesh yet.
                            if ( !mMtkMeshBulkData->is_valid( aNode ) )
                            {
                                tNodesIDeclared.insert( aNodeGlobalId );
                            }
                        }

                        // declare element in function that also declares element-node relations internally
                        aElemGlobalId = (*aMeshData.LocaltoGlobalElemMap(iET))( iElem );
                        stk::mesh::declare_element( *mMtkMeshBulkData, aPartBlocks[aOwnerPartInds( iElem )], aElemGlobalId, aCurrElemConn );
                    }
                }

                // NOTE: This implementation requires the incorporation of the elements located at the boundaries of the processors
                // (i.e., aura elements) to the connectivity table and the elements processor owner list.

                // Loop over the number of elements
                for ( uint iElem = 0; iElem < tNumElems; ++iElem )
                {
                    // Check if the element is not own by this processor
                    if ( aMeshData.NodeProcOwner[0]( iElem )!= (moris_id)tProcRank )
                    {
                        // Loop over the number of nodes in each element
                        for ( uint iNode = 0; iNode < aNumNodesPerElem; ++iNode )
                        {
                            aNodeGlobalId = (*aMeshData.ElemConn(iET))( iElem, iNode );

                            if ( tNodesIDeclared.find( aNodeGlobalId ) != tNodesIDeclared.end() )
                            {
                                // Add node if it has not been assign to the mesh yet.
                                aNode = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aNodeGlobalId );
                                mMtkMeshBulkData->add_node_sharing( aNode, aMeshData.NodeProcOwner[0]( iElem ) );
                            }
                        }
                    }
                }
            }
        }
    }
// ----------------------------------------------------------------------------
    // Access set entity ids
    moris::Matrix< IdMat >
    Mesh_STK::get_set_entity_glob_ids(
            stk::mesh::EntityRank aEntityRank,
            std::string           aSetName ) const
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


