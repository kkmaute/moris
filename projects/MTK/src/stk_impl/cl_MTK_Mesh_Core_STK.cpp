/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Core_STK.cpp
 *
 */

#include "Ioss_Region.h"      // for Region, NodeBlockContainer
#include <Teuchos_RCP.hpp>    // for RCP::RCP<T>, RCP::operator*, etc

#include <stk_mesh/base/GetEntities.hpp>          // for count_entities
#include <stk_mesh/base/Selector.hpp>             // for Selector
#include <stk_mesh/base/FEMHelpers.hpp>           // for Selector
#include "stk_io/DatabasePurpose.hpp"             // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/CoordinateSystems.hpp"    // for Cartesian
#include "stk_mesh/base/CreateFaces.hpp"          // for handling faces
#include "stk_mesh/base/CreateEdges.hpp"          // for handling faces
#include "stk_mesh/base/Bucket.hpp"               // for buckets
#include "stk_mesh/base/Field.hpp"                // for coordinates
#include "stk_mesh/base/GetEntities.hpp"          // for coordinates
#include "stk_mesh/base/FieldParallel.hpp"        // for handling parallel fields
#include <exodusII.h>

#include "fn_assert.hpp"
#include "fn_isempty.hpp"
#include "fn_find.hpp"
#include "fn_sort.hpp"
#include "fn_unique.hpp"
#include "op_equal_equal.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_Tracer.hpp"

namespace moris
{
    namespace mtk
    {
        //##############################################
        // Build mesh functionality
        //##############################################

        // ----------------------------------------------------------------------------

        Mesh_Core_STK::Mesh_Core_STK(
                std::string  aFileName,
                MtkMeshData* aSupMeshData,
                const bool   aCreateFacesAndEdges )
        {
            // Call the function that handles the communication between stk and moris
            this->build_mesh( aFileName, aSupMeshData, aCreateFacesAndEdges );
        }

        // ----------------------------------------------------------------------------

        Mesh_Core_STK::Mesh_Core_STK(
                std::shared_ptr< Mesh_Data_STK > aSTKMeshData )
                : mSTKMeshData( aSTKMeshData )
        {
        }

        // ----------------------------------------------------------------------------

        Mesh_Core_STK::~Mesh_Core_STK()
        {
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::build_mesh(
                std::string  aFileName,
                MtkMeshData* aSuppMeshData,
                const bool   aCreateFacesAndEdges )
        {
            Tracer tTracer( "STK", "Build" );
            // make a shared pointer of the stk mesh data
            mSTKMeshData = std::make_shared< Mesh_Data_STK >();

            // Start the timer
            std::clock_t start = std::clock();

            // The 'generated:' syntax in fileName makes a hex mesh to be generated in memory.
            // If an Exodus file name is used instead, then the mesh is read from file. The code
            // below remains the same in either case.

            // Declare MPI communicator
            MPI_Comm aCommunicator = MPI_COMM_WORLD;

            // Generate MetaData and Bulk Data instances (later to be pointed to member variables)
            stk::mesh::MetaData* meshMeta = new stk::mesh::MetaData;

            MORIS_LOG( "Construct Bulk Data" );
            stk::mesh::BulkData* meshBulk = new stk::mesh::BulkData( *meshMeta, aCommunicator, this->get_aura_option() );

            // Set member variables as pointers to meta_data and bulk_data
            mSTKMeshData->mMtkMeshMetaData = ( meshMeta );
            mSTKMeshData->mMtkMeshBulkData = ( meshBulk );

            // Use STK IO to populate a STK Mesh
            MORIS_LOG( "Mesh Reader Initialize" );
            mSTKMeshData->mMeshReader = new stk::io::StkMeshIoBroker( aCommunicator );

            // Create mesh database using the IO broker
            mSTKMeshData->mMeshReader->set_bulk_data( *meshBulk );
            MORIS_LOG( "Read mesh from file" );
            mSTKMeshData->mMeshReader->add_mesh_database( aFileName, stk::io::READ_MESH );
            mSTKMeshData->mMeshReader->create_input_mesh();

            mSTKMeshData->mNumDims = mSTKMeshData->mMtkMeshMetaData->spatial_dimension();

            // Include mesh fields and populate the database
            MORIS_LOG( "Add Mesh Fields" );
            mSTKMeshData->mMeshReader->add_all_mesh_fields_as_input_fields( stk::io::MeshField::CLOSEST );
            // Declare supplementary fields
            if ( aSuppMeshData != nullptr )
            {
                declare_mesh_fields( *aSuppMeshData );
                add_supplementary_fields_to_declare_at_output( *aSuppMeshData );
            }

            // Add fields to fields to declare upon output list

            // Create nodesets and sidesets
            mSTKMeshData->mSetNames.resize( 3 );    // Number of ranks for sets

            MORIS_LOG( "Populate Bulk Data" );
            mSTKMeshData->mMeshReader->populate_bulk_data();

            // Determine number of time increments on input database in region
            Teuchos::RCP< Ioss::Region > tIo_region      = mSTKMeshData->mMeshReader->get_input_io_region();
            int                          tTimestep_count = tIo_region->get_property( "state_count" ).get_int();

            MORIS_LOG( "Read Input Fields" );
            // Loop over all time increments and get field data
            for ( int step = 1; step <= tTimestep_count; step++ )
            {
                mSTKMeshData->mMeshReader->read_defined_input_fields( step );
            }

            if ( aCreateFacesAndEdges )
            {
                MORIS_LOG( "Create Edges and Faces" );
                // Create mesh edge entities
                stk::mesh::create_edges( *mSTKMeshData->mMtkMeshBulkData );
                mSTKMeshData->mCreatedEdges = true;

                // Create mesh face entities
                stk::mesh::create_faces( *mSTKMeshData->mMtkMeshBulkData, true );
                mSTKMeshData->mCreatedFaces = true;
            }

            MORIS_LOG( "Communicate fields on aura" );
            const stk::mesh::FieldVector&              fields = mSTKMeshData->mMtkMeshMetaData->get_fields();
            std::vector< const stk::mesh::FieldBase* > const_fields( fields.size() );
            for ( size_t i = 0; i < fields.size(); ++i )
            {
                const_fields[ i ] = fields[ i ];
            }

            for ( size_t ghost_i = 0; ghost_i < mSTKMeshData->mMtkMeshBulkData->ghostings().size(); ++ghost_i )
            {
                stk::mesh::communicate_field_data( *( mSTKMeshData->mMtkMeshBulkData->ghostings()[ ghost_i ] ), const_fields );
            }

            MORIS_LOG( "Create local to global maps" );

            // Create communication tables in parallel.
            // NOTE1: the information to be created in the function below duplicates communication-related data
            // already created in the STK mesh database. However, one can access that data on a 1 by 1 case scenario
            // which could be inefficient in certain operations. For that reason, communication lists that live in
            // the mesh class as member variables are used instead.
            // NOTE2: this information is already provided by the user in meshes generated from data, which implies
            // that this call is only required here (mesh from string or file).

            mSTKMeshData->mEntityLocaltoGlobalMap =
                    moris::Cell< moris::Matrix< IdMat > >( (uint)EntityRank::END_ENUM, moris::Matrix< IndexMat >( 1, 1, 0 ) );

            mSTKMeshData->mEntitySendList =
                    moris::Cell< moris::Cell< moris::Matrix< IndexMat > > >( (uint)EntityRank::END_ENUM, moris::Cell< moris::Matrix< IndexMat > >( par_size(), moris::Matrix< IdMat >( 1, 1 ) ) );

            mSTKMeshData->mEntityReceiveList =
                    moris::Cell< moris::Cell< moris::Matrix< IndexMat > > >( (uint)EntityRank::END_ENUM, moris::Cell< moris::Matrix< IndexMat > >( par_size(), moris::Matrix< IdMat >( 1, 1 ) ) );

            create_communication_lists_and_local_to_global_map( EntityRank::NODE );
            create_communication_lists_and_local_to_global_map( EntityRank::EDGE );
            create_communication_lists_and_local_to_global_map( EntityRank::FACE );
            create_communication_lists_and_local_to_global_map( EntityRank::ELEMENT );

            MORIS_LOG( "Create global to local maps" );
            // Initialize global to local map
            mSTKMeshData->mEntityGlobaltoLocalMap = moris::Cell< std::unordered_map< moris_id, moris_index > >( 4 );
            setup_entity_global_to_local_map( EntityRank::NODE );
            setup_entity_global_to_local_map( EntityRank::EDGE );
            setup_entity_global_to_local_map( EntityRank::FACE );
            setup_entity_global_to_local_map( EntityRank::ELEMENT );

            MORIS_LOG( "Setup vertex pairing" );
            // setup vertex pairing
            setup_parallel_vertex_pairing();

            MORIS_LOG( "Setup cell sharing" );
            // setup cell sharing
            setup_parallel_cell_sharing();

            MORIS_LOG( "Setup vertices and cells" );
            setup_vertices_and_cell();

            if ( mVerbose )
            {
                MORIS_LOG_INFO( "MTK: Load mesh from file completed in %f s", ( std::clock() - start ) / (double)( CLOCKS_PER_SEC ) );
            }
        }

        // ----------------------------------------------------------------------------

        Mesh_Core_STK::Mesh_Core_STK( MtkMeshData& aMeshData )
        {
            // allocate stk mesh data
            mSTKMeshData                          = std::make_shared< Mesh_Data_STK >();
            mSTKMeshData->mEntityLocaltoGlobalMap = moris::Cell< moris::Matrix< IdMat > >( 4 );
            mSTKMeshData->mEntityGlobaltoLocalMap = moris::Cell< std::unordered_map< moris_id, moris_index > >( 4 );
            mSTKMeshData->mSetRankFlags           = std::vector< bool >( { false, false, false } );

            // Set verbose flag
            mVerbose = aMeshData.Verbose;

            // start timing
            std::clock_t start = std::clock();

            // Flag for handling data generated mesh
            mSTKMeshData->mDataGeneratedMesh = true;

            // Verify the minimum required data was provided and initialize uninitialized variables if necessary
            this->check_and_update_input_data( aMeshData );

            // Call the function that handles the communication between stk and moris
            this->build_mesh( aMeshData );

            if ( mVerbose )
            {
                MORIS_LOG_INFO( "MTK: Create mesh from data completed in %f s", ( std::clock() - start ) / (double)( CLOCKS_PER_SEC ) );
            }
        }

        // ----------------------------------------------------------------------------

        //##############################################
        // General mesh information access
        //##############################################
        uint
        Mesh_Core_STK::get_num_entities(
                enum EntityRank   aEntityRank,
                const moris_index aIndex ) const
        {
            // Initialize
            stk::mesh::EntityRank requestedRank = this->get_stk_entity_rank( aEntityRank );

            // Get all entities from meta data
            stk::mesh::Selector allEntities = mSTKMeshData->mMtkMeshMetaData->universal_part();

            uint tNumSTKEntities = (uint)stk::mesh::count_selected_entities( allEntities, mSTKMeshData->mMtkMeshBulkData->buckets( requestedRank ) );

            return tNumSTKEntities;
        }

        // ----------------------------------------------------------------------------

        //##############################################
        // Access Mesh Data by index Functions
        //##############################################
        Matrix< IndexMat >
        Mesh_Core_STK::get_entity_connected_to_entity_loc_inds(
                moris_index       aEntityIndex,
                enum EntityRank   aInputEntityRank,
                enum EntityRank   aOutputEntityRank,
                const moris_index aIndex ) const
        {
            MORIS_ERROR( aInputEntityRank != aOutputEntityRank,
                    " Input and output entity rank cannot be the same (this is an invalid connectivity inside STK). Use get_element_to_element_loc_inds for element to element connectivity." );

            // Get Stk entity Id from local to global map
            moris_id            tId          = (moris_id)this->get_glb_entity_id_from_entity_loc_index( aEntityIndex, aInputEntityRank );
            stk::mesh::EntityId tStkEntityId = (stk::mesh::EntityId)tId;

            // Call function that gets the connected entities
            stk::mesh::EntityRank tStkInputRank  = this->get_stk_entity_rank( aInputEntityRank );
            stk::mesh::EntityRank tStkOutputRank = this->get_stk_entity_rank( aOutputEntityRank );
            stk::mesh::Entity     tStkEntity     = mSTKMeshData->mMtkMeshBulkData->get_entity( tStkInputRank, tStkEntityId );

            std::vector< stk::mesh::Entity > tEntitiesConnected = this->entities_connected_to_entity_stk( &tStkEntity, tStkInputRank, tStkOutputRank );

            // Get the number of entities
            uint tNumOutputEntities = tEntitiesConnected.size();

            // Declare xtk::Mat that will contain the local entities
            Matrix< IndexMat > tLocalIndices( 1, tNumOutputEntities );

            // Fill local ids to xtk::Mat
            for ( uint i = 0; i < tNumOutputEntities; ++i )
            {
                moris_id tId          = mSTKMeshData->mMtkMeshBulkData->identifier( tEntitiesConnected[ i ] );
                tLocalIndices( 0, i ) = get_loc_entity_ind_from_entity_glb_id( tId, aOutputEntityRank );
            }
            return tLocalIndices;
        }

        // ----------------------------------------------------------------------------

        Matrix< IndexMat >
        Mesh_Core_STK::get_elements_connected_to_element_and_face_ord_loc_inds( moris_index aElementIndex ) const
        {
            // First get faces connected to element
            // Get Stk entity Id from local to global map
            moris_id            tId          = this->get_glb_entity_id_from_entity_loc_index( aElementIndex, EntityRank::ELEMENT );
            stk::mesh::EntityId tStkEntityId = (stk::mesh::EntityId)tId;

            enum EntityRank tFacetRank = this->get_facet_rank();

            // Call function that gets the connected entities
            stk::mesh::EntityRank tStkInputRank  = stk::topology::ELEMENT_RANK;
            stk::mesh::EntityRank tStkOutputRank = stk::topology::FACE_RANK;
            if ( tFacetRank == EntityRank::EDGE )
            {
                tStkOutputRank = stk::topology::EDGE_RANK;
            }

            stk::mesh::Entity tStkEntity = mSTKMeshData->mMtkMeshBulkData->get_entity( tStkInputRank, tStkEntityId );

            std::vector< stk::mesh::Entity > tFacesInElem = this->entities_connected_to_entity_stk( &tStkEntity, tStkInputRank, tStkOutputRank );

            MORIS_ASSERT( ( tFacesInElem.size() != 0 ) || ( tFacesInElem.size() != 0 ),
                    "No faces connected to element found. Maybe the CreateAllEdgesAndFaces flag is set to false. Check mesh struct." );

            // Then for each face get elements connected
            uint tCounter  = 0;
            uint tNumFaces = tFacesInElem.size();

            moris::Matrix< IndexMat > tElemsConnectedToElem( 4, tNumFaces );

            for ( uint faceIt = 0; faceIt < tNumFaces; ++faceIt )
            {
                std::vector< stk::mesh::Entity > tDummyConnectivity =
                        this->entities_connected_to_entity_stk( &tFacesInElem[ faceIt ], stk::topology::FACE_RANK, stk::topology::ELEMENT_RANK );

                moris_index tFacetId = (moris_index)mSTKMeshData->mMtkMeshBulkData->local_id( tFacesInElem[ faceIt ] );

                // Faces in mesh boundaries do not have more than one element
                if ( tDummyConnectivity.size() > 0 )
                {
                    if ( mSTKMeshData->mMtkMeshBulkData->identifier( tDummyConnectivity[ 0 ] ) != mSTKMeshData->mMtkMeshBulkData->identifier( tStkEntity ) )
                    {
                        tElemsConnectedToElem( 0, tCounter ) = (moris_index)mSTKMeshData->mMtkMeshBulkData->local_id( tDummyConnectivity[ 0 ] );
                        tElemsConnectedToElem( 1, tCounter ) = faceIt;
                        tElemsConnectedToElem( 2, tCounter ) = this->get_facet_ordinal_from_cell_and_facet_loc_inds( tFacetId, tElemsConnectedToElem( 0, tCounter ) );
                        tElemsConnectedToElem( 3, tCounter ) = MORIS_INDEX_MAX;
                        tCounter++;
                    }
                }

                if ( tDummyConnectivity.size() > 1 )
                {
                    if ( mSTKMeshData->mMtkMeshBulkData->identifier( tDummyConnectivity[ 1 ] ) != mSTKMeshData->mMtkMeshBulkData->identifier( tStkEntity ) )
                    {
                        tElemsConnectedToElem( 0, tCounter ) = (moris_index)mSTKMeshData->mMtkMeshBulkData->local_id( tDummyConnectivity[ 1 ] );
                        tElemsConnectedToElem( 1, tCounter ) = faceIt;
                        tElemsConnectedToElem( 2, tCounter ) = this->get_facet_ordinal_from_cell_and_facet_loc_inds( tFacetId, tElemsConnectedToElem( 0, tCounter ) );
                        tElemsConnectedToElem( 3, tCounter ) = MORIS_INDEX_MAX;

                        tCounter++;
                    }
                }

                MORIS_ASSERT( tDummyConnectivity.size() <= 2,
                        "For some reason face has more than 2 elements connected to it... Check get_elements_connected_to_element." );
            }

            // Resize to include only ids added above and get rid of initialized extra zeros
            tElemsConnectedToElem.resize( 4, tCounter );

            return tElemsConnectedToElem;
        }

        // ----------------------------------------------------------------------------

        Matrix< IndexMat >
        Mesh_Core_STK::get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const
        {
            // First get faces connected to element
            // Get Stk entity Id from local to global map
            moris_id            tId          = this->get_glb_entity_id_from_entity_loc_index( aElementIndex, EntityRank::ELEMENT );
            stk::mesh::EntityId tStkEntityId = (stk::mesh::EntityId)tId;

            // Call function that gets the connected entities
            enum EntityRank tFacetRank = this->get_facet_rank();

            stk::mesh::EntityRank tStkInputRank  = stk::topology::ELEMENT_RANK;
            stk::mesh::EntityRank tStkOutputRank = stk::topology::FACE_RANK;
            if ( tFacetRank == EntityRank::EDGE )
            {
                tStkOutputRank = stk::topology::EDGE_RANK;
            }

            stk::mesh::Entity tStkEntity = mSTKMeshData->mMtkMeshBulkData->get_entity( tStkInputRank, tStkEntityId );

            std::vector< stk::mesh::Entity > tFacetsInElem = this->entities_connected_to_entity_stk( &tStkEntity, tStkInputRank, tStkOutputRank );

            MORIS_ASSERT( ( tFacetsInElem.size() != 0 ) || ( tFacetsInElem.size() != 0 ),
                    "No faces connected to element found. Maybe the CreateAllEdgesAndFaces flag is set to false. Check mesh struct." );

            // Then for each face get elements connected
            uint tCounter  = 0;
            uint tNumFaces = tFacetsInElem.size();

            moris::Matrix< IndexMat > tElemsConnectedToElem( 2, tNumFaces );

            for ( uint faceIt = 0; faceIt < tNumFaces; ++faceIt )
            {
                std::vector< stk::mesh::Entity > tDummyConnectivity =
                        this->entities_connected_to_entity_stk( &tFacetsInElem[ faceIt ], tStkOutputRank, stk::topology::ELEMENT_RANK );

                // Faces in mesh boundaries do not have more than one element
                if ( tDummyConnectivity.size() > 0 )
                {
                    if ( mSTKMeshData->mMtkMeshBulkData->identifier( tDummyConnectivity[ 0 ] ) != mSTKMeshData->mMtkMeshBulkData->identifier( tStkEntity ) )
                    {
                        tElemsConnectedToElem( 0, tCounter ) = (moris_index)mSTKMeshData->mMtkMeshBulkData->local_id( tDummyConnectivity[ 0 ] );
                        tElemsConnectedToElem( 1, tCounter ) = (moris_index)mSTKMeshData->mMtkMeshBulkData->local_id( tFacetsInElem[ faceIt ] );
                        ;
                        tCounter++;
                    }
                }

                if ( tDummyConnectivity.size() > 1 )
                {
                    if ( mSTKMeshData->mMtkMeshBulkData->identifier( tDummyConnectivity[ 1 ] ) != mSTKMeshData->mMtkMeshBulkData->identifier( tStkEntity ) )
                    {
                        tElemsConnectedToElem( 0, tCounter ) = (moris_index)mSTKMeshData->mMtkMeshBulkData->local_id( tDummyConnectivity[ 1 ] );
                        tElemsConnectedToElem( 1, tCounter ) = (moris_index)mSTKMeshData->mMtkMeshBulkData->local_id( tFacetsInElem[ faceIt ] );
                        ;
                        tCounter++;
                    }
                }

                MORIS_ASSERT( tDummyConnectivity.size() <= 2,
                        "For some reason face has more than 2 elements connected to it... Check get_elements_connected_to_element." );
            }

            // Resize to include only ids added above and get rid of initialized extra zeros
            tElemsConnectedToElem.resize( 2, tCounter );

            return tElemsConnectedToElem;
        }

        // ----------------------------------------------------------------------------

        moris::Cell< moris::mtk::Vertex const * >
        Mesh_Core_STK::get_all_vertices() const
        {
            enum EntityRank tEntityRank = EntityRank::NODE;

            // Get pointer to field defined by input name
            // set and set selector
            stk::mesh::Part&    tUniversal = mSTKMeshData->mMtkMeshMetaData->universal_part();
            stk::mesh::Selector tUniversalSelector( tUniversal );

            // aura and aura selector
            //        stk::mesh::Part & tAura = mSTKMeshData->mMtkMeshMetaData->aura_part();
            //        stk::mesh::Selector tAuraSelector( tAura );

            // difference of the set with the auraz
            //        stk::mesh::Selector tDifference = tUniversalSelector -= tAuraSelector;

            stk::mesh::EntityVector aEntities;
            stk::mesh::get_selected_entities( tUniversalSelector, mSTKMeshData->mMtkMeshBulkData->buckets( stk::topology::NODE_RANK ), aEntities );

            // Get entity Ids
            uint tNumEntities = aEntities.size();

            moris::Cell< moris::mtk::Vertex const * > tOutputVertices( tNumEntities );
            for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
            {
                moris::moris_index tVertexIndex =
                        this->get_loc_entity_ind_from_entity_glb_id( (moris_id)mSTKMeshData->mMtkMeshBulkData->identifier( aEntities[ iEntity ] ), tEntityRank );

                tOutputVertices( tVertexIndex ) = &this->get_mtk_vertex( tVertexIndex );
            }

            return tOutputVertices;
        }

        // ----------------------------------------------------------------------------

        //##############################################
        // Local to global functions
        //##############################################

        moris_id
        Mesh_Core_STK::get_glb_entity_id_from_entity_loc_index(
                moris_index       aEntityIndex,
                enum EntityRank   aEntityRank,
                const moris_index aIndex ) const
        {
            return mSTKMeshData->mEntityLocaltoGlobalMap( (uint)aEntityRank )( aEntityIndex );
        }

        // ----------------------------------------------------------------------------

        moris_index
        Mesh_Core_STK::get_loc_entity_ind_from_entity_glb_id(
                moris_id          aEntityId,
                enum EntityRank   aEntityRank,
                const moris_index aIndex ) const
        {
            auto tIter = mSTKMeshData->mEntityGlobaltoLocalMap( (uint)aEntityRank ).find( aEntityId );

            MORIS_ERROR( tIter != mSTKMeshData->mEntityGlobaltoLocalMap( (uint)aEntityRank ).end(),
                    "Provided Entity Id is not in the map, Has the map been initialized?: aEntityId =%u EntityRank = %u on process %u",
                    aEntityId,
                    (uint)aEntityRank,
                    par_rank() );

            return tIter->second;
        }

        // ----------------------------------------------------------------------------

        Matrix< IdMat >
        Mesh_Core_STK::get_entity_connected_to_entity_glob_ids(
                moris_id          aEntityId,
                enum EntityRank   aInputEntityRank,
                enum EntityRank   aOutputEntityRank,
                const moris_index aIndex ) const
        {

            // Call function that gets the connected entities
            stk::mesh::EntityRank tStkInputRank  = this->get_stk_entity_rank( aInputEntityRank );
            stk::mesh::EntityRank tStkOutputRank = this->get_stk_entity_rank( aOutputEntityRank );
            stk::mesh::Entity     tStkEntity     = mSTKMeshData->mMtkMeshBulkData->get_entity( tStkInputRank, (stk::mesh::EntityId)aEntityId );

            std::vector< stk::mesh::Entity > tSTKEntitiesConnectedGlobIds =
                    this->entities_connected_to_entity_stk( &tStkEntity, tStkInputRank, tStkOutputRank );

            // Get the number of entities
            uint tNumOutputEntities = tSTKEntitiesConnectedGlobIds.size();

            // Declare xtk::Mat that will contain the local entities
            Matrix< IndexMat > tEntitiesConnectedGlobIds( 1, tNumOutputEntities );

            // Fill local ids to xtk::Mat
            for ( uint i = 0; i < tNumOutputEntities; ++i )
            {
                tEntitiesConnectedGlobIds( 0, i ) = (moris_index)mSTKMeshData->mMtkMeshBulkData->identifier( tSTKEntitiesConnectedGlobIds[ i ] );
            }

            return tEntitiesConnectedGlobIds;
        }

        // ----------------------------------------------------------------------------
        void
        Mesh_Core_STK::get_elements_in_support_of_basis( const uint aMeshIndex,
                const uint                                          aBasisIndex,
                Matrix< IndexMat >&                                 aElementIndices )
        {
            aElementIndices = get_entity_connected_to_entity_loc_inds( aBasisIndex, EntityRank::NODE, EntityRank::ELEMENT );
        }

        // ----------------------------------------------------------------------------

        Matrix< IdMat >
        Mesh_Core_STK::get_element_connected_to_element_glob_ids( moris_id aElementId ) const
        {
            // First get faces connected to element
            // Get Stk entity Id from local to global map
            stk::mesh::EntityId tStkEntityId = (stk::mesh::EntityId)aElementId;

            // Call function that gets the connected entities
            stk::mesh::EntityRank tStkInputRank  = stk::topology::ELEMENT_RANK;
            stk::mesh::EntityRank tStkOutputRank = stk::topology::FACE_RANK;
            stk::mesh::Entity     tStkEntity     = mSTKMeshData->mMtkMeshBulkData->get_entity( tStkInputRank, tStkEntityId );

            std::vector< stk::mesh::Entity > tFacesInElem = this->entities_connected_to_entity_stk( &tStkEntity, tStkInputRank, tStkOutputRank );

            MORIS_ASSERT( ( tFacesInElem.size() != 0 ) || ( tFacesInElem.size() != 0 ),
                    "No faces connected to element found. Maybe the CreateAllEdgesAndFaces flag is set to false. Check mesh struct." );

            // Then for each face get elements connected
            uint tCounter  = 0;
            uint tNumFaces = tFacesInElem.size();

            moris::Matrix< IdMat > tElemsConnectedToElem( 2, tNumFaces );

            for ( uint faceIt = 0; faceIt < tNumFaces; ++faceIt )
            {
                std::vector< stk::mesh::Entity > tDummyConnectivity =
                        this->entities_connected_to_entity_stk( &tFacesInElem[ faceIt ], stk::topology::FACE_RANK, stk::topology::ELEMENT_RANK );

                // Faces in mesh boundaries do not have more than one element
                if ( tDummyConnectivity.size() > 0 )
                {
                    if ( mSTKMeshData->mMtkMeshBulkData->identifier( tDummyConnectivity[ 0 ] ) != mSTKMeshData->mMtkMeshBulkData->identifier( tStkEntity ) )
                    {
                        tElemsConnectedToElem( 0, tCounter ) = (moris_id)mSTKMeshData->mMtkMeshBulkData->identifier( tDummyConnectivity[ 0 ] );
                        tElemsConnectedToElem( 1, tCounter ) = (moris_id)faceIt;
                        tCounter++;
                    }
                }

                if ( tDummyConnectivity.size() > 1 )
                {
                    if ( mSTKMeshData->mMtkMeshBulkData->identifier( tDummyConnectivity[ 1 ] ) != mSTKMeshData->mMtkMeshBulkData->identifier( tStkEntity ) )
                    {
                        tElemsConnectedToElem( 0, tCounter ) = (moris_id)mSTKMeshData->mMtkMeshBulkData->identifier( tDummyConnectivity[ 1 ] );
                        tElemsConnectedToElem( 1, tCounter ) = (moris_id)faceIt;
                        tCounter++;
                    }
                }

                MORIS_ASSERT( tDummyConnectivity.size() <= 2,
                        "For some reason face has more than 2 elements connected to it... Check get_elements_connected_to_element." );
            }

            // Resize to include only ids added above and get rid of initialized extra zeros
            tElemsConnectedToElem.resize( 2, tCounter );

            return tElemsConnectedToElem;
        }

        // ----------------------------------------------------------------------------

        moris_index
        Mesh_Core_STK::get_facet_ordinal_from_cell_and_facet_id_glob_ids( moris_id aFaceId,
                moris_id                                                           aCellId ) const
        {
            MORIS_ASSERT( mSTKMeshData->mCreatedFaces, "Faces need to be created for this function" );

            Matrix< IdMat > tElementFaces = get_entity_connected_to_entity_glob_ids( aCellId, EntityRank::ELEMENT, this->get_facet_rank() );

            moris_index tOrdinal = MORIS_INDEX_MAX;
            for ( moris_index iOrd = 0; iOrd < (moris_index)tElementFaces.numel(); iOrd++ )
            {
                if ( tElementFaces( iOrd ) == aFaceId )
                {
                    tOrdinal = iOrd;
                    return tOrdinal;
                }
            }
            MORIS_ERROR( tOrdinal != MORIS_INDEX_MAX, " Facet ordinal not found" );
            return tOrdinal;
        }

        // ----------------------------------------------------------------------------

        Matrix< IdMat >
        Mesh_Core_STK::generate_unique_entity_ids(
                uint            aNumEntities,
                enum EntityRank aEntityRank ) const
        {
            std::vector< stk::mesh::EntityId > tAvailableNodeIDs;    // generate_new_ids requires a variable of this type
            mSTKMeshData->mMtkMeshBulkData->generate_new_ids( get_stk_entity_rank( aEntityRank ), aNumEntities, tAvailableNodeIDs );

            Matrix< IdMat > aAvailableNodeIDs( aNumEntities, 1 );

            for ( uint i = 0; i < aNumEntities; i++ )
            {
                aAvailableNodeIDs( i, 0 ) = (moris_id)tAvailableNodeIDs[ i ];
            }

            return aAvailableNodeIDs;
        }

        // ----------------------------------------------------------------------------

        moris_id
        Mesh_Core_STK::get_max_entity_id(
                enum EntityRank   aEntityRank,
                const moris_index aIndex ) const
        {
            MORIS_ASSERT( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT, "Only Elements or Nodes have max id" );

            moris::uint tNumEntities = this->get_num_entities( aEntityRank );

            moris_id tLocalMaxId = 0;

            for ( moris::uint i = 0; i < tNumEntities; i++ )
            {
                moris_id tId = this->get_glb_entity_id_from_entity_loc_index( i, aEntityRank );

                if ( tId > tLocalMaxId )
                {
                    tLocalMaxId = tId;
                }
            }

            moris_id tGlobalMaxId = moris::max_all( tLocalMaxId );
            return tGlobalMaxId;
        }

        // ----------------------------------------------------------------------------

        //##############################################
        // Coordinate Field Functions
        //##############################################

        Matrix< DDRMat >
        Mesh_Core_STK::get_node_coordinate( moris_index aNodeIndex ) const
        {
            // Get the coordinate field from stk
            stk::mesh::FieldBase const * coord = mSTKMeshData->mMtkMeshMetaData->coordinate_field();

            // Get node id from provided index
            moris_id tId = get_glb_entity_id_from_entity_loc_index( aNodeIndex, EntityRank::NODE );

            stk::mesh::EntityId tNodeId = (stk::mesh::EntityId)tId;

            // Declare node entity
            stk::mesh::Entity tNodeEntity = mSTKMeshData->mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, tNodeId );

            // Get coordinates of node n
            double* fieldValue = static_cast< double* >( stk::mesh::field_data( *coord, tNodeEntity ) );

            Matrix< DDRMat > tNodeCoord( 1, get_spatial_dim() );

            for ( uint dim = 0; dim < get_spatial_dim(); dim++ )
            {
                tNodeCoord( 0, dim ) = fieldValue[ dim ];
            }

            return tNodeCoord;
        }

        //##############################################
        // Entity Ownership Functions
        //##############################################

        uint
        Mesh_Core_STK::get_node_owner( moris_index aNodeIndex ) const
        {
            return this->get_entity_owner( aNodeIndex, EntityRank::NODE );
        }

        // ----------------------------------------------------------------------------

        uint
        Mesh_Core_STK::get_element_owner( moris_index aElementIndex ) const
        {
            return this->get_entity_owner( aElementIndex, EntityRank::ELEMENT );
        }

        // ----------------------------------------------------------------------------

        uint
        Mesh_Core_STK::get_entity_owner(
                moris_index       aEntityIndex,
                enum EntityRank   aEntityRank,
                const moris_index aIndex ) const
        {

            // Convert index to ID
            moris_id tEntityId = get_glb_entity_id_from_entity_loc_index( aEntityIndex, aEntityRank );

            // Get entity Id
            stk::mesh::Entity tEntity = mSTKMeshData->mMtkMeshBulkData->get_entity( get_stk_entity_rank( aEntityRank ), tEntityId );

            // processor rank that owns entity
            moris_id tOwningProcessor = mSTKMeshData->mMtkMeshBulkData->parallel_owner_rank( tEntity );

            return tOwningProcessor;
        }

        // ----------------------------------------------------------------------------

        Matrix< IdMat >
        Mesh_Core_STK::get_processors_whom_share_entity_glob_ids(
                moris_id        aEntityId,
                enum EntityRank aEntityRank ) const
        {
            // Convert index to ID
            stk::mesh::EntityId tEntityId = aEntityId;

            // Get entity
            stk::mesh::Entity tEntity = mSTKMeshData->mMtkMeshBulkData->get_entity( get_stk_entity_rank( aEntityRank ), tEntityId );

            // Intialize shared procs
            std::vector< int > tSharedProcs;

            // get shared processor IDs
            mSTKMeshData->mMtkMeshBulkData->comm_procs( tEntity, tSharedProcs );

            if ( tSharedProcs.size() == 0 )
            {
                tSharedProcs.push_back( mSTKMeshData->mMtkMeshBulkData->parallel_owner_rank( tEntity ) );
            }

            // Initialize output
            Matrix< IdMat > tProcsWhomShareEntity( 1, tSharedProcs.size() );

            // Cell to vector conversion
            for ( uint i = 0; i < tSharedProcs.size(); i++ )
            {
                tProcsWhomShareEntity( 0, i ) = tSharedProcs[ i ];
            }

            return tProcsWhomShareEntity;
        }

        // ----------------------------------------------------------------------------

        moris::Cell< std::string >
        Mesh_Core_STK::get_set_names( enum EntityRank aSetEntityRank ) const
        {
            const stk::mesh::PartVector& tMeshParts = mSTKMeshData->mMtkMeshMetaData->get_mesh_parts();

            stk::mesh::EntityRank tStkRank = get_stk_entity_rank( aSetEntityRank );

            moris::Cell< std::string > tSetNames;

            for ( moris::uint i = 0; i < tMeshParts.size(); i++ )
            {
                if ( tMeshParts[ i ]->primary_entity_rank() == tStkRank )
                {
                    tSetNames.push_back( tMeshParts[ i ]->name() );
                }
            }

            //        print(tSetNames,"tSetNames");
            // For whatever reason, the face sets have some additional internal sets
            //        if(aSetEntityRank == EntityRank::FACE)
            //        {
            //            moris::size_t  tNumParts = tSetNames.size()/2;
            //            std::string tDummy = "dummy";
            //            tSetNames.resize(tNumParts,tDummy);
            //        }
            return tSetNames;
        }

        // ----------------------------------------------------------------------------

        Matrix< IndexMat >
        Mesh_Core_STK::get_set_entity_loc_inds(
                enum EntityRank aSetEntityRank,
                std::string     aSetName ) const
        {
            // Get pointer to field defined by input name
            stk::mesh::Part* const tSetPart = mSTKMeshData->mMtkMeshMetaData->get_part( aSetName );

            MORIS_ASSERT( tSetPart != nullptr, "Set not found. Double check name provided." );

            // Access data through a selector
            stk::mesh::Selector     tSetSelector( *tSetPart );
            stk::mesh::EntityVector aEntities;
            stk::mesh::get_selected_entities( tSetSelector, mSTKMeshData->mMtkMeshBulkData->buckets( this->get_stk_entity_rank( aSetEntityRank ) ), aEntities );

            // Get entity Ids
            uint               tNumEntities = aEntities.size();
            Matrix< IndexMat > tOutputEntityInds( tNumEntities, 1 );
            for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
            {
                tOutputEntityInds( iEntity ) = this->get_loc_entity_ind_from_entity_glb_id( (moris_id)mSTKMeshData->mMtkMeshBulkData->identifier( aEntities[ iEntity ] ), aSetEntityRank );
            }

            return tOutputEntityInds;
        }

        // ----------------------------------------------------------------------------

        Matrix< IndexMat >
        Mesh_Core_STK::get_element_indices_in_block_set( uint aSetIndex )
        {
            // Get element IDs
            Matrix< IdMat > tElementIDs = this->get_element_ids_in_block_set( aSetIndex );

            // Element Indices
            Matrix< IndexMat > tElementIndices( tElementIDs.n_rows(), tElementIDs.n_cols() );

            // Transform into indices
            for ( uint tElementIndex = 0; tElementIndex < tElementIDs.length(); tElementIndex++ )
            {
                tElementIndices( tElementIndex ) = this->get_loc_entity_ind_from_entity_glb_id(
                        tElementIDs( tElementIndex ), EntityRank::ELEMENT );
            }

            return tElementIndices;
        }

        // ----------------------------------------------------------------------------

        Matrix< IdMat >
        Mesh_Core_STK::get_element_ids_in_block_set( uint aSetIndex )
        {
            // Get pointer to field defined by input name
            stk::mesh::Part* const tSetPart = mSTKMeshData->mMtkMeshMetaData->get_part(
                    this->get_set_names( EntityRank::ELEMENT )( aSetIndex ) );

            // Access data through a selector
            stk::mesh::Selector     tSetSelector( *tSetPart );
            stk::mesh::EntityVector aElements;

            stk::mesh::get_selected_entities(
                    tSetSelector,
                    mSTKMeshData->mMtkMeshBulkData->buckets(
                            this->get_stk_entity_rank( EntityRank::ELEMENT ) ),
                    aElements );

            // Get entity Ids
            uint            tNumElements = aElements.size();
            Matrix< IdMat > tElementIDs( tNumElements, 1 );

            for ( uint tElementIndex = 0; tElementIndex < tNumElements; tElementIndex++ )
            {
                tElementIDs( tElementIndex ) = (moris_id)mSTKMeshData->mMtkMeshBulkData->identifier( aElements[ tElementIndex ] );
            }

            return tElementIDs;
        }

        // ----------------------------------------------------------------------------

        enum CellTopology
        Mesh_Core_STK::get_blockset_topology( const std::string& aSetName )
        {
            // Get pointer to field defined by input name
            stk::mesh::Part* const tSetPart = mSTKMeshData->mMtkMeshMetaData->get_part( aSetName );

            // get part topology
            stk::topology::topology_t tTopology = tSetPart->topology();

            return stk_topo_to_moris_topo( tTopology );
        }

        // ----------------------------------------------------------------------------

        enum CellShape
        Mesh_Core_STK::get_IG_blockset_shape( const std::string& aSetName )
        {
            // get the clusters in the set
            moris::Cell< Cluster const * > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

            // init cell shape
            CellShape tCellShape = CellShape::EMPTY;

            // if the set isn't empty exist
            if ( tSetClusters.size() > 0 )
            {
                // get the cells in the first cluster
                moris::Cell< moris::mtk::Cell const * > tClusterCells = tSetClusters( 0 )->get_primary_cells_in_cluster();

                // compute the cell shape of the first cell
                tCellShape = tClusterCells( 0 )->get_cell_info()->compute_cell_shape( tClusterCells( 0 ) );
            }

            // checking all cells to make sure that they are the same Cell Shape; if not set cell shape to general
            // looping through the clusters
            for ( uint iCluster = 0; iCluster < tSetClusters.size(); iCluster++ )
            {
                // get cell of cells in the cluster
                moris::Cell< moris::mtk::Cell const * > tClusterCellsCheck = tSetClusters( iCluster )->get_primary_cells_in_cluster();

                // looping through the cells in the cluster
                for ( uint iCheckCell = 0; iCheckCell < tClusterCellsCheck.size(); iCheckCell++ )
                {
                    auto tCellInfo = tClusterCellsCheck( iCheckCell )->get_cell_info();

                    CellShape tCurrentShape = tCellInfo->compute_cell_shape( tClusterCellsCheck( iCheckCell ) );

                    if ( tCurrentShape == tCellShape )
                    {
                        tCellShape = CellShape::GENERAL;
                    }
                }
            }

            return tCellShape;
        }

        // ----------------------------------------------------------------------------

        enum CellShape
        Mesh_Core_STK::get_IP_blockset_shape( const std::string& aSetName )
        {
            // get the clusters in the set
            moris::Cell< Cluster const * > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

            // init cell shape
            CellShape tCellShape = CellShape::EMPTY;

            // if the set isn't empty exist
            if ( tSetClusters.size() > 0 )
            {
                // get the cells in the first cluster
                Cell const & tClusterCell = tSetClusters( 0 )->get_interpolation_cell();

                // compute the cell shape of the first cell
                tCellShape = tClusterCell.get_cell_info()->compute_cell_shape( &tClusterCell );
            }

            // checking all cells to make sure that they are the same Cell Shape; if not set cell shape to general
            // looping through the clusters
            for ( uint iCluster = 1; iCluster < tSetClusters.size(); iCluster++ )
            {
                auto tCellInfo = tSetClusters( iCluster )->get_interpolation_cell().get_cell_info();

                CellShape tCurrentShape = tCellInfo->compute_cell_shape( &tSetClusters( iCluster )->get_interpolation_cell() );

                if ( tCurrentShape == tCellShape )
                {
                    tCellShape = CellShape::GENERAL;
                }
            }

            return tCellShape;
        }

        // ----------------------------------------------------------------------------

        enum CellTopology
        Mesh_Core_STK::get_sideset_topology( const std::string& aSetName )
        {
            // Get pointer to field defined by input name
            stk::mesh::Part* const tSetPart = mSTKMeshData->mMtkMeshMetaData->get_part( aSetName );

            // get part topology
            stk::topology::topology_t tTopology = tSetPart->topology();

            return stk_topo_to_moris_topo( tTopology );
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::get_sideset_elems_loc_inds_and_ords(
                const std::string&  aSetName,
                Matrix< IndexMat >& aElemIndices,
                Matrix< IndexMat >& aSidesetOrdinals ) const
        {
            // Entity rank of a facet
            enum EntityRank tFacetRank = this->get_facet_rank();

            // Get facets in side set
            Matrix< IndexMat > tSidesInSideSetInds = this->get_set_entity_loc_inds( tFacetRank, aSetName );

            // allocate space in inputs
            aElemIndices.resize( 1, tSidesInSideSetInds.numel() * 2 );
            aSidesetOrdinals.resize( 1, tSidesInSideSetInds.numel() * 2 );

            // loop through sides and get cells attached, then figure out the ordinal
            uint tCount = 0;

            for ( moris::uint i = 0; i < tSidesInSideSetInds.numel(); i++ )
            {
                Matrix< IndexMat > tElementConnectedToFace =
                        this->get_entity_connected_to_entity_loc_inds( tSidesInSideSetInds( i ), tFacetRank, EntityRank::ELEMENT );

                // iterate through the cells
                for ( moris::uint j = 0; j < tElementConnectedToFace.numel(); j++ )
                {
                    aElemIndices( tCount ) = tElementConnectedToFace( j );

                    aSidesetOrdinals( tCount ) =
                            this->get_facet_ordinal_from_cell_and_facet_loc_inds( tSidesInSideSetInds( i ), tElementConnectedToFace( j ) );

                    tCount++;
                }
            }

            aElemIndices.resize( 1, tCount );
            aSidesetOrdinals.resize( 1, tCount );
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::get_sideset_cells_and_ords(
                const std::string&                aSetName,
                moris::Cell< mtk::Cell const * >& aCells,
                Matrix< IndexMat >&               aSidesetOrdinals ) const
        {
            Matrix< IndexMat > tCellInds;

            this->get_sideset_elems_loc_inds_and_ords( aSetName, tCellInds, aSidesetOrdinals );

            // moris::print(tCellInds,"tCellInds");

            aCells.resize( tCellInds.numel() );

            // iterate through cell inds and get cell ptrs
            for ( moris::uint i = 0; i < tCellInds.numel(); i++ )
            {
                aCells( i ) = &this->get_mtk_cell( tCellInds( i ) );
            }
        }

        // ----------------------------------------------------------------------------
        moris::Cell< moris::mtk::Vertex const * >
        Mesh_Core_STK::get_vertices_in_vertex_set_no_aura( std::string aSetName ) const
        {
            enum EntityRank tEntityRank = EntityRank::NODE;

            // Get pointer to field defined by input name
            // set and set selector
            stk::mesh::Part* const tSetPart = mSTKMeshData->mMtkMeshMetaData->get_part( aSetName );
            MORIS_ASSERT( tSetPart != nullptr, "Set not found. Double check name provided." );
            stk::mesh::Selector tSetSelector( *tSetPart );

            // aura and aura selector
            stk::mesh::Part&    tAura = mSTKMeshData->mMtkMeshMetaData->aura_part();
            stk::mesh::Selector tAuraSelector( tAura );

            // difference of the set with the auraz
            stk::mesh::Selector tDifference = tSetSelector -= tAuraSelector;

            stk::mesh::EntityVector aEntities;
            stk::mesh::get_selected_entities( tDifference, mSTKMeshData->mMtkMeshBulkData->buckets( stk::topology::NODE_RANK ), aEntities );

            // Get entity Ids
            uint                                      tNumEntities = aEntities.size();
            moris::Cell< moris::mtk::Vertex const * > tOutputEntityIds( tNumEntities );
            for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
            {
                moris::moris_index tVertexIndex =
                        this->get_loc_entity_ind_from_entity_glb_id( (moris_id)mSTKMeshData->mMtkMeshBulkData->identifier( aEntities[ iEntity ] ), tEntityRank );

                tOutputEntityIds( iEntity ) = &this->get_mtk_vertex( tVertexIndex );
            }

            return tOutputEntityIds;
        }

        // ----------------------------------------------------------------------------

        uint
        Mesh_Core_STK::get_num_fields(
                const enum EntityRank aEntityRank,
                const moris_index     aIndex ) const
        {
            stk::mesh::EntityRank         tEntityRank  = this->get_stk_entity_rank( aEntityRank );
            const stk::mesh::FieldVector& tFieldVector = mSTKMeshData->mMtkMeshMetaData->get_fields( tEntityRank );
            return tFieldVector.size();
        }

        // ----------------------------------------------------------------------------

        moris_index
        Mesh_Core_STK::get_field_ind(
                const std::string&    aFieldLabel,
                const enum EntityRank aEntityRank ) const
        {
            stk::mesh::EntityRank         tEntityRank  = this->get_stk_entity_rank( aEntityRank );
            const stk::mesh::FieldVector& tFieldVector = mSTKMeshData->mMtkMeshMetaData->get_fields( tEntityRank );

            moris_index tFieldIndex = MORIS_INDEX_MAX;
            for ( moris::uint i = 0; i < tFieldVector.size(); i++ )
            {
                if ( tFieldVector[ i ]->name().compare( aFieldLabel ) == 0 )
                {
                    tFieldIndex = i;
                    break;
                }
            }

            MORIS_ERROR( tFieldIndex != MORIS_INDEX_MAX, "Field not found in get_field_ind()" );
            return tFieldIndex;
        }

        // ----------------------------------------------------------------------------

        Matrix< DDRMat >
        Mesh_Core_STK::get_entity_field_value_real_scalar(
                Matrix< IndexMat > const & aEntityIndices,
                std::string const &        aFieldName,
                enum EntityRank            aFieldEntityRank ) const
        {
            // MORIS_ASSERT(aFieldEntityRank==EntityRank::NODE,"Only implemented for nodal scalar field");

            // Initialize Output
            size_t           tNumEntities = aEntityIndices.n_cols();
            Matrix< DDRMat > tFieldValues( 1, tNumEntities );

            // Get field by name and entity rank
            stk::mesh::EntityRank tEntityRank = this->get_stk_entity_rank( aFieldEntityRank );

            stk::mesh::Field< real >* tField =
                    mSTKMeshData->mMtkMeshMetaData->get_field< stk::mesh::Field< real > >( tEntityRank, aFieldName );

            // make sure that field actually exists
            if ( tField == NULL )
            {
                // select specifier for rank
                std::string tRank;

                switch ( aFieldEntityRank )
                {
                    case ( moris::EntityRank::NODE ):
                    {
                        tRank = " node ";
                        break;
                    }
                    case ( moris::EntityRank::ELEMENT ):
                    {
                        tRank = " element ";
                        break;
                    }
                    default:
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
            for ( size_t i = 0; i < tNumEntities; i++ )
            {
                moris_id          tId     = get_glb_entity_id_from_entity_loc_index( aEntityIndices( 0, i ), aFieldEntityRank );
                stk::mesh::Entity tEntity = mSTKMeshData->mMtkMeshBulkData->get_entity( tEntityRank, tId );

                // Store the coordinates of the current node
                real* tFieldData     = stk::mesh::field_data( *tField, tEntity );
                tFieldValues( 0, i ) = tFieldData[ 0 ];
            }

            return tFieldValues;
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::add_mesh_field_real_scalar_data_loc_inds(
                std::string const &      aFieldName,
                enum EntityRank const &  aFieldEntityRank,
                Matrix< DDRMat > const & aFieldData )
        {

            MORIS_ASSERT( aFieldEntityRank == EntityRank::NODE || aFieldEntityRank == EntityRank::ELEMENT,
                    "Only tested for nodal and element scalar field" );

            // Write Data to Field
            size_t tNumEntities = get_num_entities( aFieldEntityRank );

            // Get Field
            stk::mesh::EntityRank     tEntityRank = this->get_stk_entity_rank( aFieldEntityRank );
            stk::mesh::Field< real >* tField      = mSTKMeshData->mMtkMeshMetaData->get_field< stk::mesh::Field< real > >( tEntityRank, aFieldName );
            for ( size_t i = 0; i < tNumEntities; i++ )
            {
                // Get global Id of current node and create "node entity" for stk mesh
                // stk::mesh::EntityId nodeGlobalId = node_i;
                moris_id          tId     = get_glb_entity_id_from_entity_loc_index( i, aFieldEntityRank );
                stk::mesh::Entity tEntity = mSTKMeshData->mMtkMeshBulkData->get_entity( tEntityRank, tId );

                // Store the coordinates of the current node
                real* tFieldData = stk::mesh::field_data( *tField, tEntity );

                tFieldData[ 0 ] = aFieldData( i );
            }
        }

        // ----------------------------------------------------------------------------

        //##############################################
        // moris::Cell and Vertex Pointer Functions
        //##############################################

        mtk::Cell&
        Mesh_Core_STK::get_mtk_cell( moris_index aCellIndex )
        {
            MORIS_ASSERT( aCellIndex < (moris_index)mSTKMeshData->mMtkCells.size(),
                    "Provided cell index out of bounds: aCellIndex = %u , mSTKMeshData->mMtkCells.size() = %u ",
                    aCellIndex,
                    mSTKMeshData->mMtkCells.size() );

            return mSTKMeshData->mMtkCells( aCellIndex );
        }

        // ----------------------------------------------------------------------------

        mtk::Cell const &
        Mesh_Core_STK::get_mtk_cell( moris_index aCellIndex ) const
        {
            MORIS_ASSERT( aCellIndex < (moris_index)mSTKMeshData->mMtkCells.size(),
                    "Provided cell index out of bounds: aCellIndex = %u , mSTKMeshData->mMtkCells.size() = %u ",
                    aCellIndex,
                    mSTKMeshData->mMtkCells.size() );

            return mSTKMeshData->mMtkCells( aCellIndex );
        }

        // ----------------------------------------------------------------------------

        mtk::Vertex&
        Mesh_Core_STK::get_mtk_vertex( moris_index aVertexIndex )
        {
            return mSTKMeshData->mMtkVertices( aVertexIndex );
        }

        // ----------------------------------------------------------------------------

        mtk::Vertex const &
        Mesh_Core_STK::get_mtk_vertex( moris_index aVertexIndex ) const
        {
            return mSTKMeshData->mMtkVertices( aVertexIndex );
        }

        // ----------------------------------------------------------------------------

        mtk::Vertex_Core_STK&
        Mesh_Core_STK::get_mtk_vertex_stk( moris_index aVertexIndex )
        {
            return mSTKMeshData->mMtkVertices( aVertexIndex );
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::create_output_mesh(
                std::string& aFileName,
                bool         aAddElemCmap )
        {
            // start timing
            std::clock_t start = std::clock();

            if ( mSTKMeshData->mDataGeneratedMesh )
            {
                // Generate data for mesh from mesh reader
                size_t outputFileIdx = mSTKMeshData->mMeshReader->create_output_mesh( aFileName, stk::io::WRITE_RESULTS );
                mSTKMeshData->mMeshReader->write_output_mesh( outputFileIdx );

                // Get fields and initialize fields information
                const std::vector< stk::mesh::FieldBase* >& fields = mSTKMeshData->mMtkMeshMetaData->get_fields();

                // Add fields to output mesh
                std::string tFieldNoData = "dummyField";
                std::string tCoordField  = "coordinates";

                std::vector< stk::mesh::FieldBase* >::const_iterator fieldIterator = fields.begin();

                for ( ; fieldIterator != fields.end(); ++fieldIterator )
                {
                    // Get field name
                    std::string tIterFieldName = ( *fieldIterator )->name();
                    // Do not add dummy or coordinate fields to the output mesh
                    if ( ( tIterFieldName.compare( tFieldNoData ) != 0 ) && ( tIterFieldName.compare( tCoordField ) != 0 ) )
                    {
                        mSTKMeshData->mMeshReader->add_field(
                                outputFileIdx,
                                *( stk::mesh::get_field_by_name(
                                        ( *fieldIterator )->name(),
                                        *mSTKMeshData->mMtkMeshMetaData ) ) );
                    }
                }

                // Provisionally only handles static problems (hard-coded time)
                mSTKMeshData->mMeshReader->begin_output_step( outputFileIdx, mTimeStamp );
                mSTKMeshData->mMeshReader->write_defined_output_fields( outputFileIdx );
                mSTKMeshData->mMeshReader->end_output_step( outputFileIdx );
            }
            else
            {
                // Generate data for mesh from mesh reader
                size_t fh = mSTKMeshData->mMeshReader->create_output_mesh( aFileName, stk::io::WRITE_RESULTS );

                // write mesh with the information generated from the mesh reader
                mSTKMeshData->mMeshReader->write_output_mesh( fh );

                // Add fields that weren't in the file loaded in
                for ( auto iField : mSTKMeshData->mRealNodeScalarFieldsToAddToOutput )
                {
                    mSTKMeshData->mMeshReader->add_field( fh, *iField );
                }

                mSTKMeshData->mMeshReader->begin_output_step( fh, mTimeStamp );
                mSTKMeshData->mMeshReader->write_defined_output_fields( fh );
                mSTKMeshData->mMeshReader->end_output_step( fh );
            }

            if ( aAddElemCmap && par_size() > 1 )
            {
                Exodus_IO_Helper tMTKExoIO( aFileName.c_str() );
                add_element_cmap_to_exodus( aFileName, tMTKExoIO );
            }

            if ( mVerbose )
            {
                MORIS_LOG_INFO( "MTK: Exodus output completed in %f s", ( std::clock() - start ) / (double)( CLOCKS_PER_SEC ) );
                MORIS_LOG_INFO( "MTK: Exodus file: %s", aFileName.c_str() );
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::add_element_cmap_to_exodus(
                std::string&      aFileName,
                Exodus_IO_Helper& aExoIO )
        {

            MORIS_ASSERT( mSTKMeshData->mCreatedFaces, "Element CMap requires faces to be created" );

            enum EntityRank tFacetRank = this->get_facet_rank();
            // Get the sides on processor boundaries
            stk::mesh::EntityVector tSidesOnProcBoundaries;
            stk::mesh::get_selected_entities( mSTKMeshData->mMtkMeshMetaData->globally_shared_part(),
                    mSTKMeshData->mMtkMeshBulkData->buckets( mSTKMeshData->mMtkMeshMetaData->side_rank() ),
                    tSidesOnProcBoundaries );

            // Number of sides on mesh along processor boundaries.
            uint tNumSidesOnProcBoundaries = tSidesOnProcBoundaries.size();

            // Allocate vector for element ids attached to sides and side ordinals
            // attached to the (note: we limit the side to be attached to one element along the boundary)
            Matrix< IdMat > tElementIdsOnBoundaries( 1, tNumSidesOnProcBoundaries );
            Matrix< IdMat > tSideOrdinalsOnBoundaries( 1, tNumSidesOnProcBoundaries );
            Matrix< IdMat > tSideSharedProc( 1, tNumSidesOnProcBoundaries );

            // Allocate vector for elements connected to a given faces
            Matrix< IdMat > tElementToFace( 1, 2 );
            Matrix< IdMat > tSharedProcessorsOfFace( 1, 2 );
            uint            tProcRank = par_rank();

            // Populate the above vectors
            for ( uint iS = 0; iS < tNumSidesOnProcBoundaries; iS++ )
            {
                // Get the side id
                moris_id tSideId = mSTKMeshData->mMtkMeshBulkData->identifier( tSidesOnProcBoundaries[ iS ] );

                // Get the elements attached to this faces
                tElementToFace = this->get_entity_connected_to_entity_glob_ids( tSideId, tFacetRank, EntityRank::ELEMENT );

                // Figure out which other processors share this sides
                tSharedProcessorsOfFace = this->get_processors_whom_share_entity_glob_ids( tSideId, tFacetRank );
                for ( uint iP = 0; iP < tSharedProcessorsOfFace.numel(); iP++ )
                {
                    if ( tSharedProcessorsOfFace( iP ) != (moris_id)tProcRank )
                    {
                        tSideSharedProc( iS ) = tSharedProcessorsOfFace( iP );
                        break;
                    }
                }    // shared proc loop

                // iterate through elements connected to faces
                for ( uint iE = 0; iE < tElementToFace.numel(); iE++ )
                {
                    // get the element id from the vector
                    moris_id tElementId = tElementToFace( iE );

                    // check whether the element is in the aura
                    //                if(!this->is_aura_cell(tElementId))
                    //                {
                    // figure out the side ordinal of this element where tSideId is
                    moris_id tSideOrdinal = this->get_facet_ordinal_from_cell_and_facet_id_glob_ids( tSideId, tElementId );

                    // Add to vectors
                    tElementIdsOnBoundaries( iS )   = tElementId;
                    tSideOrdinalsOnBoundaries( iS ) = tSideOrdinal;
                    //                }
                }    // Element loop
            }        // side loop

            // Query the exodus file about comm maps

            aExoIO.create_new_exo_with_elem_cmaps_from_existing_exo( aFileName,
                    tElementIdsOnBoundaries,
                    tSideOrdinalsOnBoundaries,
                    tSideSharedProc );
        }

        //##############################################
        // Private functions to build mesh
        //##############################################
        void
        Mesh_Core_STK::create_communication_lists_and_local_to_global_map( enum EntityRank aEntityRank )
        {
            const int tParallelSize = mSTKMeshData->mMtkMeshBulkData->parallel_size();
            const int tParallelRank = mSTKMeshData->mMtkMeshBulkData->parallel_rank();

            // Declare vector of entity counts
            std::vector< long unsigned int > tEntityCounts;

            // Get all entities from meta data
            stk::mesh::Selector tSharedSelector = mSTKMeshData->mMtkMeshMetaData->universal_part();

            // Count entities
            stk::mesh::count_entities( tSharedSelector, *mSTKMeshData->mMtkMeshBulkData, tEntityCounts );

            uint tNumEntities = static_cast< uint >( tEntityCounts[ (uint)aEntityRank ] );

            // Resize comm lists to maximum possible
            for ( int i = 0; i < tParallelSize; i++ )
            {
                moris::Matrix< IndexMat > tSendMat( 1, tNumEntities, (uint)0 );
                moris::Matrix< IndexMat > tRecvMat( 1, tNumEntities, (uint)0 );

                mSTKMeshData->mEntitySendList( (uint)aEntityRank )( i )    = tSendMat;
                mSTKMeshData->mEntityReceiveList( (uint)aEntityRank )( i ) = tRecvMat;
            }

            mSTKMeshData->mEntityLocaltoGlobalMap( (uint)aEntityRank ) = moris::Matrix< IndexMat >( 1, tNumEntities, (moris_index)0 );

            stk::mesh::BucketVector const & shared_node_buckets =
                    mSTKMeshData->mMtkMeshBulkData->get_buckets( get_stk_entity_rank( aEntityRank ), tSharedSelector );

            uint tCurrentIndex = 0;
            // Initialize proc counter
            // moris::Cell #  = Proc rank
            moris::Cell< uint > tSendProcCounter( tParallelSize );
            moris::Cell< uint > tRecvProcCounter( tParallelSize );

            // Loop over shared nodes
            for ( uint i = 0; i < shared_node_buckets.size(); i++ )
            {
                stk::mesh::Bucket& bucket = *shared_node_buckets[ i ];

                for ( uint j = 0; j < bucket.size(); j++ )
                {
                    moris_id tEntityId      = (moris_id)mSTKMeshData->mMtkMeshBulkData->identifier( bucket[ j ] );
                    int      tOwnerProcRank = mSTKMeshData->mMtkMeshBulkData->parallel_owner_rank( bucket[ j ] );

                    // Set local to global map in mesh and STK
                    mSTKMeshData->mEntityLocaltoGlobalMap( (uint)aEntityRank )( 0, tCurrentIndex ) = tEntityId;

                    mSTKMeshData->mMtkMeshBulkData->set_local_id( bucket[ j ], tCurrentIndex );

                    std::vector< int > sharedProcs;

                    // Get shared procs Ids
                    mSTKMeshData->mMtkMeshBulkData->comm_procs( bucket[ j ], sharedProcs );

                    if ( sharedProcs.size() != 0 )
                    {

                        // Sort (if current proc owns entity then add to send comm lists)
                        if ( tOwnerProcRank == tParallelRank )
                        {
                            // loop over processors
                            for ( uint p = 0; p < sharedProcs.size(); p++ )
                            {
                                if ( sharedProcs[ p ] != mSTKMeshData->mMtkMeshBulkData->parallel_rank() )
                                {
                                    uint tSharedProcRank                                                                   = sharedProcs[ p ];
                                    uint tSendCount                                                                        = tSendProcCounter( tSharedProcRank );
                                    mSTKMeshData->mEntitySendList( (uint)aEntityRank )( tSharedProcRank )( 0, tSendCount ) = tEntityId;
                                    tSendProcCounter( tSharedProcRank )++;
                                }
                            }
                        }
                        // (if current proc does not own entity then add to recv comm lists)
                        else if ( tOwnerProcRank != tParallelRank )
                        {
                            for ( uint p = 0; p < sharedProcs.size(); p++ )
                            {
                                if ( sharedProcs[ p ] != tParallelRank )
                                {
                                    uint tSharedProc                                                                      = sharedProcs[ p ];
                                    uint tRecvCount                                                                       = tRecvProcCounter( tSharedProc );
                                    mSTKMeshData->mEntityReceiveList( (uint)aEntityRank )( tSharedProc )( 0, tRecvCount ) = tEntityId;
                                    tRecvProcCounter( tSharedProc )++;
                                }
                            }
                        }
                    }

                    tCurrentIndex++;
                }
            }

            for ( int pr = 0; pr < tParallelSize; pr++ )
            {
                uint tRecvCount = tRecvProcCounter( pr );
                uint tSendCount = tSendProcCounter( pr );
                mSTKMeshData->mEntitySendList( (uint)aEntityRank )( pr ).resize( 1, tSendCount );
                mSTKMeshData->mEntityReceiveList( (uint)aEntityRank )( pr ).resize( 1, tRecvCount );
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::setup_parallel_vertex_pairing()
        {
            int tParSize = par_size();

            if ( tParSize > 1 )
            {
                moris::Cell< moris::Matrix< IdMat > > tUniqueNodeSharedIds( tParSize );
                moris::Cell< moris::Matrix< IdMat > > tReceivedNodeSharedIds( tParSize );

                // iterate through procs and gather shared nodes between them
                for ( int pr = 0; pr < tParSize; pr++ )
                {
                    // Sizes
                    moris::uint tNumSendToProc  = mSTKMeshData->mEntitySendList( 0 )( pr ).numel();
                    moris::uint tNumRecFromProc = mSTKMeshData->mEntityReceiveList( 0 )( pr ).numel();

                    // compile a list of all ids
                    Matrix< IdMat > tNodesSharedWithProc( tNumSendToProc + tNumRecFromProc, 1 );

                    // count
                    moris::uint tCount = 0;

                    for ( moris::uint i = 0; i < tNumSendToProc; i++ )
                    {
                        tNodesSharedWithProc( tCount ) = mSTKMeshData->mEntitySendList( 0 )( pr )( i );
                        tCount++;
                    }

                    for ( moris::uint i = 0; i < tNumRecFromProc; i++ )
                    {
                        tNodesSharedWithProc( tCount ) = mSTKMeshData->mEntityReceiveList( 0 )( pr )( i );
                        tCount++;
                    }

                    // use unique here to just to be safe, although I do not foresee there is a case
                    // where a vertex shows up in both send and receive.
                    moris::unique( tNodesSharedWithProc, tUniqueNodeSharedIds( pr ) );

                    // at this point, we have a unique node list which should be in ascending order (unique sorts the matrix)

                    // allocate buffer for received communication
                    tReceivedNodeSharedIds( pr ).resize( tUniqueNodeSharedIds( pr ).n_rows(), 2 );
                }

                // add the processor local index to the communication, also remove unused procs
                mSTKMeshData->mProcsWithSharedVertex.resize( tParSize, 1 );
                moris::uint tCount = 0;
                for ( int pr = 0; pr < tParSize; pr++ )
                {
                    moris::uint tNumShared = tUniqueNodeSharedIds( pr ).numel();

                    if ( tNumShared > 0 )
                    {
                        mSTKMeshData->mProcsWithSharedVertex( tCount ) = (moris_id)pr;
                        mSTKMeshData->mVertexSharingProcsMap[ pr ]     = tCount;
                        tCount++;
                        tUniqueNodeSharedIds( pr ).resize( tNumShared, 2 );

                        // iterate through and get local index
                        for ( moris::uint i = 0; i < tNumShared; i++ )
                        {
                            tUniqueNodeSharedIds( pr )( i, 1 ) = this->get_loc_entity_ind_from_entity_glb_id( tUniqueNodeSharedIds( pr )( i, 0 ), EntityRank::NODE );
                        }
                    }
                }

                mSTKMeshData->mProcsWithSharedVertex.resize( tCount, 1 );

                // remove the empty comms
                tUniqueNodeSharedIds.data().erase( std::remove_if( tUniqueNodeSharedIds.begin(), tUniqueNodeSharedIds.end(), []( const Matrix< IdMat >& o ) { return o.numel() == 0; } ), tUniqueNodeSharedIds.end() );

                moris::communicate_mats( mSTKMeshData->mProcsWithSharedVertex, tUniqueNodeSharedIds, mSTKMeshData->mVertexSharingData );
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::setup_parallel_cell_sharing()
        {
            // construct sharing information as seen in STK
            setup_parallel_cell_sharing_without_aura_resolved();

            // get the aura cell sharing that need special considerations (cells that show up in multiple process aura require extra care to work
            // with XTK and cell clustering concept
            moris::Cell< Matrix< IdMat > > tAuraCellSharing = resolve_aura_cell_sharing();

            setup_parallel_cell_sharing_with_resolved_aura( tAuraCellSharing );
        }

        // ----------------------------------------------------------------------------

        moris::Cell< Matrix< IdMat > >
        Mesh_Core_STK::resolve_aura_cell_sharing()
        {
            moris_index tParSize = par_size();
            moris_index tParRank = par_rank();

            // received cell sharing information
            moris::Cell< Matrix< IdMat > > tReceivedCellSharing;

            if ( tParSize > 1 )
            {
                // number of cells in the mesh
                uint tNumCells = this->get_num_entities( EntityRank::ELEMENT );

                // allocate space in member data
                mSTKMeshData->mCellSharingData.resize( tNumCells );

                // locally owned  and local own selector
                stk::mesh::Part&    tLocalOwnedCells = mSTKMeshData->mMtkMeshMetaData->locally_owned_part();
                stk::mesh::Selector tLocalOwnedSelector( tLocalOwnedCells );

                // get cells in the local own part
                stk::mesh::EntityVector aEntities;
                stk::mesh::get_selected_entities( tLocalOwnedSelector, mSTKMeshData->mMtkMeshBulkData->buckets( stk::topology::ELEMENT_RANK ), aEntities );

                // cells to communicate additionally sharing
                moris::Cell< moris::Cell< Matrix< IdMat > > > tCellsWithAdditionalSharing( tParSize );

                // keep track of the maximum so we can convert the above cell to a
                moris::Cell< uint > tMaxSizeWithProc( tParSize, 0 );

                // procs to communicate with
                moris::Cell< moris_index > tProcsWithCommunication;

                // iterate through cells in aura
                for ( moris_index iCell = 0; iCell < (moris_index)aEntities.size(); iCell++ )
                {
                    // cell id
                    moris_id tCellId = (moris_id)mSTKMeshData->mMtkMeshBulkData->identifier( aEntities[ iCell ] );

                    // Cell index
                    moris_index tCellIndex = this->get_loc_entity_ind_from_entity_glb_id( tCellId, EntityRank::ELEMENT );

                    // get process rank that owns cell
                    uint tCellOwnerRank = this->get_entity_owner( tCellIndex, EntityRank::ELEMENT );

                    // sharing processors
                    Matrix< IdMat > tSharingProcs;
                    this->get_processors_whom_share_entity( tCellIndex, EntityRank::ELEMENT, tSharingProcs );

                    // check for the case where this cell I own that shows up in a few processors aura
                    if ( tParRank == (moris_index)tCellOwnerRank && tSharingProcs.numel() > 0 )
                    {
                        Matrix< IdMat > tSharingInfo( 1, tSharingProcs.numel() + 1 );

                        // add cell id to first spot
                        tSharingInfo( 0 )                                      = tCellId;
                        tSharingInfo( { 0, 0 }, { 1, tSharingProcs.numel() } ) = tSharingProcs.get_row( 0 );

                        // iterate through the shared processors
                        for ( moris_index iProc = 0; iProc < (moris_index)tSharingProcs.numel(); iProc++ )
                        {
                            // sharing processor rank
                            moris_index tShareProcRank = tSharingProcs( iProc );

                            // keep track of largest number
                            if ( tSharingProcs.numel() > tMaxSizeWithProc( tShareProcRank ) )
                            {
                                // add cell sharing information
                                tCellsWithAdditionalSharing( tShareProcRank ).push_back( tSharingInfo );

                                // if this is the first time we encountered this proc add it to the list of processors to communicate with
                                if ( tMaxSizeWithProc( tShareProcRank ) == 0 )
                                {
                                    tProcsWithCommunication.push_back( tShareProcRank );
                                }

                                // change max size
                                tMaxSizeWithProc( tShareProcRank ) = tSharingProcs.numel();
                            }
                        }
                    }
                }

                // create a cell of matrix rather than a cell of cells for communication purposes
                moris::Cell< Matrix< IdMat > > tCellsMatsWithAdditionalSharing( tProcsWithCommunication.size() );

                Matrix< IdMat > tProcsWithCommunicationMat( 1, tProcsWithCommunication.size() );

                for ( moris_index iProc = 0; iProc < (moris_index)tProcsWithCommunication.size(); iProc++ )
                {
                    // processor rank
                    moris_index tShareProcRank = tProcsWithCommunication( iProc );

                    // add to matrix of processor ranks
                    tProcsWithCommunicationMat( iProc ) = tShareProcRank;

                    // number of columns
                    uint tNumCol = tMaxSizeWithProc( tShareProcRank ) + 1;

                    // number of rows
                    uint tNumCellsWithAddShare = tCellsWithAdditionalSharing( tShareProcRank ).size();

                    // initialize matrix with a dummy value
                    tCellsMatsWithAdditionalSharing( iProc ) = Matrix< IdMat >( tNumCellsWithAddShare, tNumCol );
                    tCellsMatsWithAdditionalSharing( iProc ).fill( MORIS_INDEX_MAX );

                    // iterate through and add data to matrix
                    for ( uint iCell = 0; iCell < tNumCellsWithAddShare; iCell++ )
                    {
                        // number of sharing (size of matrix)
                        uint tSizeOfMat = tCellsWithAdditionalSharing( tShareProcRank )( iCell ).numel();

                        tCellsMatsWithAdditionalSharing( iProc )( { iCell, iCell }, { 0, tSizeOfMat - 1 } ) = tCellsWithAdditionalSharing( tShareProcRank )( iCell ).get_row( 0 );
                    }
                }

                // at this point we have collected the sharing information needed
                // so time to communicate
                moris::communicate_mats( tProcsWithCommunicationMat, tCellsMatsWithAdditionalSharing, tReceivedCellSharing );
            }

            return tReceivedCellSharing;
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::setup_parallel_cell_sharing_without_aura_resolved()
        {

            enum EntityRank tEntityRank = EntityRank::ELEMENT;

            // number of cells in mesh
            uint tNumCells = this->get_num_entities( tEntityRank );

            // allocate member data
            mSTKMeshData->mCellSharingData = moris::Cell< moris::Matrix< moris::IdMat > >( tNumCells );

            // iterate through cells and add sharing information
            for ( moris::uint i = 0; i < tNumCells; i++ )
            {

                // Convert index to ID
                stk::mesh::EntityId tEntityId = { (stk::mesh::EntityId)this->get_glb_entity_id_from_entity_loc_index( (moris_index)i, tEntityRank ) };

                // Get entity
                stk::mesh::Entity tEntity = mSTKMeshData->mMtkMeshBulkData->get_entity( get_stk_entity_rank( tEntityRank ), tEntityId );

                // Initialize shared procs
                std::vector< int > tSharedProcs;

                // get shared processor IDs
                mSTKMeshData->mMtkMeshBulkData->comm_procs( tEntity, tSharedProcs );

                if ( tSharedProcs.size() > 0 )
                {
                    // Initialize output
                    Matrix< IdMat > tSharingProcs( 1, tSharedProcs.size() );

                    // Cell to vector conversion
                    for ( uint j = 0; j < tSharedProcs.size(); j++ )
                    {
                        tSharingProcs( 0, j ) = tSharedProcs[ j ];
                    }

                    mSTKMeshData->mCellSharingData( i ) = tSharingProcs.copy();
                }
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::setup_parallel_cell_sharing_with_resolved_aura( moris::Cell< Matrix< IdMat > > const & aAuraCellSharing )
        {
            // my processor rank
            moris_id tMyProcRank = par_rank();

            // iterate through cell of matrices where each cell if from a different processor
            for ( moris::uint iProc = 0; iProc < aAuraCellSharing.size(); iProc++ )
            {
                // All cell sharing to resolve for this processor (note: which processor it came from does not matter)
                Matrix< IdMat > const & tCellSharingMat = aAuraCellSharing( iProc );

                // number of cols
                uint tNumCols = tCellSharingMat.n_cols();

                // iterate through the received information and update the cell sharing information
                for ( moris::uint iCell = 0; iCell < tCellSharingMat.n_rows(); iCell++ )
                {
                    // cell id
                    moris_id tCellId = tCellSharingMat( iCell, 0 );

                    // cell index
                    moris_index tCellIndex = this->get_loc_entity_ind_from_entity_glb_id( tCellId, EntityRank::ELEMENT );

                    // owning processor
                    moris_id tOwner = this->get_entity_owner( tCellIndex, EntityRank::ELEMENT );

                    // construct sharing matrix
                    Matrix< IdMat > tCellSharing( 1, tNumCols );

                    // add owner in first spot
                    tCellSharing( 0, 0 ) = tOwner;

                    MORIS_ASSERT( tOwner != tMyProcRank,
                            "Current processor should not be informed about aura, it should be the one telling other processors about sharing " );

                    // iterate through and construct
                    uint tCount = 1;

                    for ( moris::uint iShare = 0; iShare < tNumCols - 1; iShare++ )
                    {
                        // Process cell is shared with
                        moris_id tCellSharedWith = tCellSharingMat( iCell, iShare + 1 );

                        if ( tCellSharedWith == MORIS_ID_MAX )
                        {
                            break;
                        }

                        else if ( tCellSharedWith != tMyProcRank )
                        {
                            tCellSharing( 0, tCount ) = tCellSharedWith;
                            tCount++;
                        }
                    }

                    // size out extra space
                    tCellSharing.resize( 1, tCount );

                    // update cell sharing table
                    mSTKMeshData->mCellSharingData( tCellIndex ) = tCellSharing.copy();
                }
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::setup_vertices_and_cell()
        {
            // Get information about the mesh
            uint tNumNodes = this->get_num_entities( EntityRank::NODE );

            // Setup vertices
            mSTKMeshData->mMtkVertices = moris::Cell< Vertex_Core_STK >( tNumNodes );
            for ( uint iVertInd = 0; iVertInd < tNumNodes; iVertInd++ )
            {
                // pass global node ids, node index and a pointer to this mesh into the vertex
                mSTKMeshData->mMtkVertices( iVertInd ) =
                        Vertex_Core_STK( this->get_glb_entity_id_from_entity_loc_index( iVertInd, EntityRank::NODE ),
                                iVertInd,
                                this );
            }

            // Setup Vertices interpolation
            mSTKMeshData->mMtkVerticeInterpolation = moris::Cell< Vertex_Interpolation_STK >( tNumNodes );

            for ( moris::moris_index iVertInd = 0; iVertInd < (moris::moris_index)tNumNodes; iVertInd++ )
            {
                // pass global node ids, node index and a pointer to this mesh into the vertex
                mSTKMeshData->mMtkVerticeInterpolation( iVertInd ) = Vertex_Interpolation_STK( &this->get_mtk_vertex( iVertInd ) );

                // set in vertex
                this->get_mtk_vertex_stk( iVertInd ).set_vertex_interpolation( &mSTKMeshData->mMtkVerticeInterpolation( iVertInd ) );
            }

            // Setup Cells
            uint tNumElems = this->get_num_entities( EntityRank::ELEMENT );

            // allocate member data
            mSTKMeshData->mMtkCells = moris::Cell< mtk::Cell_STK >( tNumElems );

            // connectivity factory
            mtk::Cell_Info_Factory tFactory;

            // iterate through buckets
            const stk::mesh::BucketVector& tCellBuckets =
                    mSTKMeshData->mMtkMeshBulkData->get_buckets(
                            stk::topology::ELEMENT_RANK,
                            mSTKMeshData->mMtkMeshMetaData->universal_part() );

            for ( size_t iBucket = 0; iBucket < tCellBuckets.size(); iBucket++ )
            {
                // cell bucket
                stk::mesh::Bucket& tCellBucket = *tCellBuckets[ iBucket ];

                // stk bucket topology
                stk::topology tSTKBucketTopo = tCellBucket.topology();

                // moris topology
                enum CellTopology tBucketTopo = this->stk_topo_to_moris_topo( tSTKBucketTopo );

                // connectivity location
                moris_index tInfoIndex = mSTKMeshData->mCellInfo.size();
                mSTKMeshData->mCellInfo.push_back( tFactory.create_cell_info_sp( tBucketTopo ) );

                for ( size_t iC = 0; iC < tCellBucket.size(); iC++ )
                {
                    // get the stk cell
                    stk::mesh::Entity tSTKCell = tCellBucket[ iC ];

                    // get the id
                    moris_id tId = mSTKMeshData->mMtkMeshBulkData->identifier( tSTKCell );

                    // local index
                    moris_index tIndex = this->get_loc_entity_ind_from_entity_glb_id( tId, EntityRank::ELEMENT );

                    Matrix< IndexMat > tElementToNode =
                            get_entity_connected_to_entity_loc_inds(
                                    tIndex,
                                    EntityRank::ELEMENT,
                                    EntityRank::NODE );

                    // setup vertices of cells
                    moris::Cell< Vertex* > tElementVertices( tElementToNode.numel() );
                    for ( uint iNodes = 0; iNodes < tElementToNode.numel(); iNodes++ )
                    {
                        tElementVertices( iNodes ) = &mSTKMeshData->mMtkVertices( tElementToNode( iNodes ) );
                    }

                    // Add cell to member data
                    mSTKMeshData->mMtkCells( tIndex ) = Cell_STK( mSTKMeshData->mCellInfo( tInfoIndex ),
                            tId,
                            tIndex,
                            tElementVertices,
                            this );
                }
            }
        }
        // ----------------------------------------------------------------------------

        //##############################################
        // Private functions to access mesh information
        //##############################################
        stk::mesh::EntityRank
        Mesh_Core_STK::get_stk_entity_rank( enum EntityRank aMRSEntityRank ) const
        {
            if ( aMRSEntityRank == EntityRank::NODE )
            {
                return stk::topology::NODE_RANK;
            }
            else if ( aMRSEntityRank == EntityRank::EDGE )
            {
                return stk::topology::EDGE_RANK;
            }

            else if ( aMRSEntityRank == EntityRank::FACE )
            {
                return stk::topology::FACE_RANK;
            }
            else if ( aMRSEntityRank == EntityRank::ELEMENT )
            {
                return stk::topology::ELEMENT_RANK;
            }
            else
            {
                return stk::topology::INVALID_RANK;
            }
        }

        // ----------------------------------------------------------------------------

        enum CellTopology
        Mesh_Core_STK::stk_topo_to_moris_topo( stk::topology::topology_t aSTKTopo ) const
        {
            switch ( aSTKTopo )
            {
                case ( stk::topology::TRI_3 ):
                    return CellTopology::TRI3;
                    break;
                case ( stk::topology::TRI_3_2D ):
                    return CellTopology::TRI3;
                    break;
                case ( stk::topology::TRI_6 ):
                    return CellTopology::TRI6;
                    break;
                case ( stk::topology::TRI_6_2D ):
                    return CellTopology::TRI6;
                    break;
                case ( stk::topology::QUAD_4 ):
                    return CellTopology::QUAD4;
                    break;
                case ( stk::topology::QUAD_4_2D ):
                    return CellTopology::QUAD4;
                    break;
                case ( stk::topology::QUAD_8 ):
                    return CellTopology::QUAD8;
                    break;
                case ( stk::topology::QUAD_8_2D ):
                    return CellTopology::QUAD8;
                    break;
                case ( stk::topology::QUAD_9 ):
                    return CellTopology::QUAD9;
                    break;
                case ( stk::topology::QUAD_9_2D ):
                    return CellTopology::QUAD9;
                    break;
                case ( stk::topology::TET_4 ):
                    return CellTopology::TET4;
                    break;
                case ( stk::topology::TET_10 ):
                    return CellTopology::TET10;
                    break;
                case ( stk::topology::HEX_8 ):
                    return CellTopology::HEX8;
                    break;
                case ( stk::topology::HEX_20 ):
                    return CellTopology::HEX20;
                    break;
                case ( stk::topology::HEX_27 ):
                    return CellTopology::HEX27;
                    break;
                case ( stk::topology::WEDGE_6 ):
                    return CellTopology::PRISM6;
                    break;
                default:
                    MORIS_ERROR( 0, "Unhandled stk topology passed in, only ones currently used in MORIS have been added" );
                    return CellTopology::INVALID;
                    break;
            }
        }

        // ----------------------------------------------------------------------------

        Matrix< IdMat >
        Mesh_Core_STK::get_entities_owned_and_shared_by_current_proc(
                EntityRank aEntityRank ) const
        {
            if ( aEntityRank == EntityRank::NODE )
            {
                return mSTKMeshData->mEntityLocaltoGlobalMap( 0 );
            }
            else if ( aEntityRank == EntityRank::EDGE )
            {
                return mSTKMeshData->mEntityLocaltoGlobalMap( 1 );
            }
            else if ( aEntityRank == EntityRank::FACE )
            {
                return mSTKMeshData->mEntityLocaltoGlobalMap( 2 );
            }
            else if ( aEntityRank == EntityRank::ELEMENT )
            {
                return mSTKMeshData->mEntityLocaltoGlobalMap( 3 );
            }
            else
            {
                MORIS_ASSERT( 0, "Invalid rank provided in get_entities_owned_and_shared_by_current_proc." );
            }

            Matrix< IdMat > tDummyConn( 1, 1, 0 );
            return tDummyConn;
        }

        // ----------------------------------------------------------------------------
        // Access entities in selected portion of the mesh
        moris::Matrix< IdMat >
        Mesh_Core_STK::get_entities_in_selector_interface_glob_id(
                EntityRank          aRequestedEntityRank,
                stk::mesh::Selector aSelectedEntities ) const
        {
            // Get selected entities
            std::vector< stk::mesh::Entity > tOutputEntityIDs;
            stk::mesh::get_selected_entities(
                    aSelectedEntities, mSTKMeshData->mMtkMeshBulkData->buckets( this->get_stk_entity_rank( aRequestedEntityRank ) ), tOutputEntityIDs );

            // Interface with STK, get the Ids and return them
            uint            tNumEntity = tOutputEntityIDs.size();
            Matrix< IdMat > tOutputEntityIDMat( tNumEntity, 1, 0 );

            for ( uint i = 0; i < tNumEntity; i++ )
            {
                tOutputEntityIDMat( i, 0 ) = (moris_id)mSTKMeshData->mMtkMeshBulkData->identifier( tOutputEntityIDs[ i ] );
            }

            return tOutputEntityIDMat;
        }

        // ----------------------------------------------------------------------------

        // return processors sharing a particular entity
        moris::Matrix< IdMat >
        Mesh_Core_STK::get_procs_sharing_entity_by_id(
                moris_id        aEntityID,
                enum EntityRank aEntityRank ) const
        {
            // Initialize returning mat with UINT_MAX in case no processors are sharing the given entity
            Matrix< IdMat > tSharedProcsMat( 1, 1, INT_MAX );

            // Check if function is called in a serial run. If so, no parallel information is available.
            uint tProcSize = par_size();
            if ( tProcSize == 1 )
            {
                return tSharedProcsMat;
            }

            // Define entity
            const stk::mesh::Entity tEntity = mSTKMeshData->mMtkMeshBulkData->get_entity( this->get_stk_entity_rank( aEntityRank ), aEntityID );

            // Intialize shared procs
            std::vector< int > tSharedProcsVec;

            // get shared processor IDs
            mSTKMeshData->mMtkMeshBulkData->comm_shared_procs( mSTKMeshData->mMtkMeshBulkData->entity_key( tEntity ), tSharedProcsVec );

            uint tNumSharedProcs = tSharedProcsVec.size();

            // Check if no processors are sharing the given entity
            if ( tNumSharedProcs != 0 )
            {
                tSharedProcsMat.resize( tNumSharedProcs, 1 );
            }

            // Transform from standard to moris call
            for ( uint iProc = 0; iProc < tNumSharedProcs; ++iProc )
            {
                tSharedProcsMat( iProc ) = tSharedProcsVec[ iProc ];
            }

            return tSharedProcsMat;
        }

        // ----------------------------------------------------------------------------

        moris::Cell< moris::Cell< uint > >
        Mesh_Core_STK::get_shared_info_by_entity(
                uint            aNumActiveSharedProcs,
                enum EntityRank aEntityRank )
        {
            // Generate elements shared per processor list
            // -------------------------------------------
            Matrix< IdMat > tEntitiesShared    = this->get_entities_glb_shared_current_proc( aEntityRank );
            uint            tNumEntitiesShared = tEntitiesShared.length();
            moris_id        tParallelRank      = par_rank();

            moris::Cell< moris::Cell< uint > > tTemporaryEntityMapSharingProcs( aNumActiveSharedProcs );

            // Loop over the number of nodes shared to get the shared processors
            for ( uint iElemShared = 0; iElemShared < tNumEntitiesShared; ++iElemShared )
            {
                Matrix< IdMat > tProcsSharing = this->get_procs_sharing_entity_by_id( tEntitiesShared( iElemShared ), aEntityRank );
                for ( uint iProc = 0; iProc < tProcsSharing.length(); ++iProc )
                {
                    if ( tProcsSharing( iProc ) != (moris_id)tParallelRank )
                    {
                        tTemporaryEntityMapSharingProcs( mSTKMeshData->mProcsSharedToIndex[ tProcsSharing( iProc ) ] ).push_back( tEntitiesShared( iElemShared ) );
                    }
                }
            }

            return tTemporaryEntityMapSharingProcs;
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::get_processors_whom_share_entity(
                moris_index      aEntityIndex,
                enum EntityRank  aEntityRank,
                Matrix< IdMat >& aProcsWhomShareEntity ) const
        {
            if ( aEntityRank != EntityRank::ELEMENT )
            {

                // Convert index to ID
                stk::mesh::EntityId tEntityId = { (stk::mesh::EntityId)this->get_glb_entity_id_from_entity_loc_index( aEntityIndex, aEntityRank ) };

                // Get entity
                stk::mesh::Entity tEntity = mSTKMeshData->mMtkMeshBulkData->get_entity( get_stk_entity_rank( aEntityRank ), tEntityId );

                // Intialize shared procs
                std::vector< int > tSharedProcs;

                // get shared processor IDs
                // mSTKMeshData->mMtkMeshBulkData->comm_procs( mSTKMeshData->mMtkMeshBulkData->entity_key( tEntity ), tSharedProcs );

                if ( tSharedProcs.size() == 0 )
                {
                    tSharedProcs.push_back( mSTKMeshData->mMtkMeshBulkData->parallel_owner_rank( tEntity ) );
                }

                // Initialize output
                aProcsWhomShareEntity.resize( 1, tSharedProcs.size() );

                // Cell to vector conversion
                for ( uint i = 0; i < tSharedProcs.size(); i++ )
                {
                    aProcsWhomShareEntity( 0, i ) = tSharedProcs[ i ];
                }
            }
            else
            {
                aProcsWhomShareEntity = mSTKMeshData->mCellSharingData( aEntityIndex ).copy();
            }
        }

        // ----------------------------------------------------------------------------

        uint
        Mesh_Core_STK::get_num_of_entities_shared_with_processor(
                moris_id        aProcessorRank,
                enum EntityRank aEntityRank,
                bool            aSendFlag ) const
        {
            uint tNumEntitiesInCommList = 0;
            if ( aSendFlag )
            {
                tNumEntitiesInCommList = mSTKMeshData->mEntitySendList( (uint)aEntityRank )( aProcessorRank ).n_cols();
            }
            else
            {
                tNumEntitiesInCommList = mSTKMeshData->mEntityReceiveList( (uint)aEntityRank )( aProcessorRank ).n_cols();
            }

            return tNumEntitiesInCommList;
        }

        // ----------------------------------------------------------------------------
        //##############################################
        // Build mesh from data functions internal
        //##############################################
        // Verifications of mesh essential information provided for meshes generated from data.
        void
        Mesh_Core_STK::check_and_update_input_data(
                MtkMeshData& aMeshData )
        {
            // Mesh main variables
            uint tNumNodes = aMeshData.NodeCoords[ 0 ].n_rows();

            // Parallel initialization
            uint tProcSize = par_size();
            // Access mesh data from struc and check input values for indispensable arguments
            MORIS_ASSERT( aMeshData.SpatialDim != nullptr, "Number of spatial dimensions was not provided." );
            MORIS_ASSERT( aMeshData.ElemConn( 0 ) != nullptr, "Element connectivity was not provided." );
            MORIS_ASSERT( aMeshData.NodeCoords != nullptr, "Node coordinates were not provided." );

            // Initialize number of dimensions
            mSTKMeshData->mNumDims = aMeshData.SpatialDim[ 0 ];

            if ( ( aMeshData.NodeProcOwner == NULL ) && ( tProcSize == 1 ) )
            {
                // Do nothing. This information is not needed in serial
            }
            else if ( aMeshData.NodeProcOwner != NULL )
            {
                // Verify sizes
                MORIS_ASSERT( ( aMeshData.NodeProcOwner[ 0 ].numel() == tNumNodes ),
                        "Number of rows for EntProcOwner should match number" );
            }

            // If nodal map was not provided and the simulation is run with 1 processor,
            // just generate the map. The number of ids provided in the map should match
            // the number of coordinates given for each processor.
            if ( aMeshData.LocaltoGlobalNodeMap != nullptr )
            {
                // Verify sizes
                MORIS_ASSERT( aMeshData.LocaltoGlobalNodeMap[ 0 ].numel() == tNumNodes, "Number of rows for LocaltoGlobalNodeMap should match number nodes." );

                mSTKMeshData->mEntityLocaltoGlobalMap( 0 ) = aMeshData.LocaltoGlobalNodeMap[ 0 ];
            }
            else if ( ( aMeshData.LocaltoGlobalNodeMap == nullptr ) && ( tProcSize == 1 ) )
            {
                // Generate nodes maps
                mSTKMeshData->mEntityLocaltoGlobalMap( 0 ) = Matrix< IdMat >( tNumNodes, 1 );
                for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
                {
                    mSTKMeshData->mEntityLocaltoGlobalMap( 0 )( iNode ) = iNode + 1;
                }

                aMeshData.LocaltoGlobalNodeMap = &mSTKMeshData->mEntityLocaltoGlobalMap( 0 );
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
            MORIS_ASSERT( aMeshData.LocaltoGlobalElemMap( 0 ) != nullptr, " No Local to global element map provided" );
            // Verify sizes
            MORIS_ASSERT( aMeshData.size_local_to_global_elem_map() == aMeshData.get_num_elements(),
                    "Number of rows for LocaltoGlobalElemMap should match number of elements provided in connectivity map" );

            if ( aMeshData.LocaltoGlobalElemMap.size() == 0 )
            {
                mSTKMeshData->mEntityLocaltoGlobalMap( 3 ) = ( *aMeshData.LocaltoGlobalElemMap( 0 ) );
            }

            else
            {
                mSTKMeshData->mEntityLocaltoGlobalMap( 3 ) = aMeshData.collapse_element_map();
            }

            // No problem if field information was not provided (i.e., FieldsData, FieldsName, PartNames)
            // Only need to check if data given is consistent.
            if ( aMeshData.FieldsInfo != nullptr )
            {
                this->check_and_update_fields_data( aMeshData );
            }
            mSTKMeshData->mSetNames.resize( 3 );    // Number of ranks for set names
            // No problem if (block,node,side) sets information was not provided
            // Only need to check if data given is consistent.

            if ( aMeshData.SetsInfo != nullptr )
            {
                this->check_and_update_sets_data( aMeshData );
            }
            else
            {
                // Create a block set that contains the entire mesh by default
                mSTKMeshData->mSetNames[ 2 ].resize( 1, "block_1" );
            }
        }

        // ----------------------------------------------------------------------------

        // Verification of fields data arrangement and consistency for meshes generated from data.
        void
        Mesh_Core_STK::check_and_update_fields_data( MtkMeshData& aMeshData )
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
            //        if ( aMeshData.FieldsInfo[0].SetsOwner != nullptr )
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
            //        aMeshData.FieldsInfo[0].FieldsName.resize( mSTKMeshData->mMaxNumFields );
            //        for ( uint iField = tNumFields; iField < mSTKMeshData->mMaxNumFields; ++iField )
            //        {
            //            aMeshData.FieldsInfo[0].FieldsName( iField ) = "dummyField";
            //        }
        }
        // ----------------------------------------------------------------------------

        // Verifications of set data arrangement and consistency for meshes generated from data.
        void
        Mesh_Core_STK::check_and_update_sets_data(
                MtkMeshData& aMeshData )
        {
            ///////////////////////////
            // Checks for block sets //
            ///////////////////////////
            if ( aMeshData.has_mesh_sets() )
            {
                if ( aMeshData.SetsInfo->get_num_block_sets() != 0 )
                {
                    uint tNumBlockSets = aMeshData.SetsInfo->get_num_block_sets();

                    // Communicate with other processors and see which one has the maximum
                    uint tNumGlobalBlockSets         = gather_value_and_bcast_max( tNumBlockSets );
                    mSTKMeshData->mSetRankFlags[ 2 ] = true;
                    mSTKMeshData->mSetNames[ 2 ].resize( tNumGlobalBlockSets );

                    // Loop over the number of block sets
                    for ( uint iBSet = 0; iBSet < tNumBlockSets; ++iBSet )
                    {
                        // Get the block set
                        MtkBlockSetInfo* tBlockSet = aMeshData.SetsInfo->get_block_set( iBSet );

                        // Check if set names were provided
                        if ( !tBlockSet->blockset_has_name() )
                        {
                            mSTKMeshData->mSetNames[ 2 ][ iBSet ] = "BlockSet_" + std::to_string( iBSet );
                            tBlockSet->mBlockSetName              = "BlockSet_" + std::to_string( iBSet );
                        }
                        else
                        {
                            mSTKMeshData->mSetNames[ 2 ][ iBSet ] = tBlockSet->mBlockSetName;
                        }
                    }
                }
                else
                {
                    // Create a block set tha contains the entire mesh by default
                    mSTKMeshData->mSetNames[ 2 ].resize( 1, "block_1" );
                }

                ///////////////////////////
                // Checks for side sets //
                ///////////////////////////
                uint tNumSideSets = aMeshData.SetsInfo->get_num_side_sets();
                if ( tNumSideSets != 0 )
                {
                    mSTKMeshData->mSetRankFlags[ 1 ] = true;

                    // size the mtk mesh member data
                    mSTKMeshData->mSetNames[ 1 ].resize( tNumSideSets );

                    // Loop over the number of block sets
                    for ( uint iSSet = 0; iSSet < tNumSideSets; ++iSSet )
                    {
                        // Get the side set
                        MtkSideSetInfo* tSideSet = aMeshData.SetsInfo->get_side_set( iSSet );

                        // Check if set names were provided
                        if ( !tSideSet->sideset_has_name() )
                        {
                            mSTKMeshData->mSetNames[ 1 ][ iSSet ] = "SideSet_" + std::to_string( iSSet );
                        }
                        else
                        {
                            mSTKMeshData->mSetNames[ 1 ][ iSSet ] = tSideSet->mSideSetName;
                        }

                        // Check if side set specific info was provided
                        MORIS_ASSERT( ( tSideSet->mElemIdsAndSideOrds->n_cols() == 2 ) || ( tSideSet->mElemIdsAndSideOrds->n_cols() == 0 ),
                                "Number of columns in side set should be equal to 2. "
                                "The first column should have element Ids; and the second, side ordinals." );
                    }
                }

                ///////////////////////////
                // Checks for node sets //
                ///////////////////////////
                uint tNumNodeSets = aMeshData.SetsInfo->get_num_node_sets();
                if ( tNumNodeSets != 0 )
                {
                    mSTKMeshData->mSetRankFlags[ 0 ] = true;

                    mSTKMeshData->mSetNames[ 0 ].resize( tNumNodeSets );

                    // Loop over the number of block sets
                    for ( uint iNSet = 0; iNSet < tNumNodeSets; ++iNSet )
                    {
                        // Get node set
                        MtkNodeSetInfo* tNodeSet = aMeshData.SetsInfo->get_node_set( iNSet );

                        // Check if set names were provided
                        if ( !tNodeSet->nodeset_has_name() )
                        {
                            mSTKMeshData->mSetNames[ 0 ][ iNSet ] = "NodeSet_" + std::to_string( iNSet );
                        }
                        else
                        {
                            mSTKMeshData->mSetNames[ 0 ][ iNSet ] = tNodeSet->mNodeSetName;
                        }
                    }
                }
            }
        }

        // Main interface with STK that include calls to functions that provide specific implementation details.
        void
        Mesh_Core_STK::build_mesh(
                MtkMeshData& aMeshData )
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
            stk::mesh::MetaData* meshMeta = new stk::mesh::MetaData( mSTKMeshData->mNumDims );

            // Set member variable as pointer to meta_data
            mSTKMeshData->mMtkMeshMetaData = meshMeta;

            // set aura option
            mAutoAuraOption = aMeshData.AutoAuraOptionInSTK;

            // Declare all additional parts provided by the user (default parts are always created by STK)
            this->declare_mesh_parts( aMeshData );

            // Declare all fields (including coordinates)
            this->declare_mesh_fields( aMeshData );

            // Commit MetaData before populating the BulkData
            mSTKMeshData->mMtkMeshMetaData->commit();

            ////////////////////////////////
            // BULK DATA INITIALIZATION   //
            ////////////////////////////////

            // The BulkData component of a STK Mesh contains entities, entity ownership and ghosting
            // information, connectivity data, and field data. For efficiency, the BulkData API enables access to
            // data via buckets, in addition to via entity and rank.

            // Declare MPI communicator
            stk::ParallelMachine tPM = MPI_COMM_WORLD;

            // Create BulkData Object
            stk::mesh::BulkData* meshBulk = new stk::mesh::BulkData( *mSTKMeshData->mMtkMeshMetaData, tPM, this->get_aura_option() );

            // Set member variable as pointer and bulk_data
            mSTKMeshData->mMtkMeshBulkData = ( meshBulk );

            // Use STK IO to populate a STK Mesh
            MPI_Comm aCommunicator    = MPI_COMM_WORLD;
            mSTKMeshData->mMeshReader = new stk::io::StkMeshIoBroker( aCommunicator );

            // Create mesh database using the IO broker
            mSTKMeshData->mMeshReader->set_bulk_data( *mSTKMeshData->mMtkMeshBulkData );

            // Assign element to node connectivity
            this->populate_mesh_database( aMeshData );

            // Assign coordinates and any other field given by the user
            this->populate_mesh_fields( aMeshData );

            // Generate additional local to global maps (only for meshes generated from data).
            // Elemental and nodal information has been taken care of already in this case.
            // setup maps which include the aura
            mSTKMeshData->mEntityLocaltoGlobalMap = moris::Cell< moris::Matrix< IdMat > >( (uint)EntityRank::END_ENUM, moris::Matrix< IndexMat >( 1, 1, 0 ) );
            mSTKMeshData->mEntitySendList         = moris::Cell< moris::Cell< moris::Matrix< IndexMat > > >( (uint)EntityRank::END_ENUM, moris::Cell< moris::Matrix< IndexMat > >( par_size(), moris::Matrix< IdMat >( 1, 1 ) ) );
            mSTKMeshData->mEntityReceiveList      = moris::Cell< moris::Cell< moris::Matrix< IndexMat > > >( (uint)EntityRank::END_ENUM, moris::Cell< moris::Matrix< IndexMat > >( par_size(), moris::Matrix< IdMat >( 1, 1 ) ) );
            create_communication_lists_and_local_to_global_map( EntityRank::NODE );
            create_communication_lists_and_local_to_global_map( EntityRank::EDGE );
            create_communication_lists_and_local_to_global_map( EntityRank::FACE );
            create_communication_lists_and_local_to_global_map( EntityRank::ELEMENT );

            // Initialize global to local map
            mSTKMeshData->mEntityGlobaltoLocalMap = moris::Cell< std::unordered_map< moris_id, moris_index > >( 4 );
            setup_entity_global_to_local_map( EntityRank::NODE );
            setup_entity_global_to_local_map( EntityRank::ELEMENT );

            if ( aMeshData.CreateAllEdgesAndFaces )
            {
                this->create_additional_communication_lists_from_data();
                setup_entity_global_to_local_map( EntityRank::FACE );
                setup_entity_global_to_local_map( EntityRank::EDGE );
            }

            setup_vertices_and_cell();

            // setup vertex pairing
            setup_parallel_vertex_pairing();

            // setup the parallel cell pairs (including aura)
            setup_parallel_cell_sharing();

            // set timestamp
            mTimeStamp = aMeshData.TimeStamp;
        }

        // ----------------------------------------------------------------------------

        // First declaration to structure the database before filling the data
        void
        Mesh_Core_STK::declare_mesh_parts(
                MtkMeshData& aMeshData )
        {
            // Part is a general term for a subset of the entities in a mesh. STK Mesh automatically creates
            // four parts at startup: the universal part, the locally-owned part, the globally-shared part,
            // and the aura part. In addition, Exodus parts, such as blocks, sidesets, and nodesets, are
            // created if an Exodus file is read in. Each entity in the mesh must be a member of one or more parts.

            uint tNumElementTypes = aMeshData.ElemConn.size();

            // Loop over the different element types and declare a part
            for ( uint iET = 0; iET < tNumElementTypes; iET++ )
            {
                uint tNumNodesPerElem = aMeshData.ElemConn( iET )->n_cols();
                // Declare and initialize topology type. Also check if element type is supported

                stk::topology::topology_t tTopology;
                if ( aMeshData.CellTopology( iET ) == CellTopology::INVALID )
                {
                    tTopology = get_mesh_topology( mSTKMeshData->mNumDims, tNumNodesPerElem );
                }
                else
                {
                    tTopology = get_stk_topo( aMeshData.CellTopology( iET ) );
                }

                // Add default part if no block sets were provided
                stk::mesh::Part& tBlock = mSTKMeshData->mMtkMeshMetaData->declare_part_with_topology( "noblock_" + std::to_string( iET ), tTopology );

                if ( aMeshData.MarkNoBlockForIO )
                {
                    // Add Part to the IOBroker (needed for output).
                    stk::io::put_io_part_attribute( tBlock );
                }
            }

            if ( aMeshData.SetsInfo != nullptr )    // For all (block, node, side) sets
            {
                ////////////////////////
                // Declare block sets //
                ////////////////////////
                uint tNumBlockSets = aMeshData.SetsInfo->get_num_block_sets();
                for ( uint iSet = 0; iSet < tNumBlockSets; ++iSet )
                {

                    MtkBlockSetInfo* tBlockSet = aMeshData.SetsInfo->get_block_set( iSet );
                    // Declare part and add it to the IOBroker (needed for output).

                    stk::mesh::Part& aSetPart = mSTKMeshData->mMtkMeshMetaData->declare_part_with_topology( tBlockSet->mBlockSetName,
                            get_stk_topo( tBlockSet->mBlockSetTopo ),
                            true );

                    // mark as parallel consistent or not (needed for XTK aura to allow for inconsistent aura declaration)
                    aSetPart.entity_membership_is_parallel_consistent( tBlockSet->mParallelConsistencyReq );

                    // Add Part to the IOBroker (needed for output)
                    stk::io::put_io_part_attribute( aSetPart );
                }

                ///////////////////////
                // Declare side sets //
                ///////////////////////
                if ( mSTKMeshData->mSetRankFlags[ 1 ] )
                {
                    uint tNumSideSets = aMeshData.SetsInfo->get_num_side_sets();

                    for ( uint iSet = 0; iSet < tNumSideSets; ++iSet )
                    {
                        // get side set
                        moris::mtk::MtkSideSetInfo* tSideSet = aMeshData.SetsInfo->SideSetsInfo( iSet );

                        // Declare part and add it to the IOBroker (needed for output).
                        stk::mesh::Part& aSetPart = mSTKMeshData->mMtkMeshMetaData->declare_part( mSTKMeshData->mSetNames[ 1 ][ iSet ], mSTKMeshData->mMtkMeshMetaData->side_rank(), true );

                        // mark as parallel consistent or not (needed for XTK aura to allow for inconsistent aura declaration)
                        aSetPart.entity_membership_is_parallel_consistent( tSideSet->mParallelConsistencyReq );

                        // Add Part to the IOBroker (needed for output).
                        stk::io::put_io_part_attribute( aSetPart );
                    }
                }

                ///////////////////////
                // Declare node sets //
                ///////////////////////
                if ( mSTKMeshData->mSetRankFlags[ 0 ] )
                {
                    uint tNumNodeSets = aMeshData.SetsInfo->get_num_node_sets();

                    for ( uint iSet = 0; iSet < tNumNodeSets; ++iSet )
                    {
                        // Declare part and add it to the IOBroker (needed for output).
                        stk::mesh::Part& aSetPart = mSTKMeshData->mMtkMeshMetaData->declare_part( mSTKMeshData->mSetNames[ 0 ][ iSet ], stk::topology::NODE_RANK, true );

                        // Add Part to the IOBroker (needed for output).
                        stk::io::put_io_part_attribute( aSetPart );
                    }
                }
            }
        }

        // ----------------------------------------------------------------------------

        // Second declaration to structure the database before filling the data

        void
        Mesh_Core_STK::declare_mesh_fields(
                MtkMeshData& aMeshData )
        {
            // Fields are data associated with mesh entities. Examples include coordinates, velocity,
            // displacement, and temperature. A field in STK Mesh can hold any data type (e.g., double or int)
            // and any number of scalars per entity (e.g., nodal velocity field has three doubles per node).
            // A field can be allocated (defined) on a whole mesh or on only a subset (part) of that mesh.
            // For example, a material property can be allocated on a specified element block.

            // Declare coordinates field (only if not supplementary)
            // If supplementary we assume this is done when loading
            // the exodus file.
            if ( !aMeshData.SupplementaryToFile )
            {
                if ( aMeshData.SpatialDim[ 0 ] == 2 )
                {
                    Field2CompReal* tCoord_field = &mSTKMeshData->mMtkMeshMetaData->declare_field< Field2CompReal >( stk::topology::NODE_RANK, "coordinates" );
                    stk::mesh::put_field_on_mesh( *tCoord_field, mSTKMeshData->mMtkMeshMetaData->universal_part(), (stk::mesh::FieldTraits< stk::mesh::Field< real > >::data_type*)nullptr );
                }
                else
                {
                    Field3CompReal* tCoord_field = &mSTKMeshData->mMtkMeshMetaData->declare_field< Field3CompReal >( stk::topology::NODE_RANK, "coordinates" );
                    stk::mesh::put_field_on_mesh( *tCoord_field, mSTKMeshData->mMtkMeshMetaData->universal_part(), (stk::mesh::FieldTraits< stk::mesh::Field< real > >::data_type*)nullptr );
                }
            }

            // Declare all additional fields provided by the user
            if ( aMeshData.FieldsInfo != nullptr )
            {
                // Iterate over real scalar fields amd declare them
                uint tNumRealScalarFields = aMeshData.FieldsInfo->get_num_real_scalar_fields();

                for ( uint iF = 0; iF < tNumRealScalarFields; iF++ )
                {
                    Scalar_Field_Info< DDRMat >* tRealScalarField = ( aMeshData.FieldsInfo->mRealScalarFields )( iF );
                    enum EntityRank              tFieldEntityRank = tRealScalarField->get_field_entity_rank();
                    std::string                  tFieldName       = tRealScalarField->get_field_name();
                    stk::mesh::Selector          tFieldPart;

                    if ( !tRealScalarField->field_has_part_name() )
                    {
                        tFieldPart = mSTKMeshData->mMtkMeshMetaData->universal_part();
                    }
                    else
                    {
                        tFieldPart = *mSTKMeshData->mMtkMeshMetaData->get_part( tRealScalarField->get_part_name() );
                    }

                    stk::mesh::Field< real >& tSTKRealScalarField = mSTKMeshData->mMtkMeshMetaData->declare_field< stk::mesh::Field< real > >( this->get_stk_entity_rank( tFieldEntityRank ), tFieldName, 1 );

                    stk::mesh::put_field_on_mesh( tSTKRealScalarField, tFieldPart, (stk::mesh::FieldTraits< stk::mesh::Field< real > >::data_type*)nullptr );

                    stk::io::set_field_role( tSTKRealScalarField, Ioss::Field::TRANSIENT );
                }

                // iterate over real matrix fields and declare them

                // Iterate over real scalar fields amd declare them
                uint tNumRealMatrixFields = aMeshData.FieldsInfo->get_num_real_matrix_fields();
                for ( uint iF = 0; iF < tNumRealMatrixFields; iF++ )
                {
                    Matrix_Field_Info< DDRMat >* tRealMatrixField = ( aMeshData.FieldsInfo->mRealMatrixFields )( iF );
                    enum EntityRank              tFieldEntityRank = tRealMatrixField->get_field_entity_rank();
                    std::string                  tFieldName       = tRealMatrixField->get_field_name();
                    const uint                   tNumRows         = tRealMatrixField->get_num_rows();
                    const uint                   tNumCols         = tRealMatrixField->get_num_cols();

                    stk::mesh::Selector tFieldPart;
                    if ( !tRealMatrixField->field_has_part_name() )
                    {
                        tFieldPart = mSTKMeshData->mMtkMeshMetaData->universal_part();
                    }
                    else
                    {
                        tFieldPart = *mSTKMeshData->mMtkMeshMetaData->get_part( tRealMatrixField->get_part_name() );
                    }

                    internal_declare_mesh_real_matrix_fields( tFieldName, tFieldEntityRank, tFieldPart, tNumRows, tNumCols );
                }
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::add_supplementary_fields_to_declare_at_output( MtkMeshData& aMeshData )
        {
            MORIS_ASSERT( aMeshData.SupplementaryToFile,
                    "add_supplementary_fields_to_declare_at_output only should be called during load mesh from file construction" );

            uint tNumRealScalarFields = aMeshData.FieldsInfo->get_num_real_scalar_fields();

            mSTKMeshData->mRealNodeScalarFieldsToAddToOutput = moris::Cell< Field1CompReal* >( tNumRealScalarFields );

            for ( uint iF = 0; iF < tNumRealScalarFields; iF++ )
            {
                Scalar_Field_Info< DDRMat >* tRealScalarField = ( aMeshData.FieldsInfo->mRealScalarFields )( iF );

                enum EntityRank tFieldEntityRank = tRealScalarField->get_field_entity_rank();

                std::string tFieldName = tRealScalarField->get_field_name();

                stk::mesh::EntityRank     tEntityRank = this->get_stk_entity_rank( tFieldEntityRank );
                stk::mesh::Field< real >* tField      = mSTKMeshData->mMtkMeshMetaData->get_field< stk::mesh::Field< real > >( tEntityRank, tFieldName );

                mSTKMeshData->mRealNodeScalarFieldsToAddToOutput( iF ) = tField;
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::internal_declare_mesh_real_matrix_fields(
                std::string         aFieldName,
                enum EntityRank     aFieldRank,
                stk::mesh::Selector aFieldPart,
                uint                aNumRows,
                uint                aNumCols )
        {
            uint tNumFieldComp = aNumCols * aNumRows;

            switch ( tNumFieldComp )
            {
                case 1:    // Scalar Field
                {
                    // Declare fields
                    mSTKMeshData->mField1CompVecsReal.push_back( &mSTKMeshData->mMtkMeshMetaData->declare_field< Field1CompReal >( this->get_stk_entity_rank( aFieldRank ), aFieldName ) );
                    stk::mesh::put_field_on_mesh( *mSTKMeshData->mField1CompVecsReal.back(), aFieldPart, (stk::mesh::FieldTraits< Field1CompReal >::data_type*)nullptr );
                    stk::io::set_field_role( *mSTKMeshData->mField1CompVecsReal.back(), Ioss::Field::TRANSIENT );
                    break;
                }
                case 2:    // Vector Field with 2 components
                {
                    // Declare fields
                    mSTKMeshData->mField2CompVecsReal.push_back( &mSTKMeshData->mMtkMeshMetaData->declare_field< Field2CompReal >( this->get_stk_entity_rank( aFieldRank ), aFieldName ) );
                    stk::mesh::put_field_on_mesh( *mSTKMeshData->mField2CompVecsReal.back(), aFieldPart, (stk::mesh::FieldTraits< Field2CompReal >::data_type*)nullptr );
                    stk::io::set_field_role( *mSTKMeshData->mField2CompVecsReal.back(), Ioss::Field::TRANSIENT );
                    break;
                }
                case 3:    // Vector Field with 3 components
                {
                    // Declare fields
                    mSTKMeshData->mField3CompVecsReal.push_back( &mSTKMeshData->mMtkMeshMetaData->declare_field< Field3CompReal >( this->get_stk_entity_rank( aFieldRank ), aFieldName ) );
                    stk::mesh::put_field_on_mesh( *mSTKMeshData->mField3CompVecsReal.back(), aFieldPart, (stk::mesh::FieldTraits< Field3CompReal >::data_type*)nullptr );
                    stk::io::set_field_role( *mSTKMeshData->mField3CompVecsReal.back(), Ioss::Field::TRANSIENT );
                    break;
                }
                case 4:    // Tensor Field with 4 components
                {
                    // Declare fields
                    mSTKMeshData->mField4CompVecsReal.push_back( &mSTKMeshData->mMtkMeshMetaData->declare_field< Field4CompReal >( this->get_stk_entity_rank( aFieldRank ), aFieldName ) );
                    stk::mesh::put_field_on_mesh( *mSTKMeshData->mField4CompVecsReal.back(), aFieldPart, (stk::mesh::FieldTraits< Field4CompReal >::data_type*)nullptr );
                    stk::io::set_field_role( *mSTKMeshData->mField4CompVecsReal.back(), Ioss::Field::TRANSIENT );
                    break;
                }
                case 9:    // Tensor Field with 9 components
                {
                    // Declare fields
                    mSTKMeshData->mField9CompVecsReal.push_back( &mSTKMeshData->mMtkMeshMetaData->declare_field< Field9CompReal >( this->get_stk_entity_rank( aFieldRank ), aFieldName ) );
                    stk::mesh::put_field_on_mesh( *mSTKMeshData->mField9CompVecsReal.back(), aFieldPart, (stk::mesh::FieldTraits< Field9CompReal >::data_type*)nullptr );
                    stk::io::set_field_role( *mSTKMeshData->mField9CompVecsReal.back(), Ioss::Field::TRANSIENT );
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

        // Declare size of a field (per entity) and throw an error if it is not supported
        void
        Mesh_Core_STK::internal_declare_mesh_field(
                MtkMeshData& aMeshData,
                uint         iField )
        {
            //        // Get field variables
            //        uint tNumFieldComp     = aMeshData.FieldsInfo[0].FieldsData[0]( iField ).n_cols();
            //        std::string tFieldName = aMeshData.FieldsInfo[0].FieldsName( iField );
            //        EntityRank tFieldRank  = aMeshData.FieldsInfo[0].FieldsRank( iField );
            //
            //        stk::mesh::EntityRank tStkFieldRank = this->get_stk_entity_rank( tFieldRank );
            //        stk::mesh::Selector aFieldPart      = mSTKMeshData->mMtkMeshMetaData->universal_part();
            //
            //        if ( aMeshData.FieldsInfo[0].SetsOwner != nullptr )
            //        {
            //            if ( !aMeshData.FieldsInfo[0].SetsOwner[0]( iField ).empty() )
            //            {
            //                aFieldPart = *mSTKMeshData->mMtkMeshMetaData->get_part( aMeshData.FieldsInfo[0].SetsOwner[0]( iField ) );
            //            }
            //        }
            //
            //        switch ( tNumFieldComp )
            //        {
            //        case 1: // Scalar Field
            //        {
            //            // Declare fields
            //            mField1CompVec.push_back( & mSTKMeshData->mMtkMeshMetaData->declare_field< Field1Comp >( tStkFieldRank, tFieldName ) );
            //            stk::mesh::put_field( *mField1CompVec.back(), aFieldPart, 1 );
            //            break;
            //        }
            //        case 2: // Vector Field with 2 components
            //        {
            //            // Declare fields
            //            mField2CompVec.push_back( & mSTKMeshData->mMtkMeshMetaData->declare_field< Field2Comp >( tStkFieldRank, tFieldName ) );
            //            stk::mesh::put_field( *mField2CompVec.back(), aFieldPart );
            //            break;
            //        }
            //        case 3: // Vector Field with 3 components
            //        {
            //            // Declare fields
            //            mField3CompVec.push_back( & mSTKMeshData->mMtkMeshMetaData->declare_field< Field3Comp >( tStkFieldRank, tFieldName ) );
            //            stk::mesh::put_field( *mField3CompVec.back(), aFieldPart );
            //            break;
            //        }
            //        case 4: // Tensor Field with 4 components
            //        {
            //            // Declare fields
            //            mField4CompVec.push_back( & mSTKMeshData->mMtkMeshMetaData->declare_field< Field4Comp >( tStkFieldRank, tFieldName ) );
            //            stk::mesh::put_field( *mField4CompVec.back(), aFieldPart );
            //            break;
            //        }
            //        case 9: // Tensor Field with 9 components
            //        {
            //            // Declare fields
            //            mField9CompVec.push_back( & mSTKMeshData->mMtkMeshMetaData->declare_field< Field9Comp >( tStkFieldRank, tFieldName ) );
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
        Mesh_Core_STK::populate_mesh_database(
                MtkMeshData& aMeshData )
        {
            ///////////////////////////////
            // Begin modification cycle  //
            ///////////////////////////////
            mSTKMeshData->mMtkMeshBulkData->modification_begin();

            // Setup global to local element map
            setup_cell_global_to_local_map( aMeshData );

            // Setup global to local node map
            setup_vertex_global_to_local_map( aMeshData );

            // setup nodes
            //        this->process_nodes(aMeshData);

            // Generate basic mesh information
            this->process_block_sets( aMeshData );

            // declare nodes

            // Declare node sets to mesh if they exist
            if ( mSTKMeshData->mSetRankFlags[ 0 ] )
            {
                this->process_node_sets( aMeshData );
            }

            // Add node sharing
            this->process_node_sharing_information( aMeshData );

            ///////////////////////////////
            // Close modification cycle  //
            ///////////////////////////////
            mSTKMeshData->mMtkMeshBulkData->modification_end();

            if ( aMeshData.CreateAllEdgesAndFaces )
            {    // If the user requires create all additional entities, use the functions below.
                // Note that this could potentially increase significantly memory usage and time for
                // creating the mesh for big amounts of data.

                // Create mesh edge entities
                stk::mesh::create_edges( *mSTKMeshData->mMtkMeshBulkData );
                mSTKMeshData->mCreatedEdges = true;

                // Create mesh face entities
                stk::mesh::create_faces( *mSTKMeshData->mMtkMeshBulkData, true );
                mSTKMeshData->mCreatedFaces = true;
            }

            // Declare node sets to mesh
            if ( mSTKMeshData->mSetRankFlags[ 1 ] )
            {
                // If side sets were provided, generate only edges and/or faces of the
                // corresponding sets of the elements containing such entities. The elements
                // of the side sets entities will be moved to another part.
                this->process_side_sets( aMeshData );
            }
        }

        //------------------------------------------------------------------------------

        void
        Mesh_Core_STK::setup_cell_global_to_local_map(
                MtkMeshData& aMeshData )
        {
            uint tNumElemTypes = aMeshData.LocaltoGlobalElemMap.size();

            // Loop over element types
            moris_index tCount = 0;
            for ( uint iET = 0; iET < tNumElemTypes; iET++ )
            {
                uint tNumElements = aMeshData.LocaltoGlobalElemMap( iET )->numel();

                for ( uint iElem = 0; iElem < tNumElements; iElem++ )
                {
                    uint tElemId = ( *aMeshData.LocaltoGlobalElemMap( iET ) )( iElem );
                    if ( mSTKMeshData->mEntityGlobaltoLocalMap( 3 ).find( tElemId ) == mSTKMeshData->mEntityGlobaltoLocalMap( 3 ).end() )
                    {
                        mSTKMeshData->mEntityGlobaltoLocalMap( 3 )[ tElemId ] = tCount;
                        tCount++;
                    }
                    else
                    {
                        MORIS_ERROR( 0, "Element Id already in map, does your LocaltoGlobalElemMap have the same id twice?" );
                    }
                }
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::setup_vertex_global_to_local_map(
                MtkMeshData& aMeshData )
        {
            uint tNumNodes = aMeshData.LocaltoGlobalNodeMap->numel();

            // Loop over nodes
            moris_index tCount = 0;
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                uint tNodeId = ( *aMeshData.LocaltoGlobalNodeMap )( iNode );

                if ( mSTKMeshData->mEntityGlobaltoLocalMap( 0 ).find( tNodeId ) == mSTKMeshData->mEntityGlobaltoLocalMap( 0 ).end() )
                {
                    mSTKMeshData->mEntityGlobaltoLocalMap( 0 )[ tNodeId ] = tCount;
                    tCount++;
                }
                else
                {
                    MORIS_ERROR( 0, "Node Id already in map, does your LocaltoGlobalNodeMap have the same id twice?" );
                }
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::setup_entity_global_to_local_map( enum EntityRank aEntityRank )
        {
            moris::Matrix< IdMat > const & tEntityIds = mSTKMeshData->mEntityLocaltoGlobalMap( (uint)aEntityRank );

            uint     tNumEntities = tEntityIds.numel();
            moris_id tCount       = 0;
            for ( uint i = 0; i < tNumEntities; i++ )
            {
                if ( mSTKMeshData->mEntityGlobaltoLocalMap( (uint)aEntityRank ).find( tEntityIds( i ) ) == mSTKMeshData->mEntityGlobaltoLocalMap( (uint)aEntityRank ).end() )
                {
                    mSTKMeshData->mEntityGlobaltoLocalMap( (uint)aEntityRank )[ tEntityIds( i ) ] = tCount;
                    tCount++;
                }
                else
                {
                    MORIS_ERROR( 0, "Node Id already in map, does your LocaltoGlobalNodeMap have the same id twice?" );
                }
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::process_nodes( MtkMeshData& aMeshData )
        {
            // number of nodes
            moris::uint tNumNodes = aMeshData.get_num_nodes();

            for ( moris::uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                moris_id          tNodeGlobalId = ( *aMeshData.LocaltoGlobalNodeMap )( iNode );
                stk::mesh::Entity tNode         = mSTKMeshData->mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, tNodeGlobalId );
                tNode                           = mSTKMeshData->mMtkMeshBulkData->declare_entity( stk::topology::NODE_RANK, tNodeGlobalId, mSTKMeshData->mMtkMeshMetaData->universal_part() );
            }
        }

        // ----------------------------------------------------------------------------
        // Add all blocks information to database
        void
        Mesh_Core_STK::process_block_sets(
                MtkMeshData& aMeshData )
        {

            // Get all blocks an element belongs to
            Matrix< IndexMat > tElementToBlock = process_cell_block_membership( aMeshData );

            uint tNumBlocks = 0;
            if ( aMeshData.has_mesh_sets() )
            {
                tNumBlocks = aMeshData.SetsInfo->get_num_block_sets();
            }

            // Get part vector ordered by the order found in aMeshData block sets
            stk::mesh::PartVector tBlockParts = get_block_set_part_vector( aMeshData );

            // iterate over element types
            for ( uint iET = 0; iET < aMeshData.ElemConn.size(); iET++ )
            {
                // iterate over elements of this type
                for ( uint iElem = 0; iElem < aMeshData.LocaltoGlobalElemMap( iET )->numel(); iElem++ )
                {
                    moris_id    tCellId    = ( *aMeshData.LocaltoGlobalElemMap( iET ) )( iElem );
                    moris_index tCellIndex = get_loc_entity_ind_from_entity_glb_id( tCellId, EntityRank::ELEMENT );

                    // get the part vector associated with this cell
                    stk::mesh::PartVector tCellParts;
                    for ( uint iBl = 0; iBl < tNumBlocks; iBl++ )
                    {
                        if ( tElementToBlock( tCellIndex, iBl ) != std::numeric_limits< moris_index >::max() )
                        {
                            uint             tBlockIndex = tElementToBlock( tCellIndex, iBl );
                            stk::mesh::Part* tPart       = tBlockParts[ tBlockIndex ];
                            tCellParts.push_back( tPart );
                        }
                        else
                        {
                            break;
                        }
                    }

                    if ( tCellParts.size() == 0 )
                    {
                        stk::mesh::Part* tBlock = mSTKMeshData->mMtkMeshMetaData->get_part( "noblock_" + std::to_string( iET ) );
                        tCellParts.push_back( tBlock );
                    }

                    // Declare element
                    stk::mesh::Entity tElement = mSTKMeshData->mMtkMeshBulkData->declare_entity( stk::topology::ELEM_RANK, tCellId, tCellParts );

                    for ( uint node_i = 0; node_i < aMeshData.ElemConn( iET )->n_cols(); ++node_i )
                    {
                        moris_id tNodeGlobalId = ( stk::mesh::EntityId )( *aMeshData.ElemConn( iET ) )( iElem, node_i );

                        if ( tNodeGlobalId != 0 )
                        {
                            stk::mesh::Entity tNode = mSTKMeshData->mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, tNodeGlobalId );

                            if ( !mSTKMeshData->mMtkMeshBulkData->is_valid( tNode ) )
                            {
                                tNode = mSTKMeshData->mMtkMeshBulkData->declare_entity( stk::topology::NODE_RANK, tNodeGlobalId, mSTKMeshData->mMtkMeshMetaData->universal_part() );
                            }

                            mSTKMeshData->mMtkMeshBulkData->declare_relation( tElement, tNode, node_i );
                        }
                    }
                }
            }
        }

        // ----------------------------------------------------------------------------

        Matrix< IndexMat >
        Mesh_Core_STK::process_cell_block_membership(
                MtkMeshData& aMeshData )
        {

            uint               tNumBlocks = 0;
            Matrix< IndexMat > tElementToBlock;
            if ( aMeshData.has_mesh_sets() )
            {
                tNumBlocks      = aMeshData.SetsInfo->get_num_block_sets();
                tElementToBlock = Matrix< IndexMat >( aMeshData.get_num_elements(),
                        tNumBlocks,
                        std::numeric_limits< moris_index >::max() );
            }

            else
            {
                tElementToBlock = Matrix< IndexMat >( aMeshData.get_num_elements(),
                        1,
                        std::numeric_limits< moris_index >::max() );
            }

            // counter
            Matrix< IndexMat > tElementToBlockCounter( aMeshData.get_num_elements(), 1, 0 );

            // Iterate through blocks
            for ( uint iBlock = 0; iBlock < tNumBlocks; iBlock++ )
            {
                // get pointer to block set
                MtkBlockSetInfo* tBlockSet = aMeshData.SetsInfo->get_block_set( iBlock );

                // Iterate through elements in block
                uint tNumElemsInBlock = tBlockSet->mCellIdsInSet->numel();
                for ( uint iElem = 0; iElem < tNumElemsInBlock; iElem++ )
                {
                    // Cell index from cell id
                    //                stk::mesh::Entity aElemEntity       = mSTKMeshData->mMtkMeshBulkData->get_entity( stk::topology::ELEMENT_RANK, (*tBlockSet->mCellIdsInSet)(iElem) );

                    moris_index tElemIndex = get_loc_entity_ind_from_entity_glb_id( ( *tBlockSet->mCellIdsInSet )( iElem ), EntityRank::ELEMENT );

                    // number of blocks this element belongs to
                    uint tElemCount = tElementToBlockCounter( tElemIndex );

                    // Add block index to element to block matrix
                    tElementToBlock( tElemIndex, tElemCount ) = iBlock;

                    tElementToBlockCounter( tElemIndex )++;
                }
            }

            return tElementToBlock;
        }

        // ----------------------------------------------------------------------------

        stk::mesh::PartVector
        Mesh_Core_STK::get_block_set_part_vector( MtkMeshData& aMeshData )
        {
            if ( aMeshData.has_mesh_sets() )
            {
                uint                  tNumBlocks = aMeshData.SetsInfo->get_num_block_sets();
                stk::mesh::PartVector tBlockSetParts( tNumBlocks );

                // Iterate over block sets and get block part from the block set name
                for ( uint iBlock = 0; iBlock < tNumBlocks; iBlock++ )
                {
                    std::string      tBlockName = aMeshData.SetsInfo->get_block_set( iBlock )->mBlockSetName;
                    stk::mesh::Part* tBlock     = mSTKMeshData->mMtkMeshMetaData->get_part( tBlockName );
                    tBlockSetParts[ iBlock ]    = { tBlock };
                }

                return tBlockSetParts;
            }

            else
            {
                return stk::mesh::PartVector( 0 );
            }
        }

        // ----------------------------------------------------------------------------

        void
        Mesh_Core_STK::process_node_sharing_information( MtkMeshData& aMeshData )
        {
            if ( aMeshData.has_node_sharing_info() )
            {
                //           moris_id tParRank = par_rank();
                uint        tNumNodes          = aMeshData.get_num_nodes();
                uint        tMaxNumProcsShared = aMeshData.NodeProcsShared->n_cols();
                moris_index tMyRank            = par_rank();

                for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                {
                    stk::mesh::Entity aNode = mSTKMeshData->mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, ( *aMeshData.LocaltoGlobalNodeMap )( iNode ) );

                    if ( mSTKMeshData->mMtkMeshBulkData->is_valid( aNode ) )
                    {
                        for ( uint iShare = 0; iShare < tMaxNumProcsShared; iShare++ )
                        {
                            moris_id tShareProcRank = ( *aMeshData.NodeProcsShared )( iNode, iShare );

                            if ( tShareProcRank != MORIS_ID_MAX && tShareProcRank != tMyRank )
                            {
                                mSTKMeshData->mMtkMeshBulkData->add_node_sharing( aNode, tShareProcRank );
                            }
                        }
                    }
                }
            }
        }

        // ----------------------------------------------------------------------------
        // Add all fields information to database
        void
        Mesh_Core_STK::populate_mesh_fields(
                MtkMeshData& aMeshData )
        {
            // Get the coordinates field from Stk
            stk::mesh::FieldBase const * aCoord_field_i = mSTKMeshData->mMtkMeshMetaData->coordinate_field();
            uint                         tNumNodes      = aMeshData.NodeCoords->n_rows();

            // Loop over the number of nodes
            for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
            {
                // Get global Id of current node and create "node entity" for stk mesh
                uint              aId   = ( *aMeshData.LocaltoGlobalNodeMap )( iNode );
                stk::mesh::Entity aNode = mSTKMeshData->mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aId );

                // Store the coordinates of the current node
                if ( mSTKMeshData->mMtkMeshBulkData->is_valid( aNode ) )
                {
                    // Add coordinates information to the BulkData
                    double* tCoord_data = static_cast< double* >( stk::mesh::field_data( *aCoord_field_i, aNode ) );
                    for ( uint iDim = 0; iDim < mSTKMeshData->mNumDims; ++iDim )
                    {
                        tCoord_data[ iDim ] = ( *aMeshData.NodeCoords )( iNode, iDim );
                    }
                }
            }

            if ( aMeshData.FieldsInfo != nullptr )
            {
                // Get the number of fields
                uint tNumRealScalarFields = aMeshData.FieldsInfo->get_num_real_scalar_fields();
                for ( uint iF = 0; iF < tNumRealScalarFields; iF++ )
                {
                    Scalar_Field_Info< DDRMat >* tRealScalarField = ( aMeshData.FieldsInfo->mRealScalarFields )( iF );

                    if ( tRealScalarField->field_has_data() )
                    {
                        populate_field_data_scalar_field( tRealScalarField );
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
                for ( uint iF = 0; iF < tNumRealMatrixFields; iF++ )
                {
                    Matrix_Field_Info< DDRMat >* tRealMatrixField = ( aMeshData.FieldsInfo->mRealMatrixFields )( iF );
                    if ( tRealMatrixField->field_has_data() )
                    {
                        populate_field_data_matrix_field( tRealMatrixField );
                    }
                }
            }
        }

        // ----------------------------------------------------------------------------

        // Provide element type (Hex8, Tri3, etc) and throw error if element is not supported yet.
        stk::topology::topology_t
        Mesh_Core_STK::get_mesh_topology(
                uint aModelDim,
                uint aNumNodesInElem )
        {
            stk::topology::topology_t tTopology = stk::topology::INVALID_TOPOLOGY;

            // MTK supports the following 1D, 2D and 3D element topology temporarily

            if ( aModelDim == 1 )    // 1D
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
                        MORIS_ASSERT( 0, "MTK mesh build from data currently handles only LINE_2 for 1D elements." );
                        break;
                    }
                }
            }
            else if ( aModelDim == 2 )    // 2D
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
                        MORIS_ASSERT( 0, "MTK mesh build from data currently handles only TRI_3, TRI_6_2D, QUAD_4 and QUAD_9_2D for 2D elements." );
                        break;
                    }
                }
            }
            else if ( aModelDim == 3 )    // 3D
            {
                switch ( aNumNodesInElem )
                {
                    case 3:
                    {
                        tTopology = stk::topology::TRI_3;
                        break;
                    }
                    case 4:
                    {
                        tTopology = stk::topology::TET_4;
                        break;
                    }
                    case 6:
                        tTopology = stk::topology::WEDGE_6;
                        break;
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
                        MORIS_ASSERT( 0, "MTK mesh build from data currently handles only TET_4, WEDGE_6, HEX8, HEX_20 and HEX_27 for 3D elements." );
                        break;
                    }
                }
            }

            return tTopology;
        }

        // ----------------------------------------------------------------------------

        //    stk::topology::topology_t
        //    Mesh_Core_STK::get_stk_entity_rank(enum EntityRank aMTKEntityRank)
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

        // ----------------------------------------------------------------------------

        // Provide element type (Hex8, Tri3, etc) and throw error if element is not supported yet.
        stk::topology::topology_t
        Mesh_Core_STK::get_stk_topo( enum CellTopology aMTKCellTopo )
        {
            stk::topology::topology_t tTopology = stk::topology::INVALID_TOPOLOGY;

            moris::uint tSpatialDim = this->get_spatial_dim();

            // MTK supports the following 1D, 2D and 3D element topology temporarily
            switch ( aMTKCellTopo )
            {
                case CellTopology::TRI3:
                    if ( tSpatialDim == 2 )
                    {
                        tTopology = stk::topology::TRI_3_2D;
                    }
                    else if ( tSpatialDim == 3 )
                    {
                        tTopology = stk::topology::TRI_3;
                    }
                    break;
                case CellTopology::QUAD4:
                    if ( tSpatialDim == 2 )
                    {
                        tTopology = stk::topology::QUAD_4_2D;
                    }
                    else if ( tSpatialDim == 3 )
                    {
                        tTopology = stk::topology::QUAD_4;
                    }
                    break;
                case CellTopology::QUAD9:
                    if ( tSpatialDim == 2 )
                    {
                        tTopology = stk::topology::QUAD_9_2D;
                    }
                    else if ( tSpatialDim == 3 )
                    {
                        tTopology = stk::topology::QUAD_9;
                    }
                    break;
                case CellTopology::TET4:
                    tTopology = stk::topology::TET_4;
                    break;
                case CellTopology::TET10:
                    tTopology = stk::topology::TET_10;
                    break;
                case CellTopology::HEX8:
                    tTopology = stk::topology::HEX_8;
                    break;
                case CellTopology::HEX27:
                    tTopology = stk::topology::HEX_27;
                    break;
                case CellTopology::PRISM6:
                    tTopology = stk::topology::WEDGE_6;
                    break;
                default:
                    MORIS_ERROR( 0, "Unsupported element type in get_mesh_topology" );
                    break;
            }

            return tTopology;
        }

        // ----------------------------------------------------------------------------

        // Function to create edges and faces communication lists in parallel for meshes generated from data
        void
        Mesh_Core_STK::create_additional_communication_lists_from_data()
        {
            this->create_facets_communication_lists();
        }
        // ----------------------------------------------------------------------------

        // Function to create edges and faces communication lists in parallel
        void
        Mesh_Core_STK::create_facets_communication_lists()
        {
            // Create maps for edges if the problem is not 1D, and maps for faces if it is 3D.
            // NOTE: Not supporting 1D elements in 2D nor 3D space.
            if ( mSTKMeshData->mNumDims > 1 )
            {
                uint tNumEdges = this->get_num_edges();

                // Access entities stored in mesh database
                Matrix< IdMat > tEdgeIds = this->get_entities_universal_glob_id( EntityRank::EDGE );

                // resize member variable to its right size
                mSTKMeshData->mEntityLocaltoGlobalMap( 1 ) = Matrix< IdMat >( tNumEdges, 1 );

                // Populate internal member variable that contains the local index to
                // global id node communication information
                for ( uint iEdge = 0; iEdge < tNumEdges; ++iEdge )
                {
                    // local to global and owner processor
                    mSTKMeshData->mEntityLocaltoGlobalMap( 1 )( iEdge ) = tEdgeIds( iEdge );
                }
            }

            if ( mSTKMeshData->mNumDims > 2 )
            {
                uint tNumFaces = this->get_num_faces();

                // Access entities stored in mesh database
                Matrix< IdMat > tFaceIds = this->get_entities_universal_glob_id( EntityRank::FACE );

                // resize member variable to its right size
                mSTKMeshData->mEntityLocaltoGlobalMap( 2 ) = Matrix< IdMat >( tNumFaces, 1 );

                // Populate internal member variable that contains the local index to
                // global id element communication information
                for ( uint iFace = 0; iFace < tNumFaces; ++iFace )
                {
                    // local to global and owner processor
                    mSTKMeshData->mEntityLocaltoGlobalMap( 2 )( iFace ) = tFaceIds( iFace );
                }
            }
        }

        // ----------------------------------------------------------------------------

        // Function to create edges and faces communication lists in parallel
        void
        Mesh_Core_STK::create_shared_communication_lists()
        {
            // Get basic mesh information
            Matrix< IdMat > tNodesShared    = this->get_entities_glb_shared_current_proc( EntityRank::NODE );
            uint            tNumNodesShared = tNodesShared.length();

            // Generate list of processors sharing information
            // -----------------------------------------------

            std::vector< uint > tActiveSharedProcs;

            // Loop over the number of nodes shared to get the shared processors
            for ( uint iNodeShared = 0; iNodeShared < tNumNodesShared; ++iNodeShared )
            {
                Matrix< IdMat > tProcsSharing = this->get_procs_sharing_entity_by_id( tNodesShared( iNodeShared ), EntityRank::NODE );

                for ( uint iProc = 0; iProc < tProcsSharing.length(); ++iProc )
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
            for ( uint iProcShared = 0; iProcShared < aNumActiveSharedProcs; ++iProcShared )
            {
                mSTKMeshData->mProcsSharedToIndex.insert( std::pair< uint, uint >( tActiveSharedProcs.at( iProcShared ), iProcShared ) );
            }

            // Generate nodes shared per processor list
            // ----------------------------------------
            //        mSTKMeshData->mNodeMapToSharingProcs = this->get_shared_info_by_entity( aNumActiveSharedProcs, EntityRank::NODE );
            //
            //        // Generate elements shared per processor list (because of aura)
            //        // -------------------------------------------
            //        mSTKMeshData->mElemMapToSharingProcs = this->get_shared_info_by_entity( aNumActiveSharedProcs, EntityRank::ELEMENT );
            //
            //        // Generate edges shared per processor list
            //        // ----------------------------------------
            //        mSTKMeshData->mEdgeMapToSharingProcs = this->get_shared_info_by_entity( aNumActiveSharedProcs, EntityRank::EDGE );
            //
            //        // Generate faces shared per processor list
            //        // ----------------------------------------
            //        mSTKMeshData->mFaceMapToSharingProcs = this->get_shared_info_by_entity( aNumActiveSharedProcs, EntityRank::FACE );
        }

        // ----------------------------------------------------------------------------

        // Add all node sets information to database
        void
        Mesh_Core_STK::process_node_sets(
                MtkMeshData aMeshData )
        {
            // Declare basic node set information
            uint                  tNumNodeSets = aMeshData.SetsInfo->get_num_node_sets();
            stk::mesh::EntityRank aStkSetRank  = stk::topology::NODE_RANK;

            for ( uint iSet = 0; iSet < tNumNodeSets; ++iSet )
            {
                // STK interface variables declaration
                stk::mesh::Part*      aSetPart = mSTKMeshData->mMtkMeshMetaData->get_part( mSTKMeshData->mSetNames[ 0 ][ iSet ] );
                stk::mesh::PartVector aAddPart( 1, aSetPart );
                stk::mesh::EntityId   aGlobalId;
                stk::mesh::Entity     aEntity;

                // Get node set and size of node set
                MtkNodeSetInfo* tNodeSet       = aMeshData.SetsInfo->get_node_set( iSet );
                uint            tNumNodesInSet = tNodeSet->mNodeIds->numel();

                // Populate node sets (change entity parts if nodes were declared already)
                for ( uint iEntity = 0; iEntity < tNumNodesInSet; ++iEntity )
                {
                    // Declare new entity or add existing entity to declared part
                    aGlobalId = ( *tNodeSet->mNodeIds )( iEntity );
                    aEntity   = mSTKMeshData->mMtkMeshBulkData->get_entity( aStkSetRank, aGlobalId );

                    if ( !mSTKMeshData->mMtkMeshBulkData->is_valid( aEntity ) )
                    {
                        //                    aEntity = mSTKMeshData->mMtkMeshBulkData->declare_entity( aStkSetRank, aGlobalId, aAddPart );
                    }
                    else
                    {
                        mSTKMeshData->mMtkMeshBulkData->change_entity_parts( aEntity, aAddPart );
                    }
                }    // end of node sets declarations
            }
        }

        // ----------------------------------------------------------------------------

        // Add all side sets information to database
        void
        Mesh_Core_STK::process_side_sets(
                MtkMeshData aMeshData )
        {
            // Get all sets provided by the user
            uint tNumSideSets = aMeshData.SetsInfo->get_num_side_sets();

            ///////////////////////////////
            // Begin modification cycle  //
            ///////////////////////////////
            mSTKMeshData->mMtkMeshBulkData->modification_begin();

            for ( uint iSet = 0; iSet < tNumSideSets; ++iSet )
            {
                // STK interface variables declaration
                stk::mesh::Part*      aSetPart = mSTKMeshData->mMtkMeshMetaData->get_part( mSTKMeshData->mSetNames[ 1 ][ iSet ] );
                stk::mesh::PartVector aAddPart( 1, aSetPart );
                stk::mesh::EntityId   aGlobalElemId;
                stk::mesh::Entity     aElemEntity;
                uint                  tRequestedSideOrd;

                // Get side set and size of node set
                MtkSideSetInfo* tSideSet       = aMeshData.SetsInfo->get_side_set( iSet );
                uint            tNumSidesInSet = tSideSet->mElemIdsAndSideOrds->n_rows();

                // moris::print((*tSideSet->mElemIdsAndSideOrds),"(*tSideSet->mElemIdsAndSideOrds)(");

                for ( uint iEntity = 0; iEntity < tNumSidesInSet; ++iEntity )
                {
                    // First column contains element ids that will later be match with faces
                    aGlobalElemId     = ( *tSideSet->mElemIdsAndSideOrds )( iEntity, 0 );
                    aElemEntity       = mSTKMeshData->mMtkMeshBulkData->get_entity( stk::topology::ELEMENT_RANK, aGlobalElemId );
                    tRequestedSideOrd = ( *tSideSet->mElemIdsAndSideOrds )( iEntity, 1 );

                    if ( mSTKMeshData->mMtkMeshBulkData->is_valid( aElemEntity ) )
                    {
                        // Create side entity
                        mSTKMeshData->mMtkMeshBulkData->declare_element_side( aElemEntity, tRequestedSideOrd, aAddPart );
                    }
                }    // end of side sets declarations
            }

            ///////////////////////////////
            // Close modification cycle  //
            ///////////////////////////////
            mSTKMeshData->mMtkMeshBulkData->modification_end();
        }

        // ----------------------------------------------------------------------------
        // Parallel specific implementation for blocks in database

        void
        Mesh_Core_STK::populate_mesh_database_serial(
                moris::uint                          aElementTypeInd,
                MtkMeshData                          aMeshData,
                std::vector< stk::mesh::PartVector > aElemParts,
                Matrix< IdMat >                      aOwnerPartInds )
        {

            uint tNumElems        = aMeshData.ElemConn( aElementTypeInd )->n_rows();
            uint aNumNodesPerElem = aMeshData.ElemConn( aElementTypeInd )->n_cols();

            // Declare variables to access connectivity
            Matrix< IdMat >           tDummyMat( 1, aNumNodesPerElem );
            stk::mesh::EntityIdVector aCurrElemConn( aNumNodesPerElem );
            stk::mesh::EntityId       aElemGlobalId;

            // Loop over the number of elements and interface between MORIS and Stk for connectivity
            for ( uint iElem = 0; iElem < tNumElems; ++iElem )
            {
                // Get row of nodes connected in moris variable and assign to STK variable
                tDummyMat.set_row( 0, aMeshData.ElemConn( aElementTypeInd )->get_row( iElem ) );
                aCurrElemConn.assign( tDummyMat.data(), tDummyMat.data() + aNumNodesPerElem );

                // Declare element in function that also declares element-node relations internally
                aElemGlobalId = ( *aMeshData.LocaltoGlobalElemMap( aElementTypeInd ) )( iElem );
                stk::mesh::declare_element( *mSTKMeshData->mMtkMeshBulkData, aElemParts[ aOwnerPartInds( iElem ) ], aElemGlobalId, aCurrElemConn );
            }
        }

        // ----------------------------------------------------------------------------

        // Access set entity ids
        moris::Matrix< IdMat >
        Mesh_Core_STK::get_set_entity_glob_ids(
                stk::mesh::EntityRank aEntityRank,
                std::string           aSetName ) const
        {
            // Get pointer to field defined by input name
            stk::mesh::Part* const tSetPart = mSTKMeshData->mMtkMeshMetaData->get_part( aSetName );

            MORIS_ASSERT( tSetPart != nullptr, "Set not found. Double check name provided." );

            // Access data through a selector
            stk::mesh::Selector     tSetSelector( *tSetPart );
            stk::mesh::EntityVector aEntities;
            stk::mesh::get_selected_entities( tSetSelector, mSTKMeshData->mMtkMeshBulkData->buckets( aEntityRank ), aEntities );

            // Get entity Ids
            uint            tNumEntities = aEntities.size();
            Matrix< IdMat > tOutputEntityIds( tNumEntities, 1 );
            for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
            {
                tOutputEntityIds( iEntity ) = (moris_id)mSTKMeshData->mMtkMeshBulkData->identifier( aEntities[ iEntity ] );
            }

            return tOutputEntityIds;
        }

        // ----------------------------------------------------------------------------

        //##############################################
        // internal id functions
        //##############################################
        std::vector< stk::mesh::Entity >
        Mesh_Core_STK::entities_connected_to_entity_stk(
                stk::mesh::Entity* const    aInputEntity,
                stk::mesh::EntityRank const aInputEntityRank,
                stk::mesh::EntityRank const aOutputEntityRank ) const
        {
            // Declare the object where we are going to store the shared faces and handlers
            std::vector< stk::mesh::Entity > tDesiredEntitiesConnectedToInputEntities;

            // Check if the connectivity exists (i.e., was already generated and is stored in mesh data)
            //        if (mSTKMeshData->mMtkMeshBulkData->connectivity_map().valid(aInputEntityRank, aOutputEntityRank))
            //        {

            switch ( aOutputEntityRank )
            {
                case stk::topology::NODE_RANK:

                    // Fill entities connected
                    if ( mSTKMeshData->mMtkMeshBulkData->num_nodes( aInputEntity[ 0 ] ) > 0 )
                    {
                        // Get pointers to the location of the connected nodes
                        stk::mesh::Entity const * tDesiredEntityStart = mSTKMeshData->mMtkMeshBulkData->begin_nodes( aInputEntity[ 0 ] );
                        stk::mesh::Entity const * tDesiredEntityEnd   = mSTKMeshData->mMtkMeshBulkData->end_nodes( aInputEntity[ 0 ] );

                        // Store faces in output vector
                        tDesiredEntitiesConnectedToInputEntities.assign( tDesiredEntityStart, tDesiredEntityEnd );
                    }
                    break;

                case stk::topology::EDGE_RANK:

                    // Fill entities connected
                    if ( mSTKMeshData->mMtkMeshBulkData->num_edges( aInputEntity[ 0 ] ) > 0 )
                    {
                        // Get pointers to the location of the connected edges
                        stk::mesh::Entity const * tDesiredEntityStart = mSTKMeshData->mMtkMeshBulkData->begin_edges( aInputEntity[ 0 ] );
                        stk::mesh::Entity const * tDesiredEntityEnd   = mSTKMeshData->mMtkMeshBulkData->end_edges( aInputEntity[ 0 ] );

                        // Store faces in output vector
                        tDesiredEntitiesConnectedToInputEntities.assign( tDesiredEntityStart, tDesiredEntityEnd );
                    }
                    break;

                case stk::topology::FACE_RANK:

                    // Fill entities connected
                    if ( mSTKMeshData->mMtkMeshBulkData->num_faces( aInputEntity[ 0 ] ) > 0 )
                    {
                        // Get pointers to the location of the connected faces
                        stk::mesh::Entity const * tDesiredEntityStart = mSTKMeshData->mMtkMeshBulkData->begin_faces( aInputEntity[ 0 ] );
                        stk::mesh::Entity const * tDesiredEntityEnd   = mSTKMeshData->mMtkMeshBulkData->end_faces( aInputEntity[ 0 ] );

                        // Store faces in output vector
                        tDesiredEntitiesConnectedToInputEntities.assign( tDesiredEntityStart, tDesiredEntityEnd );
                    }
                    break;

                case stk::topology::ELEMENT_RANK:

                    // Fill entities connected
                    if ( mSTKMeshData->mMtkMeshBulkData->num_elements( aInputEntity[ 0 ] ) > 0 )
                    {
                        // Get pointers to the location of the connected elements
                        stk::mesh::Entity const * tDesiredEntityStart = mSTKMeshData->mMtkMeshBulkData->begin_elements( aInputEntity[ 0 ] );
                        stk::mesh::Entity const * tDesiredEntityEnd   = mSTKMeshData->mMtkMeshBulkData->end_elements( aInputEntity[ 0 ] );

                        // Store faces in output vector
                        tDesiredEntitiesConnectedToInputEntities.assign( tDesiredEntityStart, tDesiredEntityEnd );
                    }
                    break;

                default:
                    std::cerr << " wrong topology in entities_connected_to_entity_stk ";
                    break;
            }
            //        }
            //        else
            //        {
            //            std::cerr << " STK already has valid connectivity maps. Check if you are trying to access invalid connectivity (e.g., edge to edge)";
            //        }
            return tDesiredEntitiesConnectedToInputEntities;
        }

        // ----------------------------------------------------------------------------

        bool
        Mesh_Core_STK::is_aura_cell( moris_id aElementId ) const
        {
            return false;
        }

        // ----------------------------------------------------------------------------

    }    // namespace mtk
}    // namespace moris
