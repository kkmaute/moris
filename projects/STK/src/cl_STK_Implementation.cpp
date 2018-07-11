
// MORIS library header files.
#include "cl_STK_Implementation.hpp" // STK/src
#include "ios.hpp"
#include "ioUtils.hpp"
#include "fn_find.hpp" // LNA/src

//#include <Teuchos_RCP.hpp>              // for RCP::RCP<T>, RCP::operator*, etc
//#include "Teuchos_RCPDecl.hpp"          // for RCP

#include <stk_util/environment/FileUtils.hpp>
#include <stk_io/InputFile.hpp>   // for DatabasePurpose
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose
#include <stk_io/StkMeshIoBroker.hpp>        // for stk broker
#include <stk_io/IossBridge.hpp>        // for FieldAndName, STKIORequire, etc
#include <Ionit_Initializer.h>     // for Initializer

//stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp
//#include <stk_mesh/base/Entity.hpp>
//#include <stk_mesh/base/EntityComm.hpp>

#include <Ioss_ElementBlock.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_SideSet.h>
#include <Ioss_SideBlock.h>
#include <Ioss_NodeSet.h>

#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region, NodeBlockContainer

// Class member functions
/**
 * mesh_io reads a string with the name of the file containing the mesh
 * and creates an stk mesh from it.
 *
 */

// Destructor
// ----------
moris::STK_Implementation::~STK_Implementation()
{
    // Delete member variables that contain the mesh database
    delete mMeshReader;
    delete mMtkMeshBulkData;
    delete mMtkMeshMetaData;
}

// -----------------------------------------------------------------------------------

// Constructor for stk meshes generated internally by stk or from input (exodus) files
moris::STK_Implementation::STK_Implementation(
        std::string    aFileName,
        MtkSetsInfo*   aSetsInfo ): database(), mSetRankFlags({ false, false, false})
{
    // Define the path prefix as an empty string
    std::string tPrefix;

    // Check if aFileName contains a ':' which indicates that the mesh will be generated internally
    size_t tColon = aFileName.find(':');
    if ( tColon == std::string::npos )
    {
        tPrefix = "./";
    }

    // Call the function that handles the communication between stk and moris
    this->build_mesh( tPrefix+aFileName, aSetsInfo );
}

// ----------------------------------------------------------------------------

// Interface with STK to generate the mesh database
void
moris::STK_Implementation::build_mesh(
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

    // Get number of spatial dimensions for this mesh
    mNumDims = mMtkMeshMetaData->spatial_dimension();

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
    uint tProcSize = par_size();
    if ( tProcSize > 1 )
    {
        this->create_communication_lists_from_string();
    }
}

// -----------------------------------------------------------------

// Constructors for stk meshes generated from data given by the user
moris::STK_Implementation::STK_Implementation(
        MtkMeshData   aMeshData ):
                                  database(), mSetRankFlags( { false, false, false} )
{
    // Flag for handling data generated mesh
    mDataGeneratedMesh = true;

    // Verify the minimum required data was provided and initialize uninitialized variables if necessary
    this->check_and_update_input_data( aMeshData );

    // Call the function that handles the communication between stk and moris
    this->build_mesh( aMeshData );
}

// ----------------------------------------------------------------------------

// Verifications of mesh essential information provided for meshes generated from data.
void
moris::STK_Implementation::check_and_update_input_data(
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

        mLocalToGlobalNodeMap = aMeshData.LocaltoGlobalNodeMap[0];
    }
    else if ( ( aMeshData.LocaltoGlobalNodeMap == NULL ) && ( tProcSize == 1 ) )
    {
        // Generate nodes maps
        mLocalToGlobalNodeMap.resize( tNumNodes, 1 );
        for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
        {
            mLocalToGlobalNodeMap( iNode ) = iNode + 1;
        }

        aMeshData.LocaltoGlobalNodeMap = &mLocalToGlobalNodeMap;
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

        mLocalToGlobalElemMap = aMeshData.LocaltoGlobalElemMap[0];
    }
    else
    {
        // Threshold value and number of elements in current processor
        uint tThresholdValue = 1000000;

        // Check that no one processor exceeds the threshold value.
        MORIS_ASSERT( tNumElems <= tThresholdValue, "Constructor generates a maximum of 1000000 elements per processor in current constructor. "
                                     "Change Threshold value or provide element local to global map." );

        // Generate IDs assuming there are no more than tThresholdValue elements per processor
        mLocalToGlobalElemMap.resize( tNumElems, 1 );
        for ( uint iElem = 0; iElem < tNumElems; ++iElem )
        {
            mLocalToGlobalElemMap( iElem ) = tThresholdValue*tProcRank + iElem + 1;
        }

        aMeshData.LocaltoGlobalElemMap = &mLocalToGlobalElemMap;
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

// Verifications of set data arrangement and consistency for meshes generated from data.
void
moris::STK_Implementation::check_and_update_sets_data(
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

// ----------------------------------------------------------------------------

// Verification of fields data arrangement and consistency for meshes generated from data.
void
moris::STK_Implementation::check_and_update_fields_data(
        MtkMeshData&   aMeshData )
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
        Mat< uint > tSupFieldComp = { {1, 2, 3, 4, 9} };
        uint tNumFieldComp = aMeshData.FieldsInfo[0].FieldsData[0]( iField ).n_cols();
        Mat< uint > tDummy = ( tSupFieldComp == tNumFieldComp );
        Mat< uint > tCompFound = find ( tDummy );

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

// Main interface with STK that include calls to functions that provide specific implementation details.
void
moris::STK_Implementation::build_mesh(
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
moris::STK_Implementation::declare_mesh_parts(
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

// Provide element type (Hex8, Tri3, etc) and throw error if element is not supported yet.
stk::topology::topology_t
moris::STK_Implementation::get_mesh_topology(
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

// Second declaration to structure the database before filling the data
void
moris::STK_Implementation::declare_mesh_fields(
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

// ----------------------------------------------------------------------------

// Declare size of a field (per entity) and throw an error if it is not supported
void
moris::STK_Implementation::internal_declare_mesh_field(
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
moris::STK_Implementation::populate_mesh_database(
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
moris::STK_Implementation::process_block_sets(
        MtkMeshData   aMeshData )
{
    // Get all sets provided by the user and go to the block set
    uint tNumElems     = aMeshData.ElemConn[0].size( 0 );
    uint tNumBlockSets = 1;

    Mat< uint > aOwnerPartInds( tNumElems, 1, 0 );
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

// Parallel specific implementation for blocks in database
void
moris::STK_Implementation::populate_mesh_database_serial(
        MtkMeshData                            aMeshData,
        std::vector< stk::mesh::PartVector >   aElemParts,
        Mat< uint >                            aOwnerPartInds )
{
    uint tNumElems        = aMeshData.ElemConn[0].size( 0 );
    uint aNumNodesPerElem = aMeshData.ElemConn[0].size( 1 );

    // Declare variables to access connectivity
    Mat< uint > tDummyMat( 1, aNumNodesPerElem );
    stk::mesh::EntityIdVector aCurrElemConn( aNumNodesPerElem );
    stk::mesh::EntityId aElemGlobalId;

    // Loop over the number of elements and interface between MORIS and Stk for connectivity
    for ( uint iElem = 0; iElem < tNumElems; ++iElem )
    {
        // Get row of nodes connected in moris variable and assign to STK variable
        tDummyMat.row( 0 ) = aMeshData.ElemConn->row( iElem );
        aCurrElemConn.assign( mem_pointer( tDummyMat ), mem_pointer( tDummyMat ) + aNumNodesPerElem );

        // Declare element in function that also declares element-node relations internally
        aElemGlobalId = aMeshData.LocaltoGlobalElemMap[0]( iElem );
        stk::mesh::declare_element( *mMtkMeshBulkData, aElemParts[aOwnerPartInds( iElem )], aElemGlobalId, aCurrElemConn );
    }
}

// ----------------------------------------------------------------------------

// Parallel specific implementation for blocks in database
void
moris::STK_Implementation::populate_mesh_database_parallel(
        MtkMeshData                           aMeshData,
        std::vector< stk::mesh::PartVector >   aPartBlocks,
        Mat< uint >                            aOwnerPartInds )
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
    Mat< uint > tDummyMat( 1, aNumNodesPerElem );
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
            tDummyMat.row( 0 ) = aMeshData.ElemConn->row( iElem );
            aCurrElemConn.assign( mem_pointer( tDummyMat ), mem_pointer( tDummyMat ) + aNumNodesPerElem );
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
            if ( aMeshData.EntProcOwner[0]( iNode ) != tProcRank )
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
            if ( aMeshData.EntProcOwner[0]( iElem ) == tProcRank )
            {
                // Get row of nodes connected in moris variable and assign to STK variable
                tDummyMat.row( 0 ) = aMeshData.ElemConn->row( iElem );
                aCurrElemConn.assign( mem_pointer( tDummyMat ), mem_pointer( tDummyMat ) + aNumNodesPerElem );

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
            if ( aMeshData.EntProcOwner[0]( iElem )!= tProcRank )
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

// Add all node sets information to database
void
moris::STK_Implementation::process_node_sets(
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
moris::STK_Implementation::process_side_sets(
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

// Add all fields information to database
void
moris::STK_Implementation::populate_mesh_fields(
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
            Mat< uint > tFieldIds;
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

// Tell the STK I/O broker to create an output mesh
void
moris::STK_Implementation::create_output_mesh(
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
        double tTime = 0.0;
        mMeshReader->begin_output_step( outputFileIdx, tTime );
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

// ----------------------------------------------------------------------------

// Access entities in selected portion of the mesh
moris::Mat< uint >
moris::STK_Implementation::get_entities_in_selector_interface(
        EntityRank            aRequestedEntityRank,
        stk::mesh::Selector   aSelectedEntities ) const
{
    // Get selected entities
    std::vector<stk::mesh::Entity>tOutputEntityIDs;
    stk::mesh::get_selected_entities(
            aSelectedEntities, mMtkMeshBulkData->buckets( this->get_stk_entity_rank( aRequestedEntityRank ) ), tOutputEntityIDs );

    // Interface with STK, get the Ids and return them
    uint tNumEntity = tOutputEntityIDs.size();
    Mat< uint > tOutputEntityIDMat( tNumEntity, 1, 0 );

    for ( uint i = 0; i<tNumEntity; i++)
    {
        tOutputEntityIDMat( i, 0 ) = ( uint ) mMtkMeshBulkData->identifier( tOutputEntityIDs[i] );
    }

    return tOutputEntityIDMat;
}

// ----------------------------------------------------------------------------

// Access number of entities of a particular rank for the element being used (e.g., Hex8 has 6 sides)
uint
moris::STK_Implementation::get_entity_topology_num_entities(
        stk::mesh::EntityRank   aEntityRank,
        stk::mesh::EntityRank   aRequestedRank,
        uint                    aEntityId ) const
{
    // Get topology of the current entity
    stk::topology tEntityTopo = mMtkMeshBulkData->bucket( mMtkMeshBulkData->get_entity( aEntityRank, aEntityId ) ).topology();

    // Check if the topology is valid
    MORIS_ASSERT( tEntityTopo.is_valid(), "Invalid entity topology.");

    // Return number of nodes extracted from entity counts vector
    uint tNumEntities = 0;
    switch ( aRequestedRank )
    {
        case stk::topology::NODE_RANK: // for node sets
        {
            tNumEntities = tEntityTopo.num_nodes();
            break;
        }
        case stk::topology::EDGE_RANK: // for side sets
        {
            tNumEntities = tEntityTopo.num_edges();
            break;
        }
        case stk::topology::FACE_RANK: // for side sets
        {
            tNumEntities = tEntityTopo.num_faces();
            break;
        }
        default:
        {
            MORIS_ASSERT( 0, "Invalid rank provided in get_entity_topology_num_entities." );
            break;
        }
    }
    return tNumEntities;
}

// ----------------------------------------------------------------------------

// General implementation for accessing connectivity of any entity to another entity type.
// Note: Entity ranks should be different to be able to access connectivity (e.g., face to face not supported).
moris::Mat< uint >
moris::STK_Implementation::entities_connected_to_entity_generic(
        const stk::mesh::Entity*   aInputEntity,
        stk::mesh::EntityRank      aNeedConnectivityOfType ) const
{
    // Declare the object where we are going to store the shared faces and handlers
    std::vector<stk::mesh::Entity>  tDesiredEntitiesConnectedToInputEntities;

    // Define the type of connectivity
    stk::mesh::EntityRank aInputEntityType = mMtkMeshBulkData->entity_rank( *aInputEntity );

    // Check if the connectivity exists (i.e., was already generated and is stored in mesh)
    if ( mMtkMeshBulkData->connectivity_map().valid( aInputEntityType, aNeedConnectivityOfType ) )
    {

        switch ( aNeedConnectivityOfType )
        {
        case stk::topology::NODE_RANK:
        {
            // Fill entities connected
            if ( mMtkMeshBulkData->num_nodes( aInputEntity[0] ) > 0 )
            {
                // Get pointers to the location of the connected nodes
                stk::mesh::Entity const* tDesiredEntityStart = mMtkMeshBulkData->begin_nodes( aInputEntity[0] );
                stk::mesh::Entity const* tDesiredEntityEnd   = mMtkMeshBulkData->end_nodes( aInputEntity[0] );

                // Store faces in output vector
                tDesiredEntitiesConnectedToInputEntities.assign ( tDesiredEntityStart, tDesiredEntityEnd );
            }
            break;
        }
        case stk::topology::EDGE_RANK:
        {
            // Fill entities connected
            if ( mMtkMeshBulkData->num_edges( aInputEntity[0] ) > 0 )
            {
                // Get pointers to the location of the connected edges
                stk::mesh::Entity const* tDesiredEntityStart = mMtkMeshBulkData->begin_edges( aInputEntity[0] );
                stk::mesh::Entity const* tDesiredEntityEnd   = mMtkMeshBulkData->end_edges( aInputEntity[0] );

                // Store faces in output vector
                tDesiredEntitiesConnectedToInputEntities.assign ( tDesiredEntityStart, tDesiredEntityEnd );
            }
            break;
        }
        case stk::topology::FACE_RANK:
        {
            // Fill entities connected
            if ( mMtkMeshBulkData->num_faces( aInputEntity[0] ) > 0 )
            {
                // Get pointers to the location of the connected faces
                stk::mesh::Entity const* tDesiredEntityStart = mMtkMeshBulkData->begin_faces( aInputEntity[0] );
                stk::mesh::Entity const* tDesiredEntityEnd   = mMtkMeshBulkData->end_faces( aInputEntity[0] );

                // Store faces in output vector
                tDesiredEntitiesConnectedToInputEntities.assign ( tDesiredEntityStart, tDesiredEntityEnd );
            }
            break;
        }
        case stk::topology::ELEMENT_RANK:
        {
            // Fill entities connected
            if ( mMtkMeshBulkData->num_elements( aInputEntity[0] ) > 0 )
            {
                // Get pointers to the location of the connected elements
                stk::mesh::Entity const* tDesiredEntityStart = mMtkMeshBulkData->begin_elements( aInputEntity[0] );
                stk::mesh::Entity const* tDesiredEntityEnd   = mMtkMeshBulkData->end_elements( aInputEntity[0] );

                // Store faces in output vector
                tDesiredEntitiesConnectedToInputEntities.assign ( tDesiredEntityStart, tDesiredEntityEnd );
            }
            break;
        }
        default:
        {
            MORIS_ASSERT( 0, "Incorrect topology in entities_connected_to_entity_generic_interface." );
            break;
        }
        }
    }
    else
    {
        MORIS_ASSERT( 0, "Check if you are trying to access invalid connectivity (e.g., edge to edge)" );
    }

    uint tNumNodes = tDesiredEntitiesConnectedToInputEntities.size();
    Mat< uint > tEntitiesConnectedToGivenEntity( tNumNodes, 1 );

    // Store the entity ID in output matrix
    for ( uint iNode=0; iNode < tNumNodes; ++iNode )
    {
        tEntitiesConnectedToGivenEntity( iNode ) = ( uint ) mMtkMeshBulkData->identifier( tDesiredEntitiesConnectedToInputEntities[iNode] );
    }

    return tEntitiesConnectedToGivenEntity;
}

// ----------------------------------------------------------------------------

// The generic call of entities connected to entities is not valid for entities of the same type.
// For those situations, a two step process is needed. If we want to access, for example, elements
// connected to elements, first element nodes need to be pulled to later ask these nodes what
// elements are connected to them.
moris::Mat< uint >
moris::STK_Implementation::get_elements_connected_to_element(
        uint   const aElemId ) const
{
    // First get faces connected to element
    const stk::mesh::Entity aElemEntity = mMtkMeshBulkData->get_entity( stk::topology::ELEMENT_RANK, aElemId );
    Mat< uint > tFacesInElem            = this->entities_connected_to_entity_generic( &aElemEntity, stk::topology::FACE_RANK );

    MORIS_ASSERT( ( tFacesInElem.n_rows() != 0 ) || ( tFacesInElem.n_cols() != 0 ),
            "No faces connected to element found. Maybe the CreateAllEdgesAndFaces flag is set to false. Check mesh struct." );

    // Then for each face get elements connected
    uint tCounter  = 0;
    uint tNumFaces = tFacesInElem.length();

    Mat< uint > tElemsConnectedToElem ( 2 * tNumFaces, 1, 0 );

    for ( uint faceIt = 0; faceIt < tNumFaces; ++faceIt )
    {
        uint aFaceId = tFacesInElem( faceIt );
        const stk::mesh::Entity aFaceEntity = mMtkMeshBulkData->get_entity( stk::topology::FACE_RANK, aFaceId );
        Mat< uint > tDummyConnectivity      = this->entities_connected_to_entity_generic( &aFaceEntity, stk::topology::ELEMENT_RANK );

        // Faces in mesh boundaries do not have more than one element
        if ( tDummyConnectivity.length() > 0 )
        {
            if ( tDummyConnectivity( 0 ) != aElemId )
            {
                tElemsConnectedToElem( tCounter ) = tDummyConnectivity( 0 );
                tCounter++;
            }
        }

        if ( tDummyConnectivity.length() > 1 )
        {
            if ( tDummyConnectivity( 1 ) != aElemId )
            {
                tElemsConnectedToElem( tCounter ) = tDummyConnectivity( 1 );
                tCounter++;
            }
        }

        MORIS_ASSERT( tDummyConnectivity.length() <= 2,
                "For some reason face has more than 2 elements connected to it... Check get_elements_connected_to_element." );
    }

    // Resize to include only ids added above and get rid of initialized extra zeros
    tElemsConnectedToElem.resize( tCounter, 1 );

    return tElemsConnectedToElem;
}

// ----------------------------------------------------------------------------

// Access entity ordinal
moris::Mat< uint >
moris::STK_Implementation::entity_ordinals_connected_to_given_entity(
        uint         const aEntityId,
        EntityRank   const aInputEntityRank,
        EntityRank   const aOutputEntityRank ) const
{
    // Define entity
    stk::mesh::EntityRank aInputRank   = this->get_stk_entity_rank( aInputEntityRank );
    const stk::mesh::Entity aStkEntity = mMtkMeshBulkData->get_entity( aInputRank, aEntityId );

    // Call function that gets the connected entities
    stk::mesh::EntityRank aOutputRank = this->get_stk_entity_rank( aOutputEntityRank );
    return this->entity_ordinals_connected_to_entity_generic( &aStkEntity, aOutputRank );
}

// ----------------------------------------------------------------------------

// similar to get_elements_connected_to_element, but with ordinals instead
moris::Mat< uint >
moris::STK_Implementation::entity_ordinals_connected_to_entity_generic(
        const stk::mesh::Entity*   aInputEntity,
        stk::mesh::EntityRank      aNeedConnectivityOfType ) const
{
    // Declare the object where we are going to store the shared faces and handlers
    std::vector< stk::mesh::ConnectivityOrdinal >  tDesiredEntityOrdinalsConnectedToInputEntities;

    // Define the type of connectivity
    stk::mesh::EntityRank aInputEntityType = mMtkMeshBulkData->entity_rank( *aInputEntity );

    // Check if the connectivity exists (i.e., was already generated and is stored in mesh)
    if ( mMtkMeshBulkData->connectivity_map().valid( aInputEntityType, aNeedConnectivityOfType ) )
    {
        switch ( aNeedConnectivityOfType )
        {
        case stk::topology::NODE_RANK:
        {
            // Fill entities connected
            if ( mMtkMeshBulkData->num_nodes( aInputEntity[0] ) > 0 )
            {
                // Get pointers to the location of the connected nodes
                const stk::mesh::ConnectivityOrdinal* tDesiredEntityStart = mMtkMeshBulkData->begin_node_ordinals( aInputEntity[0] );
                const stk::mesh::ConnectivityOrdinal* tDesiredEntityEnd   = mMtkMeshBulkData->end_node_ordinals( aInputEntity[0] );

                // Store faces in output vector
                tDesiredEntityOrdinalsConnectedToInputEntities.assign ( tDesiredEntityStart, tDesiredEntityEnd );
            }
            break;
        }
        case stk::topology::EDGE_RANK:
        {
            // Fill entities connected
            if ( mMtkMeshBulkData->num_edges( aInputEntity[0] ) > 0 )
            {
                // Get pointers to the location of the connected edges
                const stk::mesh::ConnectivityOrdinal* tDesiredEntityStart = mMtkMeshBulkData->begin_edge_ordinals( aInputEntity[0] );
                const stk::mesh::ConnectivityOrdinal* tDesiredEntityEnd   = mMtkMeshBulkData->end_edge_ordinals( aInputEntity[0] );

                // Store faces in output vector
                tDesiredEntityOrdinalsConnectedToInputEntities.assign ( tDesiredEntityStart, tDesiredEntityEnd );
            }
            break;
        }
        case stk::topology::FACE_RANK:
        {
            // Fill entities connected
            if ( mMtkMeshBulkData->num_faces( aInputEntity[0] ) > 0 )
            {
                // Get pointers to the location of the connected faces
                const stk::mesh::ConnectivityOrdinal* tDesiredEntityStart = mMtkMeshBulkData->begin_face_ordinals( aInputEntity[0] );
                const stk::mesh::ConnectivityOrdinal* tDesiredEntityEnd   = mMtkMeshBulkData->end_face_ordinals( aInputEntity[0] );

                // Store faces in output vector
                tDesiredEntityOrdinalsConnectedToInputEntities.assign ( tDesiredEntityStart, tDesiredEntityEnd );
            }
            break;
        }
        case stk::topology::ELEMENT_RANK:
        {
            // Fill entities connected
            if ( mMtkMeshBulkData->num_elements( aInputEntity[0] ) > 0 )
            {
                // Get pointers to the location of the connected elements
                const stk::mesh::ConnectivityOrdinal* tDesiredEntityStart = mMtkMeshBulkData->begin_element_ordinals( aInputEntity[0] );
                const stk::mesh::ConnectivityOrdinal* tDesiredEntityEnd   = mMtkMeshBulkData->end_element_ordinals( aInputEntity[0] );

                // Store faces in output vector
                tDesiredEntityOrdinalsConnectedToInputEntities.assign ( tDesiredEntityStart, tDesiredEntityEnd );
            }
            break;
        }
        default:
        {
            MORIS_ASSERT( 0, "Incorrect topology in entities_connected_to_entity_generic " );
            break;
        }
        }
    }

    uint tNumNodes = tDesiredEntityOrdinalsConnectedToInputEntities.size();
    Mat< uint > entityOrdinalsConnectedToGivenEntity( tNumNodes, 1 );

    // Store the entity ID in output celltDesiredEntityOrdinalsConnectedToInputEntities
    for ( uint iNode=0; iNode < tNumNodes; ++iNode )
    {
        entityOrdinalsConnectedToGivenEntity( iNode ) = ( uint ) tDesiredEntityOrdinalsConnectedToInputEntities[iNode];
    }

    return entityOrdinalsConnectedToGivenEntity;
}

// ----------------------------------------------------------------------------

// Access set entities
moris::Mat< uint >
moris::STK_Implementation::get_entities_in_set_generic(
        stk::mesh::EntityRank   aRequiredEntityRank,
        stk::mesh::EntityRank   aSetRank,
        uint                    const aSetId ) const
{
    // Generate the name of the required set
    std::string aSetName;

    switch ( aSetRank )
    {
    case stk::topology::NODE_RANK: // for node sets
    {
        // Check if sets were defined by the user
        if ( mSetRankFlags[0] )
        {
            aSetName = ( mSetNames[0] )[aSetId-1];
        }
        else
        {
            aSetName = "nodelist_" + std::to_string( aSetId );
        }

        break;
    }
    case stk::topology::EDGE_RANK: // for side sets
    case stk::topology::FACE_RANK: // for side sets
    {
        // Check if sets were defined by the user
        if ( mSetRankFlags[1])
        {
            aSetName = ( mSetNames[1] )[aSetId-1];
        }
        else
        {
            aSetName = "surface_" + std::to_string( aSetId ); // NOTE not sure about the name by default in Exodus for 2D
        }

        break;
    }
    case stk::topology::ELEMENT_RANK: // for node sets
    {
        // Check if sets were defined by the user
        if ( mSetRankFlags[2])
        {
            aSetName = ( mSetNames[2] )[aSetId-1];
        }
        else
        {
            aSetName = "block_" + std::to_string( aSetId
                    );
        }

        break;
    }
    default:
    {
        MORIS_ASSERT( 0, "Invalid rank provided in get_entities_in_set." );
    }
    }

    return this->get_set_entity_ids( aRequiredEntityRank, aSetName );
}

// ----------------------------------------------------------------------------

// Access set entity ids
moris::Mat < uint >
moris::STK_Implementation::get_set_entity_ids(
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
    Mat< uint > tOutputEntityIds ( tNumEntities, 1 );
    for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
    {
        tOutputEntityIds( iEntity ) = ( uint ) mMtkMeshBulkData->identifier( aEntities[iEntity] );
    }

    return tOutputEntityIds;
}

// ----------------------------------------------------------------------------

// Access field entities of selected portion of the mesh
moris::Mat < uint >
moris::STK_Implementation::get_intersected_entities_field_set(
        enum EntityRank  aEntityRank,
        std::string      aFieldName,
        std::string      aSetName ) const
{
    // Get pointer to field and set defined by input names
    stk::mesh::Part* const tSetPart = mMtkMeshMetaData->get_part( aSetName );
    MORIS_ASSERT( tSetPart != NULL, "Set not found. Double check name provided." );

    // Get pointer to field defined by input name
    stk::mesh::Part* const tFieldPart = mMtkMeshMetaData->get_part( aFieldName );
    MORIS_ASSERT( tFieldPart != NULL, "Field not found. Double check name provided." );
    
    // Access data through a selector
    stk::mesh::Selector tSetSelector = *tSetPart & *tFieldPart;
    stk::mesh::EntityVector aEntities;
    stk::mesh::get_selected_entities( tSetSelector, mMtkMeshBulkData->buckets( this->get_stk_entity_rank( aEntityRank ) ), aEntities );
    
    // Get entity Ids
    uint tNumEntities = aEntities.size();
    Mat< uint > tOutputEntityIds ( tNumEntities, 1);
    for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
    {
        tOutputEntityIds( iEntity ) = ( uint ) mMtkMeshBulkData->identifier( aEntities[iEntity] );
    }
    
    return tOutputEntityIds;
}

// ----------------------------------------------------------------------------

// Access field data of selected portion of the mesh
moris::Mat < moris::real >
moris::STK_Implementation::get_intersected_data_field_set(
        enum EntityRank  aEntityRank,
        std::string      aFieldName,
        std::string      aSetName ) const
{
    stk::mesh::EntityRank aStkFieldRank  = this->get_stk_entity_rank( aEntityRank );

    // Get pointer to field and set defined by input names
    stk::mesh::Part* const tSetPart = mMtkMeshMetaData->get_part( aSetName );
    MORIS_ASSERT( tSetPart != NULL, "Set not found. Double check name provided." );

    // Get pointer to field defined by input name
    stk::mesh::FieldBase const* aField_base = mMtkMeshMetaData->get_field( aStkFieldRank, aFieldName );
    MORIS_ASSERT( aField_base != NULL, "Field not found. Double check name provided." );

    MORIS_ASSERT( this->get_mtk_entity_rank( aField_base->entity_rank() ) == aEntityRank, "No field with name and rank provided was found." );

    // Access data through a selector
    stk::mesh::Selector tSetSelector = *tSetPart & *aField_base;
    stk::mesh::EntityVector aEntities;
    stk::mesh::get_selected_entities( tSetSelector, mMtkMeshBulkData->buckets( aStkFieldRank ), aEntities );

    // Get field variables
    uint tNumFieldComp = stk::mesh::field_scalars_per_entity( *aField_base, aEntities[0] );

    // Get entity Ids
    uint tNumFieldEntities = aEntities.size();
    Mat< real > tOutputField( tNumFieldEntities, tNumFieldComp, -1.0 );

    // Loop over number of entities
    for ( uint iEntity = 0; iEntity < tNumFieldEntities; ++iEntity )
    {
        // Store the coordinates of the current node
        if ( mMtkMeshBulkData->is_valid( aEntities[iEntity] ) )
        {
            double* tEntityData = static_cast <double*> ( stk::mesh::field_data( *aField_base, aEntities[iEntity] ) );

            for ( uint iComp = 0; iComp < tNumFieldComp; ++iComp )
            {
                tOutputField( iEntity, iComp ) =  tEntityData[iComp];
            }
        }
    }

    return tOutputField;
}

// ----------------------------------------------------------------------------

// Get local ids (only for debugging)
moris::Mat< uint >
moris::STK_Implementation::get_entity_local_ids_connected_to_entity(
        uint        const aEntityIndex,
        EntityRank   const aInputEntityRank,
        EntityRank   const aOutputEntityRank ) const
{
    // Call function that gets the connected entities
    uint aEntityId = { mLocalToGlobalNodeMap( aEntityIndex ) };

    if ( aInputEntityRank == EntityRank::ELEMENT )
    {
        aEntityId = {mLocalToGlobalElemMap( aEntityIndex )};
    }
    else
    {
        MORIS_ASSERT( 0, "procs_who_share_entity_interface only processes nodal and elemental indices for the moment." );
    }

    Mat< uint > aEntitiesConnected = this->entities_connected_to_given_entity( aEntityId, aInputEntityRank, aOutputEntityRank );

    // Declare moris::Mat that will contain the local entities
    uint tNumOutputEntities = aEntitiesConnected.length();
    Mat< uint > tLocalIds( 1, tNumOutputEntities );

    // Get stk entity rank form MORIS enum
    stk::mesh::EntityRank tOutputRank = this->get_stk_entity_rank( aOutputEntityRank );

    // Fill local ids to moris::Mat
    for ( uint iEntity = 0; iEntity < tNumOutputEntities; ++iEntity )
    {
        tLocalIds( iEntity ) = mMtkMeshBulkData->local_id( mMtkMeshBulkData->get_entity( tOutputRank, aEntitiesConnected( iEntity ) ) );
    }

    return tLocalIds;
}

// ----------------------------------------------------------------------------

// Function to access node coordinates for any given set of node ids
moris::Mat< moris::real >
moris::STK_Implementation::get_nodes_coords_generic(
        Mat< uint >   aNodeIds ) const
{
    // Get number of spatial dimensions and entities provided
    uint tNumIds = aNodeIds.length();

    // Initialize output matrix and coordinate field
    Mat< real > tSelectedNodesCoords( tNumIds, mNumDims );
    stk::mesh::FieldBase const* aCoordField = mMtkMeshMetaData->coordinate_field();

    // Loop over all nodes to get their coordinates
    for ( uint iNode=0; iNode < tNumIds; ++iNode )
    {
        // Declare node entity
        stk::mesh::Entity aNodeEnitity = mMtkMeshBulkData->get_entity( stk::topology::NODE_RANK, aNodeIds( iNode ) );

        // Get node coordinates and store them into output matrix
        double* tCoords = static_cast<double *>( stk::mesh::field_data( *aCoordField, aNodeEnitity ) );
        for ( uint iDim = 0; iDim < mNumDims ; ++iDim )
        {
            tSelectedNodesCoords( iNode, iDim ) = tCoords[iDim];
        }
    }

    return tSelectedNodesCoords;
}

// ----------------------------------------------------------------------------

// Function to access node coordinates for any given set of node ids
moris::Mat< moris::real >
moris::STK_Implementation::get_selected_nodes_coords_lcl_ind(
        Mat< uint >   aNodeInds ) const
{
    Mat< uint > aNodeIds( aNodeInds.length(), 1 );
    for ( uint iNode=0; iNode < aNodeInds.length(); ++iNode )
    {
        aNodeIds( iNode ) = mLocalToGlobalNodeMap( aNodeInds( iNode ) );
    }

    return this->get_nodes_coords_generic( aNodeIds );
}

// ----------------------------------------------------------------------------

// Ask STK to provide new unique ids
moris::Mat< uint >
moris::STK_Implementation::generate_unique_entity_ids(
        uint                    aNumEntities,
        stk::mesh::EntityRank   aNeedConnectivityOfType ) const
{
    std::vector<stk::mesh::EntityId> tAvailableNodeIDs; // generate_new_ids requires a variable of this type
    mMtkMeshBulkData->generate_new_ids( aNeedConnectivityOfType, aNumEntities, tAvailableNodeIDs );

    moris::Mat<moris::uint> aAvailableNodeIDs ( aNumEntities, 1 );

    for ( uint i = 0; i < aNumEntities; i++ )
    {
        aAvailableNodeIDs( i, 0 ) = tAvailableNodeIDs[i];
    }

    return aAvailableNodeIDs;
}

// ----------------------------------------------------------------------------

// Dummy function intended to be general in a future commit.
void
moris::STK_Implementation::add_shared_field_interface( ) const
{
    stk::mesh::Field<double>*  universal_scalar_edge_field  =
            stk::mesh::get_field_by_name<stk::mesh::Field<double> >( "universal_scalar_edge_field", *mMtkMeshMetaData );
    stk::mesh::Selector shared_sel = mMtkMeshMetaData->globally_shared_part();
    stk::mesh::BucketVector const& shared_entity_buckets = mMtkMeshBulkData->get_buckets( stk::topology::EDGE_RANK, shared_sel );

    uint p_rank = this->get_parallel_rank();
    for ( size_t b = 0, be = shared_entity_buckets.size(); b < be; ++b )
    {
        stk::mesh::Bucket& bucket = *shared_entity_buckets[b];
        for ( size_t e = 0, ee = bucket.size(); e < ee; ++e )
        {
            stk::mesh::Entity edge = bucket[e];
            stk::mesh::EntityId edge_id = mMtkMeshBulkData->identifier( edge );

            double* fieldPtr = stk::mesh::field_data(*universal_scalar_edge_field, edge );
            *fieldPtr = p_rank + 10.0 + edge_id;
        }
    }

    stk::mesh::Selector locally_owned_sel             = mMtkMeshMetaData->locally_owned_part();
    stk::mesh::BucketVector const& owned_edge_buckets = mMtkMeshBulkData->get_buckets( stk::topology::EDGE_RANK, locally_owned_sel );

    for ( size_t b = 0, be = owned_edge_buckets.size(); b < be; ++b )
    {
        stk::mesh::Bucket& bucket = *owned_edge_buckets[b];
        for ( size_t e = 0, ee = bucket.size(); e < ee; ++e )
        {
            stk::mesh::Entity edge      = bucket[e];
            stk::mesh::EntityId edge_id = mMtkMeshBulkData->identifier( edge );

            double* fieldPtr = stk::mesh::field_data( *universal_scalar_edge_field, edge );
            *fieldPtr = p_rank + 10.0 + edge_id;
        }
    }
}

// ----------------------------------------------------------------------------

// Translate ranks from MTK to STK format
stk::mesh::EntityRank
moris::STK_Implementation::get_stk_entity_rank(
        EntityRank  aEntityRank ) const
{
    if ( aEntityRank == EntityRank::NODE )
    {
        return stk::topology::NODE_RANK;
    }
    else if ( aEntityRank == EntityRank::EDGE )
    {
        return stk::topology::EDGE_RANK;
    }
    else if ( aEntityRank == EntityRank::FACE )
    {
        return stk::topology::FACE_RANK;
    }
    else if ( aEntityRank == EntityRank::ELEMENT )
    {
        return stk::topology::ELEMENT_RANK;
    }
    else
    {
        return stk::topology::INVALID_RANK;
    }
}

// ----------------------------------------------------------------------------

// Translate ranks from STK to MTK format
EntityRank
moris::STK_Implementation::get_mtk_entity_rank(
        stk::mesh::EntityRank  aEntityRank ) const
{
    if ( aEntityRank == stk::topology::NODE_RANK )
    {
        return EntityRank::NODE;
    }
    else if ( aEntityRank == stk::topology::EDGE_RANK )
    {
        return  EntityRank::EDGE;
    }
    else if ( aEntityRank == stk::topology::FACE_RANK )
    {
        return EntityRank::FACE;
    }
    else if ( aEntityRank == stk::topology::ELEMENT_RANK )
    {
        return EntityRank::ELEMENT;
    }
    else
    {
        return EntityRank::INVALID;
    }
}

// ----------------------------------------------------------------------------

// Get the processor that owns a given entity
uint
moris::STK_Implementation::parallel_owner_rank_by_entity_id(
        uint              aEntityId,
        enum EntityRank   aEntityRank ) const
{
    // Check if function is called in a serial run. If so, no parallel information is available.
    uint tProcSize = par_size();
    if ( tProcSize == 1 )
    {
        return UINT_MAX;
    }

    //Get entity
    stk::mesh::Entity aEntity = mMtkMeshBulkData->get_entity( this->get_stk_entity_rank( aEntityRank ), aEntityId );

    //processor rank that owns entity
    uint tProc = mMtkMeshBulkData->parallel_owner_rank( aEntity );

    return tProc;
}

// ----------------------------------------------------------------------------

// Get the processor that owns a given entity
uint
moris::STK_Implementation::parallel_owner_rank_by_entity_index(
        uint              aEntityIndex,
        enum EntityRank   aEntityRank ) const
{
    // Convert index to ID
    uint tEntityId = 0;
    if ( aEntityRank == EntityRank::NODE )
    {
        tEntityId = { mLocalToGlobalNodeMap( aEntityIndex ) };
    }
    else if ( aEntityRank == EntityRank::ELEMENT )
    {
        tEntityId = { mLocalToGlobalElemMap( aEntityIndex ) };
    }
    else if ( aEntityRank == EntityRank::EDGE )
    {
        tEntityId = { mLocalToGlobalEdgeMap( aEntityIndex ) };
    }
    else if ( aEntityRank == EntityRank::FACE )
    {
        tEntityId = { mLocalToGlobalFaceMap( aEntityIndex ) };
    }

    return this->parallel_owner_rank_by_entity_id( tEntityId, aEntityRank );
}

// ----------------------------------------------------------------------------

// return processors sharing a particular entity
moris::Mat< uint >
moris::STK_Implementation::get_procs_sharing_entity_by_id(
        uint              aEntityID,
        enum EntityRank   aEntityRank ) const
{
    // Initialize returning mat with UINT_MAX in case no processors are sharing the given entity
    Mat< uint > tSharedProcsMat( 1, 1, UINT_MAX );

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

// return processors sharing a particular entity
moris::Mat< uint >
moris::STK_Implementation::get_procs_sharing_entity_by_index(
        uint              aEntityIndex,
        enum EntityRank   aEntityRank ) const
{
    // Initialize returning mat with UINT_MAX in case no processors are sharing the given entity
    Mat< uint > tSharedProcsMat ( 1, 1, UINT_MAX );

    // Check if function is called in a serial run. If so, no parallel information is available.
    uint tProcSize = par_size();
    if ( tProcSize == 1 )
    {
        return tSharedProcsMat;
    }

    // Convert index to ID
    uint tEntityId = 0;
    if ( aEntityRank == EntityRank::NODE )
    {
        tEntityId = { mLocalToGlobalNodeMap( aEntityIndex ) };
    }
    else if ( aEntityRank == EntityRank::ELEMENT )
    {
        tEntityId = { mLocalToGlobalElemMap( aEntityIndex ) };
    }
    else if ( aEntityRank == EntityRank::EDGE )
    {
        tEntityId = { mLocalToGlobalEdgeMap( aEntityIndex ) };
    }
    else if ( aEntityRank == EntityRank::FACE )
    {
        tEntityId = { mLocalToGlobalFaceMap( aEntityIndex ) };
    }

    return this->get_procs_sharing_entity_by_id( tEntityId, aEntityRank );
}

// ----------------------------------------------------------------------------

// Get node ids in local map
moris::Mat < uint >
moris::STK_Implementation::get_node_ids_from_local_map(
        Mat< uint >   aLocalInds ) const
{
    uint tNumIdsRequested = aLocalInds.length();
    Mat < uint > tRequestedIds;

    for ( uint iNode = 0; iNode < tNumIdsRequested; ++iNode )
    {
        tRequestedIds( iNode ) = mLocalToGlobalNodeMap( aLocalInds( iNode ) );
    }

    return tRequestedIds;
}

// ----------------------------------------------------------------------------

// Get element ids from local map
moris::Mat < uint >
moris::STK_Implementation::get_element_ids_from_local_map_interface(
        Mat< uint >   aLocalInds ) const
{
    uint tNumIdsRequested = aLocalInds.length();
    Mat < uint > tRequestedIds;

    for ( uint iElem = 0; iElem < tNumIdsRequested; ++iElem )
    {
        tRequestedIds( iElem ) = mLocalToGlobalElemMap( aLocalInds( iElem ) );
    }

    return tRequestedIds;
}

// ----------------------------------------------------------------------------

// Access field data
moris::Mat < moris::real >
moris::STK_Implementation::get_field_values(
        enum EntityRank   aEntityRank,
        std::string       aFieldName )
{
    stk::mesh::EntityRank aStkFieldRank  = this->get_stk_entity_rank( aEntityRank );

    // Get pointer to field defined by input name
    stk::mesh::FieldBase const* aField_base = mMtkMeshMetaData->get_field( aStkFieldRank, aFieldName );
    MORIS_ASSERT( aField_base != NULL, "Field not found. Double check name provided." );

    MORIS_ASSERT( this->get_mtk_entity_rank( aField_base->entity_rank() ) == aEntityRank, "No field with name and rank provided was found." );

    // Get field variables
    Mat< uint > aEntityIds = this->get_field_entities( aEntityRank, aFieldName );
    stk::mesh::Entity aDummEntity = mMtkMeshBulkData->get_entity( aStkFieldRank, aEntityIds(0,0) );
    uint tNumFieldComp = stk::mesh::field_scalars_per_entity( *aField_base, aDummEntity );
    uint tNumFieldEntities = aEntityIds.length();

    // Declare output moris::mat
    Mat< real > tOutputField( tNumFieldEntities, tNumFieldComp, -1.0 );
    // Loop over number of entities
    for ( uint iEntity = 0; iEntity < tNumFieldEntities; ++iEntity )
    {
        // Get current entity from global Id
        stk::mesh::Entity aEntity = mMtkMeshBulkData->get_entity( aStkFieldRank, aEntityIds( iEntity ) );

        // Store the coordinates of the current node
        if ( mMtkMeshBulkData->is_valid( aEntity ) )
        {
            double* tEntityData = static_cast <double*> ( stk::mesh::field_data( *aField_base, aEntity ) );

            for ( uint iComp = 0; iComp < tNumFieldComp; ++iComp )
            {
                tOutputField( iEntity, iComp ) =  tEntityData[iComp];
            }
        }
    }

    return tOutputField;
}

// ----------------------------------------------------------------------------

// Access field entities
moris::Mat < uint >
moris::STK_Implementation::get_field_entities(
        stk::mesh::FieldBase const* aField )
{
    // Get selector with field and get entities
    stk::mesh::Selector aSelectField = stk::mesh::selectField( *aField ) &
            ( mMtkMeshMetaData->globally_shared_part() | mMtkMeshMetaData->locally_owned_part() );
    std::vector< stk::mesh::Entity > aEntities;
    stk::mesh::get_selected_entities( aSelectField, mMtkMeshBulkData->buckets( aField->entity_rank() ), aEntities );

    // Interface stk and moris variables
    uint tNumEntities = aEntities.size();
    Mat< uint > tOutputEntities( tNumEntities, 1 );

    for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
    {
        tOutputEntities( iEntity ) =  mMtkMeshBulkData->identifier( aEntities[iEntity] );
    }

    return tOutputEntities;
}

// ----------------------------------------------------------------------------

// Access entities of given portion of the mesh
moris::Mat < uint >
moris::STK_Implementation::get_part_entities(
        enum EntityRank   aEntityRank,
        std::string       aPartName )
{
    stk::mesh::Part &aPart = *mMtkMeshMetaData->get_part( aPartName );
    stk::mesh::Selector aPartSelector( aPart );
    std::vector< stk::mesh::Entity > aEntities;
    stk::mesh::get_selected_entities( aPartSelector, mMtkMeshBulkData->buckets( this->get_stk_entity_rank( aEntityRank ) ), aEntities );

    // Interface stk and moris variables
    uint tNumEntities = aEntities.size();
    Mat< uint > tOutputEntities( tNumEntities, 1 );
    for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
    {
        tOutputEntities( iEntity ) = mMtkMeshBulkData->identifier( aEntities[iEntity] );
    }

    return tOutputEntities;
}

// ----------------------------------------------------------------------------

// Print information to screen (for debugging purposes only)
void
moris::STK_Implementation::print_parts_to_screen()
{
    // Print parts (only for debugging purposes )
    const stk::mesh::PartVector & all_parts = mMtkMeshMetaData->get_mesh_parts();

    for ( stk::mesh::PartVector::const_iterator ip = all_parts.begin(); ip != all_parts.end(); ++ip )
    {
        stk::mesh::Part * const part = *ip;
        std::string PartName         = part->name();

        std::cout<<"\nPart name: "<<PartName<<std::endl;
    }
}

// ----------------------------------------------------------------------------

// Print information to screen (for debugging purposes only)
void
moris::STK_Implementation::print_fields_to_screen()
{
    // Print fields (only for debugging purposes )
    const std::vector<stk::mesh::FieldBase*>&fields = mMtkMeshMetaData->get_fields();

    std::vector<stk::mesh::FieldBase *>::const_iterator fieldIterator = fields.begin();
    for ( ; fieldIterator != fields.end(); ++fieldIterator )
    {
        std::string tIterFieldName = ( *fieldIterator )->name();
        std::cout<<"\nField name: "<<tIterFieldName<<std::endl;
    }
}

// ----------------------------------------------------------------------------

// This function computes a global coordinate for a node created on a parent entity
// give a local coordinate.First it checks to see if the local coordinates provided
// are sufficient for the parent entity rank. Depending on the number of local coordinates
// an appropriate interpolation scheme is selected.
moris::Mat< moris::real >
moris::STK_Implementation::interpolate_to_location_on_entity(
        enum EntityRank  aParentEntityRank,
        uint             aParentEntityIndex,
        Mat< real >        aLclCoord )
{
    MORIS_ASSERT( ( uint ) aParentEntityRank == aLclCoord.n_cols(), "The number of local coordinates provided need to match the dimension of parent entity -1");

    Mat< uint > tNodes      = this->get_entity_local_ids_connected_to_entity( aParentEntityIndex, aParentEntityRank, EntityRank::NODE );
    Mat< real > tNodeCoords = this->get_selected_nodes_coords_lcl_ind( tNodes );

    if ( aParentEntityRank == EntityRank::ELEMENT )
    {
        return Interpolation::trilinear_interpolation( tNodeCoords, aLclCoord );
    }
    else if ( aParentEntityRank == EntityRank::FACE )
    {
        return Interpolation::bilinear_interpolation( tNodeCoords, aLclCoord );
    }
    else if ( aParentEntityRank == EntityRank::EDGE )
    {
        return Interpolation::linear_interpolation_location( tNodeCoords, aLclCoord );
    }
    else
    {
        MORIS_ASSERT( 0, "Interpolation must occur on an element, face or edge." );
        Mat< real > dummy;
        return dummy;
    }
}

// ----------------------------------------------------------------------------

// Check for duplicate coordinates
moris::Mat< uint >
moris::STK_Implementation::duplicate_node_coord_check()
{
    Mat< real > tCoord = get_all_nodes_coords();

    // Check for duplicates in the coordinate list
    Mat< uint > duplicate_list = Debug::duplicate_row_check( tCoord );

    return duplicate_list;
}

// ----------------------------------------------------------------------------

// Check for duplicate coordinates
moris::Mat< uint >
moris::STK_Implementation::duplicate_node_coord_check(
        moris::Mat< real >& aCoord )
{
    // Check for duplicates in the coordinate list
    Mat< uint > duplicate_list = Debug::duplicate_row_check( aCoord );

    return duplicate_list;
}

// ----------------------------------------------------------------------------

// Check for duplicate coordinates and ids
moris::Mat< uint >
moris::STK_Implementation::duplicate_node_coord_and_id_check(
        Mat< real >& aCoord,
        Mat< uint >& aId )
{
    // Check for duplicates in the coordinate and id list
    Mat< uint > duplicate_coord_list = Debug::duplicate_row_check( aCoord );
    Mat< uint > duplicate_id_list    = Debug::duplicate_row_check( aId );
    Mat< uint > duplicate_list       = Debug::duplicate_row_check( duplicate_coord_list,duplicate_id_list );

    return duplicate_list;
}

// ----------------------------------------------------------------------------

// Check for duplicate coordinates and ids from a Cell through MPI
moris::Mat< uint >
moris::STK_Implementation::duplicate_node_coord_and_id_check(
        Cell< Mat< real > >& aCoord,
        Cell< Mat< uint > >& aId )
{
    uint tNumCell   = aCoord.size();
    uint tCoordrows = 0;
    uint tIdrows    = 0;
    Mat< real > aAllCoord( aCoord( 0 ).n_rows(), aCoord( 0 ).n_cols(), UINT_MAX );
    Mat< uint > aAllId( aId( 0 ).n_rows(), aId( 0 ).n_cols(), UINT_MAX );

    for ( uint i = 0; i<tNumCell; ++i )
    {
        Mat< real > tCoord = aCoord( i );
        Mat< uint > tId    = aId( i );
        tCoordrows        += tCoord.n_rows();
        tIdrows           += tId.n_rows();

        aAllCoord.resize( tCoordrows, tCoord.n_cols() );
        aAllId.resize( tIdrows,tId.n_cols() );
        aAllCoord.rows( tCoordrows - tCoord.n_rows(), tCoordrows - 1 ) = tCoord.rows( 0, tCoord.n_rows() - 1 );
        aAllId.rows( tIdrows - tId.n_rows(), tIdrows - 1 )             = tId.rows( 0, tId.n_rows() - 1 );
    }

    // Check for duplicates in the coordinate and id list
    Mat< uint > duplicate_coord_list = Debug::duplicate_row_check( aAllCoord );
    Mat< uint > duplicate_id_list    = Debug::duplicate_row_check( aAllId );
    Mat< uint > duplicate_list       = Debug::duplicate_row_check( duplicate_coord_list, duplicate_id_list );

    return duplicate_list;
}

// ----------------------------------------------------------------------------

// Check for duplicate coordinates and ids and return only those, where coordinate and id do not belong together
moris::Mat< uint >
moris::STK_Implementation::duplicate_node_coord_and_id_check_problems(
        Mat< real >& aCoord,
        Mat< uint >& aId )
{
    // Check for duplicates in the coordinate and id list
    Mat< uint > duplicate_coord_list    = Debug::duplicate_row_check( aCoord );
    Mat< uint > duplicate_id_list       = Debug::duplicate_row_check( aId );
    Mat< uint > duplicate_coord_id_list = Debug::duplicate_row_check( duplicate_coord_list,duplicate_id_list );
    Mat< uint > problem_list            = Debug::duplicate_row_check_problems( duplicate_coord_list,duplicate_coord_id_list );

    return problem_list;
}

// ----------------------------------------------------------------------------

// Check for duplicate coordinates and ids from a Cell through MPI and return only those, where coordinate and id do not belong together
moris::Mat< uint >
moris::STK_Implementation::duplicate_node_coord_and_id_check_problems(
        Cell< Mat< real > >& aCoord,
        Cell< Mat< uint > >& aId )
{
    uint tNumCell   = aCoord.size();
    uint tCoordrows = 0;
    uint tIdrows    = 0;
    Mat< real > aAllCoord( aCoord( 0 ).n_rows(), aCoord( 0 ).n_cols(), UINT_MAX );
    Mat< uint > aAllId( aId( 0 ).n_rows(), aId( 0 ).n_cols(), UINT_MAX );

    for ( uint i = 0; i<tNumCell; i++)
    {
        Mat< real > tCoord = aCoord( i );
        Mat< uint > tId    = aId( i );
        tCoordrows        += tCoord.n_rows();
        tIdrows           += tId.n_rows();
        aAllCoord.resize( tCoordrows, tCoord.n_cols() );
        aAllId.resize( tIdrows, tId.n_cols() );
        aAllCoord.rows( tCoordrows - tCoord.n_rows(), tCoordrows - 1 ) = tCoord.rows( 0, tCoord.n_rows() - 1 );
        aAllId.rows( tIdrows - tId.n_rows(), tIdrows - 1 ) = tId.rows( 0, tId.n_rows() - 1 );
    }

    // Check for duplicates in the coordinate and id list
    Mat< uint > duplicate_coord_list    = Debug::duplicate_row_check( aAllCoord );
    Mat< uint > duplicate_id_list       = Debug::duplicate_row_check( aAllId );
    Mat< uint > duplicate_coord_id_list = Debug::duplicate_row_check( duplicate_coord_list,duplicate_id_list );
    Mat< uint > problem_list            = Debug::duplicate_row_check_problems( duplicate_coord_list, duplicate_coord_id_list );

    return problem_list;
}

// -----------------------------------------------------------------

// Function to create communication lists in parallel
void
moris::STK_Implementation::create_communication_lists_from_string()
{
    // Get basic mesh information
    uint tNumElems = this->get_num_elems();
    uint tNumNodes = this->get_num_nodes();

    // Access entities stored in mesh database
    Mat< uint > tElemIds = this->get_entities_universal( EntityRank::ELEMENT );
    Mat< uint > tNodeIds = this->get_entities_universal( EntityRank::NODE );

    // resize member variable to its right size
    mLocalToGlobalNodeMap.resize( tNumNodes, 1 );
    mNodeMapToOwnerProc.  resize( tNumNodes, 1 );

    // Populate internal member variable that contains the local index to
    // global id node communication information
    for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
    {
        // local to global and owner processor
        mLocalToGlobalNodeMap( iNode ) = tNodeIds( iNode );
        mNodeMapToOwnerProc( iNode )   = this->parallel_owner_rank_by_entity_index( iNode, EntityRank::NODE );
    }

    // resize member variable to its right size
    mLocalToGlobalElemMap.resize( tNumElems, 1 );
    mElemMapToOwnerProc.  resize( tNumElems, 1 );

    // Populate internal member variable that contains the local index to
    // global id element communication information
    for ( uint iElem = 0; iElem < tNumElems; ++iElem )
    {
        // local to global and owner processor
        mLocalToGlobalElemMap( iElem ) = tNodeIds( iElem );
        mElemMapToOwnerProc( iElem )   = this->parallel_owner_rank_by_entity_index( iElem, EntityRank::ELEMENT );
    }

    // Create maps for edges if the problem is not 1D, and maps for faces if it is 3D.
    // NOTE: Not supporting 1D elements in 2D nor 3D space.
    if ( mNumDims > 1 )
    {
        uint tNumEdges = this->get_num_edges();

        // Access entities stored in mesh database
        Mat< uint > tEdgeIds = this->get_entities_universal( EntityRank::EDGE );

        // resize member variable to its right size
        mLocalToGlobalEdgeMap.resize( tNumEdges, 1 );
        mEdgeMapToOwnerProc.resize( tNumEdges, 1 );

        // Populate internal member variable that contains the local index to
        // global id node communication information
        for ( uint iEdge = 0; iEdge < tNumEdges; ++iEdge )
        {
            // local to global and owner processor
            mLocalToGlobalEdgeMap( iEdge ) = tEdgeIds( iEdge );
            mEdgeMapToOwnerProc( iEdge )   = this->parallel_owner_rank_by_entity_index( iEdge, EntityRank::EDGE );
        }
    }

    if ( mNumDims > 2 )
    {
        uint tNumFaces   = this->get_num_faces();

        // Access entities stored in mesh database
        Mat< uint > tFaceIds = this->get_entities_universal( EntityRank::FACE );

        // resize member variable to its right size
        mLocalToGlobalFaceMap.resize( tNumFaces, 1 );
        mFaceMapToOwnerProc.resize( tNumFaces, 1 );

        // Populate internal member variable that contains the local index to
        // global id element communication information
        for ( uint iFace = 0; iFace < tNumFaces; ++iFace )
        {
            // local to global and owner processor
            mLocalToGlobalFaceMap( iFace ) = tFaceIds( iFace );
            mFaceMapToOwnerProc( iFace )   = this->parallel_owner_rank_by_entity_index( iFace, EntityRank::FACE );
        }
    }

    this->create_shared_communication_lists();
}

// -----------------------------------------------------------------

// Function to create edges and faces communication lists in parallel for meshes generated from data
void
moris::STK_Implementation::create_additional_communication_lists_from_data()
{
    this->create_facets_communication_lists();
    this->create_owners_communication_lists();
    this->create_shared_communication_lists();
}

// -----------------------------------------------------------------

// Function to create edges and faces communication lists in parallel
void
moris::STK_Implementation::create_facets_communication_lists()
{
    // Create maps for edges if the problem is not 1D, and maps for faces if it is 3D.
    // NOTE: Not supporting 1D elements in 2D nor 3D space.
    if ( mNumDims > 1 )
    {
        uint tNumEdges = this->get_num_edges();

        // Access entities stored in mesh database
        Mat< uint > tEdgeIds = this->get_entities_universal( EntityRank::EDGE );

        // resize member variable to its right size
        mLocalToGlobalEdgeMap.resize( tNumEdges, 1 );
        mEdgeMapToOwnerProc.resize( tNumEdges, 1 );

        // Populate internal member variable that contains the local index to
        // global id node communication information
        for ( uint iEdge = 0; iEdge < tNumEdges; ++iEdge )
        {
            // local to global and owner processor
            mLocalToGlobalEdgeMap( iEdge ) = tEdgeIds( iEdge );
            mEdgeMapToOwnerProc( iEdge )   = this->parallel_owner_rank_by_entity_index( iEdge, EntityRank::EDGE );
        }
    }

    if ( mNumDims > 2 )
    {
        uint tNumFaces   = this->get_num_faces();

        // Access entities stored in mesh database
        Mat< uint > tFaceIds = this->get_entities_universal( EntityRank::FACE );

        // resize member variable to its right size
        mLocalToGlobalFaceMap.resize( tNumFaces, 1 );
        mFaceMapToOwnerProc.resize( tNumFaces, 1 );

        // Populate internal member variable that contains the local index to
        // global id element communication information
        for ( uint iFace = 0; iFace < tNumFaces; ++iFace )
        {
            // local to global and owner processor
            mLocalToGlobalFaceMap( iFace ) = tFaceIds( iFace );
            mFaceMapToOwnerProc( iFace )   = this->parallel_owner_rank_by_entity_index( iFace, EntityRank::FACE );
        }
    }
}

// -----------------------------------------------------------------

// Function to create edges and faces communication lists in parallel
void
moris::STK_Implementation::create_owners_communication_lists()
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
        mNodeMapToOwnerProc( iNode ) = this->parallel_owner_rank_by_entity_index( iNode, EntityRank::NODE );
    }

    // resize member variable to its right size
    mElemMapToOwnerProc.resize( tNumElems, 1 );

    // Populate internal member variable that contains the local index to
    // global id element communication information
    for ( uint iElem = 0; iElem < tNumElems; ++iElem )
    {
        // Owner processor
        mElemMapToOwnerProc( iElem ) = this->parallel_owner_rank_by_entity_index( iElem, EntityRank::ELEMENT );
    }
}

// -----------------------------------------------------------------

// Function to create edges and faces communication lists in parallel
void
moris::STK_Implementation::create_shared_communication_lists()
{
    // Get basic mesh information
    Mat < uint > tNodesShared = this->get_entities_glb_shared_current_proc( EntityRank::NODE );
    uint tNumNodesShared      = tNodesShared.length();

    // Generate list of processors sharing information
    // -----------------------------------------------

    std::vector < uint > tActiveSharedProcs;

    // Loop over the number of nodes shared to get the shared processors
    for( uint iNodeShared = 0; iNodeShared < tNumNodesShared; ++iNodeShared )
    {
        Mat < uint > tProcsSharing = this->get_procs_sharing_entity_by_id( tNodesShared( iNodeShared ), EntityRank::NODE );

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

// -----------------------------------------------------------------

// Function to provide entities owned and shared based on rank
moris::Mat< uint >
moris::STK_Implementation::get_entities_owned_and_shared_by_current_proc( EntityRank   aEntityRank ) const
{
    if ( aEntityRank == EntityRank::NODE )
    {
        return mLocalToGlobalNodeMap;
    }
    else if ( aEntityRank == EntityRank::EDGE )
    {
        return mLocalToGlobalEdgeMap;
    }
    else if ( aEntityRank == EntityRank::FACE )
    {
        return mLocalToGlobalFaceMap;
    }
    else if ( aEntityRank == EntityRank::ELEMENT )
    {
        return mLocalToGlobalElemMap;
    }
    else
    {
        MORIS_ASSERT( 0, "Invalid rank provided in get_entities_owned_and_shared_by_current_proc." );
    }

    Mat< uint > tDummyConn( 1, 1, UINT_MAX );
    return tDummyConn;
}

// -----------------------------------------------------------------

// Function to create edges and faces communication lists in parallel
moris::Cell < moris::Cell < uint > >
moris::STK_Implementation::get_shared_info_by_entity( uint aNumActiveSharedProcs, enum EntityRank  aEntityRank )
{
    // Generate elements shared per processor list
    // -------------------------------------------
    Mat < uint > tEntitiesShared = this->get_entities_glb_shared_current_proc( aEntityRank );
    uint tNumEntitiesShared      = tEntitiesShared.length();
    uint tParallelRank           = this->get_parallel_rank();

    Cell < Cell < uint > > tTemporaryEntityMapSharingProcs( aNumActiveSharedProcs );

    // Loop over the number of nodes shared to get the shared processors
    for( uint iElemShared = 0; iElemShared < tNumEntitiesShared; ++iElemShared )
    {
        Mat < uint > tProcsSharing = this->get_procs_sharing_entity_by_id( tEntitiesShared( iElemShared ), aEntityRank );
        for( uint iProc = 0; iProc < tProcsSharing.length(); ++iProc )
        {
            if( tProcsSharing(iProc) != tParallelRank )
            {
                tTemporaryEntityMapSharingProcs( mProcsSharedToIndex[ tProcsSharing( iProc ) ] ).push_back( tEntitiesShared( iElemShared ) );
            }
        }
    }

    return tTemporaryEntityMapSharingProcs;
}
