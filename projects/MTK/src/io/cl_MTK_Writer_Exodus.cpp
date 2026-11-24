/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Writer_Exodus.cpp
 *
 */

#include <exodusII.h>

#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"

#include <iostream>
#include <utility>
#include <filesystem>

namespace moris::mtk
{
    //--------------------------------------------------------------------------------------------------------------
    // Public
    //--------------------------------------------------------------------------------------------------------------

    Writer_Exodus::Writer_Exodus( mtk::Mesh* aMeshPointer )
            : mMesh( aMeshPointer )
    {
        this->set_error_options( true, true, true );
    }

    //--------------------------------------------------------------------------------------------------------------

    Writer_Exodus::Writer_Exodus()
    {
        mMesh = nullptr;
        this->set_error_options( true, true, true );
    }

    //--------------------------------------------------------------------------------------------------------------

    Writer_Exodus::~Writer_Exodus()
    {
        if ( mExoID >= 0 )
        {
            ex_close( mExoID );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::set_error_options(
            bool abort,
            bool debug,
            bool verbose )
    {
        ex_opts( abort * EX_ABORT | debug * EX_DEBUG | verbose * EX_VERBOSE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::open_file(
            std::string& aExodusFileName,
            bool         aReadOnly,
            float        aVersion )
    {
        MORIS_ERROR( mExoID == -1, "Exodus file is currently open, call close_file() before opening a new one." );

        int tCPUWordSize = sizeof( real ), tIOWordSize = 0;

        if ( aReadOnly )
        {
            mExoID = ex_open( aExodusFileName.c_str(), EX_READ, &tCPUWordSize, &tIOWordSize, &aVersion );
        }
        else
        {
            mExoID = ex_open( aExodusFileName.c_str(), EX_WRITE, &tCPUWordSize, &tIOWordSize, &aVersion );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::close_file( bool aRename )
    {
        // check that mesh is open
        MORIS_ERROR( mExoID > 0, "Exodus cannot be saved as it is not open\n." );

        ex_close( mExoID );
        mExoID = -1;

        if ( aRename )
        {
            MORIS_ERROR( std::rename( mTempFileName.c_str(), mPermFileName.c_str() ) == 0,
                    "Cannot save exodus file: %s as %s",
                    mTempFileName.c_str(),
                    mPermFileName.c_str() );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::write_mesh(
            std::string        aFilePath,
            const std::string& aFileName,
            std::string        aTempPath,
            const std::string& aTempName )
    {
        MORIS_ERROR( mMesh != nullptr, "No mesh has been given to the Exodus Writer!" );

        this->create_init_mesh_file(
                std::move( aFilePath ),
                aFileName,
                std::move( aTempPath ),
                aTempName );

        this->write_nodes();
        this->write_node_sets();
        this->write_blocks();
        this->write_side_sets();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::write_points(
            std::string        aFilePath,
            const std::string& aFileName,
            std::string        aTempPath,
            const std::string& aTempName,
            Matrix< DDRMat >   aCoordinates )
    {
        // Create the actual file
        this->create_file(
                std::move( aFilePath ),
                aFileName,
                std::move( aTempPath ),
                aTempName );

        // Initialize database
        int tNumDimensions = aCoordinates.n_cols();
        int tNumPoints     = aCoordinates.n_rows();

        ex_put_init( mExoID, "MTK", tNumDimensions, tNumPoints, 1, 1, 0, 0 );

        // Set the point coordinates
        int  tSpatialDim = aCoordinates.n_cols();
        bool tYDim       = tSpatialDim >= 2;
        bool tZDim       = tSpatialDim >= 3;

        // Set up coordinate and node map arrays based on the number of vertices
        MORIS_ERROR( aCoordinates.n_rows() > 0, "Points need to be given to create a point field" );

        Matrix< IdMat > tNodeMap( aCoordinates.n_rows(), 1, 0 );

        // Coordinate arrays
        Matrix< DDRMat > tXCoordinates( aCoordinates.n_rows(), 1, 0.0 );
        Matrix< DDRMat > tYCoordinates( aCoordinates.n_rows(), 1, 0.0 );
        Matrix< DDRMat > tZCoordinates( aCoordinates.n_rows(), 1, 0.0 );

        for ( uint tNodeIndex = 0; tNodeIndex < aCoordinates.n_rows(); tNodeIndex++ )
        {
            // Place in coordinate arrays
            tXCoordinates( tNodeIndex ) = aCoordinates( tNodeIndex, 0 );
            tYCoordinates( tNodeIndex ) = aCoordinates( tNodeIndex, 1 * tYDim ) * tYDim;
            tZCoordinates( tNodeIndex ) = aCoordinates( tNodeIndex, 2 * tZDim ) * tZDim;

            // Get global ids for id map
            tNodeMap( tNodeIndex ) = tNodeIndex + 1;
        }

        // Write coordinates
        ex_put_coord( mExoID, tXCoordinates.data(), tYCoordinates.data(), tZCoordinates.data() );

        // Write node id map
        ex_put_id_map( mExoID, EX_NODE_MAP, tNodeMap.data() );

        // Create single block
        if ( tNumDimensions <= 2 )
        {
            ex_put_block( mExoID, EX_ELEM_BLOCK, 1, "CIRCLE", tNumPoints, 1, 0, 0, 0 );
        }
        else
        {
            ex_put_block( mExoID, EX_ELEM_BLOCK, 1, "SPHERE", tNumPoints, 1, 0, 0, 0 );
        }
        std::string tBlockName( "Points" );
        ex_put_name( mExoID, EX_ELEM_BLOCK, 1, tBlockName.c_str() );

        // Create point elements
        Matrix< IndexMat > tConnectivityArray( tNumPoints, 1 );
        Matrix< IndexMat > tElementIdMap( tNumPoints, 1 );

        for ( int tNodeIndex = 0; tNodeIndex < tNumPoints; tNodeIndex++ )
        {
            tConnectivityArray( tNodeIndex ) = tNodeIndex + 1;
            tElementIdMap( tNodeIndex )      = tNodeIndex + 1;
        }
        ex_put_conn( mExoID, EX_ELEM_BLOCK, 1, tConnectivityArray.data(), nullptr, nullptr );

        // Write the element map
        ex_put_id_map( mExoID, EX_ELEM_MAP, tElementIdMap.data() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::set_point_fields( Vector< std::string > aFieldNames )
    {
        // Set the field names
        if ( aFieldNames.size() > 0 )
        {
            // Write the number of nodal fields
            ex_put_variable_param( mExoID, EX_NODAL, aFieldNames.size() );

            // Write the nodal field names and store as a map
            for ( uint tFieldIndex = 0; tFieldIndex < aFieldNames.size(); tFieldIndex++ )
            {
                ex_put_variable_name( mExoID, EX_NODAL, tFieldIndex + 1, aFieldNames( tFieldIndex ).c_str() );
                mNodalFieldNamesMap[ aFieldNames( tFieldIndex ) ] = tFieldIndex;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::set_nodal_fields( Vector< std::string > aFieldNames )
    {
        if ( aFieldNames.size() > 0 && mNumNodes > 0 )
        {
            // Write the number of nodal fields
            ex_put_variable_param( mExoID, EX_NODAL, aFieldNames.size() );

            // Write the nodal field names and store as a map
            for ( uint tFieldIndex = 0; tFieldIndex < aFieldNames.size(); tFieldIndex++ )
            {
                ex_put_variable_name( mExoID, EX_NODAL, tFieldIndex + 1, aFieldNames( tFieldIndex ).c_str() );
                mNodalFieldNamesMap[ aFieldNames( tFieldIndex ) ] = tFieldIndex;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::set_elemental_fields( Vector< std::string > aFieldNames )
    {
        if ( aFieldNames.size() > 0 && mNumUniqueExodusElements > 0 )
        {
            // Write the number of elemental fields
            ex_put_variable_param( mExoID, EX_ELEM_BLOCK, aFieldNames.size() );

            // Write the elemental field names and store as a map
            for ( uint tFieldIndex = 0; tFieldIndex < aFieldNames.size(); tFieldIndex++ )
            {
                ex_put_variable_name( mExoID, EX_ELEM_BLOCK, tFieldIndex + 1, aFieldNames( tFieldIndex ).c_str() );
                mElementalFieldNamesMap[ aFieldNames( tFieldIndex ) ] = tFieldIndex;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::set_side_set_fields( Vector< std::string > aFieldNames )
    {
        if ( aFieldNames.size() > 0 && mNumUniqueExodusElements > 0 )
        {
            // Write the number of side set fields
            ex_put_variable_param( mExoID, EX_SIDE_SET, aFieldNames.size() );

            // Write the side set field names and store as a map
            for ( uint tFieldIndex = 0; tFieldIndex < aFieldNames.size(); tFieldIndex++ )
            {
                ex_put_variable_name( mExoID, EX_SIDE_SET, tFieldIndex + 1, aFieldNames( tFieldIndex ).c_str() );
                mSideSetFieldNamesMap[ aFieldNames( tFieldIndex ) ] = tFieldIndex;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::set_global_variables( Vector< std::string > aVariableNames )
    {
        if ( aVariableNames.size() > 0 )
        {
            // Write the number of global fields
            ex_put_variable_param( mExoID, EX_GLOBAL, aVariableNames.size() );

            // Write the global field names and store as a map
            for ( uint tVariableIndex = 0; tVariableIndex < aVariableNames.size(); tVariableIndex++ )
            {
                ex_put_variable_name( mExoID, EX_GLOBAL, tVariableIndex + 1, aVariableNames( tVariableIndex ).c_str() );
                mGlobalVariableNamesMap[ aVariableNames( tVariableIndex ) ] = tVariableIndex;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::save_mesh()
    {
        namespace fs = std::filesystem;

        // check that mesh is open
        MORIS_ERROR( mExoID > 0,
                "Writer_Exodus::save_mesh() - Exodus cannot be saved as it is not open\n." );

        // close mesh
        ex_close( mExoID );

        // write log information
        MORIS_LOG( "Copying %s to %s.", mTempFileName.c_str(), mPermFileName.c_str() );

        // copy temporary file on permanent file
        fs::path source      = mTempFileName;
        fs::path destination = mPermFileName;

        try
        {
            fs::copy_file( source, destination, fs::copy_options::overwrite_existing );
        } catch ( fs::filesystem_error& e )
        {
            MORIS_ERROR( false, "Writer_Exodus::save_mesh - copying %s to %s failed.", mTempFileName.c_str(), mPermFileName.c_str() );
        }

        // open mesh file again
        int   tCPUWordSize = sizeof( real ), tIOWordSize = 0;
        float tVersion;

        mExoID = ex_open(
                mTempFileName.c_str(),
                EX_WRITE,
                &tCPUWordSize,
                &tIOWordSize,
                &tVersion );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::set_time( real aTimeValue )
    {
        ex_put_time( mExoID, ++mTimeStep, &aTimeValue );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::write_point_field(
            const std::string&      aFieldName,
            const Matrix< DDRMat >& aFieldValues )
    {
        // skip if no nodal values exist
        if ( aFieldValues.numel() == 0 )
        {
            return;
        }

        // Field name to index
        int tFieldIndex = mNodalFieldNamesMap[ aFieldName ];

        MORIS_ASSERT( mNodalFieldNamesMap.size() == mNodalFieldNamesMap.size(),
                "%s is not a point field name on this mesh!",
                aFieldName.c_str() );

        // Write the field
        int tErrMsg = ex_put_var(
                mExoID,
                mTimeStep,
                EX_NODAL,
                tFieldIndex + 1,
                1,
                aFieldValues.numel(),
                aFieldValues.data() );

        // Check for error
        MORIS_ERROR( tErrMsg == 0,
                "Point field %s could not be written to exodus file.",
                aFieldName.c_str() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::write_nodal_field(
            const std::string&      aFieldName,
            const Matrix< DDRMat >& aFieldValues )
    {
        // skip if no nodal values exist
        if ( aFieldValues.numel() == 0 )
        {
            MORIS_ASSERT( mNumNodes == 0,
                    "Number of field values is zero but not number of nodes.\n" );

            return;
        }

        // Field name to index
        int tFieldIndex = mNodalFieldNamesMap[ aFieldName ];

        MORIS_ASSERT( mNodalFieldNamesMap.size() == mNodalFieldNamesMap.size(),
                "%s is not a nodal field name on this mesh!",
                aFieldName.c_str() );

        // Check number of field values = number of nodes
        MORIS_ERROR( aFieldValues.numel() == mMesh->get_num_nodes(),
                "%s field was attempted to be written with %li values, but there are %i nodes in this mesh.",
                aFieldName.c_str(),
                aFieldValues.numel(),
                mMesh->get_num_nodes() );

        // Ensure that time step is larger than or equal 1
        int tTimeStep = mTimeStep > 1 ? mTimeStep : 1;

        // Write the field
        int tErrMsg = ex_put_var(
                mExoID,
                tTimeStep,
                EX_NODAL,
                tFieldIndex + 1,
                1,
                aFieldValues.numel(),
                aFieldValues.data() );

        // Check for error
        MORIS_ERROR( tErrMsg == 0,
                "Nodal field %s could not be written to exodus file.",
                aFieldName.c_str() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::write_elemental_field(
            const std::string&      aBlockName,
            const std::string&      aFieldName,
            const Matrix< DDRMat >& aFieldValues )
    {
        // skip if no elemental values exist
        if ( aFieldValues.numel() == 0 )
        {
            return;
        }

        MORIS_ERROR( mBlockNamesMap.key_exists( aBlockName ),
                "%s is not a block name on this mesh!",
                aBlockName.c_str() );

        // Block name to local index of non-empty blocks
        int tBlockIndex = mBlockNamesMap[ aBlockName ];

        // Check that block index is valid
        int tNumBlocks = ex_inquire_int( mExoID, EX_INQ_ELEM_BLK );
        MORIS_ERROR(
                tNumBlocks > tBlockIndex,
                "Writer_Exodus::write_elemental_field() - Index of block set is larger than number of blocks" );

        // Field name to index
        int tFieldIndex = mElementalFieldNamesMap[ aFieldName ];
        MORIS_ERROR(
                (uint)tFieldIndex < mElementalFieldNamesMap.size(),
                "Writer_Exodus::write_elemental_field() - '%s' is not an elemental field name on this mesh!",
                aFieldName.c_str() );

        // Check number of field values = number of elements
        MORIS_ERROR( aFieldValues.numel() == mMesh->get_set_cells( aBlockName ).size(),
                "%s field was attempted to be written with %li values, but there are %li  elements in block %s",
                aFieldName.c_str(),
                aFieldValues.numel(),
                mMesh->get_set_cells( aBlockName ).size(),
                aBlockName.c_str() );

        // Check that number of field values = number of element stored in mesh
        ex_block tBlockInfo;
        tBlockInfo.id   = tBlockIndex + 1;
        tBlockInfo.type = EX_ELEM_BLOCK;

        ex_get_block_param( mExoID, &tBlockInfo );

        MORIS_ERROR( tBlockInfo.num_entry == (int)aFieldValues.numel(),
                "Number of entries in field does not match number of elements stored in mesh for current block." );

        // Ensure that time step is larger than or equal 1
        int tTimeStep = mTimeStep > 1 ? mTimeStep : 1;

        // Write the field
        int tErrMsg = ex_put_var(
                mExoID,
                tTimeStep,
                EX_ELEM_BLOCK,
                tFieldIndex + 1,
                tBlockIndex + 1,
                aFieldValues.numel(),
                aFieldValues.data() );

        // Check for error
        MORIS_ERROR( tErrMsg == 0,
                "Elemental field %s could not be written for element block %s to exodus file.",
                aFieldName.c_str(),
                aBlockName.c_str() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::write_side_set_field(
            const std::string&      aSideSetName,
            const std::string&      aFieldName,
            const Matrix< DDRMat >& aFieldValues )
    {
        // skip if no elemental values exist
        uint tNumFieldEntries = aFieldValues.numel();
        if ( tNumFieldEntries == 0 )
        {
            return;
        }

        // Check that side set is valid
        MORIS_ERROR(
                mSideSetNamesMap.key_exists( aSideSetName ),
                "Writer_Exodus::write_side_set_field() - "
                "The side set (%s) the data (size: %i) is supposed to be written to does not exist in exodus mesh.",
                aSideSetName.c_str(),
                tNumFieldEntries );

        // Block name to local index of non-empty blocks
        int tSideSetIndex = mSideSetNamesMap[ aSideSetName ];

        // Field name to index
        int tFieldIndex = mSideSetFieldNamesMap[ aFieldName ];
        MORIS_ERROR(
                (uint)tFieldIndex < mSideSetFieldNamesMap.size(),
                "Writer_Exodus::write_side_set_field() - '%s' is not an elemental field name on this mesh.",
                aFieldName.c_str() );

        // Ensure that time step is larger than or equal 1
        int tTimeStep = mTimeStep > 1 ? mTimeStep : 1;

        // remove values associated with side sets that are not part of the exodus mesh on the current proc's domain
        // Note: this can happen when a sideset snaps to a background facet on the domain boundary where the leader side is in the aura.

        // count number of facets from the current side set that are part of the exodus mesh
        uint tNumUsedFacets = 0;
        for ( uint iFacet = 0; iFacet < mFacetUsedInExodus( tSideSetIndex ).size(); iFacet++ )
        {
            tNumUsedFacets += mFacetUsedInExodus( tSideSetIndex )( iFacet );
        }

        // if no facets are part of the exodus mesh, skip writing the field
        if ( tNumUsedFacets == 0 )
        {
            return;
        }

        // clean out field values which are not used
        Matrix< DDRMat > tUsedFieldValues( tNumUsedFacets, 1 );
        uint             tCounter = 0;

        for ( uint iFieldValue = 0; iFieldValue < tNumFieldEntries; iFieldValue++ )
        {
            if ( mFacetUsedInExodus( tSideSetIndex )( iFieldValue ) )
            {
                tUsedFieldValues( tCounter ) = aFieldValues( iFieldValue );
                tCounter++;
            }
        }

        // Write the field
        int tErrMsg = ex_put_var(
                mExoID,
                tTimeStep,
                EX_SIDE_SET,
                tFieldIndex + 1,
                tSideSetIndex + 1,
                tUsedFieldValues.numel(),
                tUsedFieldValues.data() );

        // Check for error
        MORIS_ERROR( tErrMsg == 0,
                "Writer_Exodus::write_side_set_field() - "
                "Side set field '%s' could not be written for side set '%s' to exodus file.",
                aFieldName.c_str(),
                aSideSetName.c_str() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::write_global_variables(
            Vector< std::string >&  aVariableNames,
            const Matrix< DDRMat >& aVariableValues )
    {
        // number of global variables
        uint tNumVariables = aVariableNames.size();

        // check that number of variables names is not larger than number of values
        MORIS_ASSERT( tNumVariables <= aVariableValues.numel(),
                "Number of global variables names larger than number of values." );

        // allocated vector of sorted global variables
        Matrix< DDRMat > tSortedValues( tNumVariables, 1 );

        for ( uint ig = 0; ig < tNumVariables; ++ig )
        {
            // Variable name to index
            uint tVariableIndex = mGlobalVariableNamesMap[ aVariableNames( ig ) ];

            // check that index is valid
            MORIS_ASSERT( tVariableIndex < tNumVariables, "Global variable index too large." );

            // store global variable in sorted list
            tSortedValues( tVariableIndex ) = aVariableValues( ig );
        }

        // Write the variables
        int tErrMsg = ex_put_var(
                mExoID,
                mTimeStep,
                EX_GLOBAL,
                1,
                0,
                tNumVariables,
                tSortedValues.data() );

        // check for error
        MORIS_ERROR( tErrMsg == 0, "Global variables could not be written to exodus file." );
    }

    //--------------------------------------------------------------------------------------------------------------
    // Private
    //--------------------------------------------------------------------------------------------------------------

    void
    Writer_Exodus::create_file(
            std::string        aFilePath,
            const std::string& aFileName,
            std::string        aTempPath,
            const std::string& aTempName )
    {
        MORIS_ERROR( mExoID == -1,
                "Exodus file is currently open, call close_file() before creating a new one." );

        // check that filename does not contain any path information, except for ./
        MORIS_ERROR( isFileNameOnly( aFileName ),
                "Exodus file name (%s) must not contain path information.",
                aFileName.c_str() );

        MORIS_ERROR( isFileNameOnly( aTempName ),
                "Exodus temporary file name (%s) must not contain path information.",
                aTempName.c_str() );

        // Add temporary and permanent file names to file paths
        if ( !aFilePath.empty() )
        {
            aFilePath += "/";
        }

        if ( !aTempPath.empty() )
        {
            aTempPath += "/";
        }

        // check and if necessary create the temporary path
        create_directory( aTempPath );
        create_directory( aFilePath );

        mTempFileName = aTempPath + aTempName;
        mPermFileName = aFilePath + aFileName;

        // Make file name parallel, if necessary
        if ( par_size() > 1 )
        {
            // Get par size and rank as strings
            std::string tParSizeStr     = std::to_string( par_size() );
            std::string tParRankBaseStr = std::to_string( par_rank() );

            // Make sure all par rank strings have the same number of leading zeros
            std::string tParRankStr =
                    std::string( tParSizeStr.length() - tParRankBaseStr.length(), '0' ).append( tParRankBaseStr );

            // Append to temporary and permanent file names
            mTempFileName += "." + tParSizeStr + "." + tParRankStr;
            mPermFileName += "." + tParSizeStr + "." + tParRankStr;
        }

        // Create the database
        int cpu_ws = sizeof( real );    // word size in bytes of the floating point variables used in moris
        int io_ws  = sizeof( real );    // word size as stored in exodus

        mExoID = ex_create( mTempFileName.c_str(), EX_CLOBBER, &cpu_ws, &io_ws );

        MORIS_ERROR( mExoID > -1, "Exodus file cannot be created: %s", mTempFileName.c_str() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Writer_Exodus::create_init_mesh_file(
            std::string        aFilePath,
            const std::string& aFileName,
            std::string        aTempPath,
            const std::string& aTempName )
    {
        // Create the actual file
        this->create_file(
                std::move( aFilePath ),
                aFileName,
                std::move( aTempPath ),
                aTempName );

        // Number of dimensions
        int tNumDimensions = mMesh->get_spatial_dim();

        // Number of nodes
        mNumNodes = mMesh->get_num_nodes();

        // Number of elements in mesh (may be different from number of elements in blocks)
        mNumMtkElements = mMesh->get_num_elems();

        // Determine number of non-empty sets
        this->get_node_sets();
        this->get_side_sets();
        this->get_block_sets();

        // Initialize database
        ex_put_init(
                mExoID,
                "MTK",
                tNumDimensions,
                mNumNodes,
                mNumUniqueExodusElements,
                (int)mElementBlockIndices.size(),
                (int)mNodeSetIndices.size(),
                (int)mSideSetIndices.size() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Writer_Exodus::get_node_sets()
    {
        Vector< std::string > tNodeSetNames = mMesh->get_set_names( EntityRank::NODE );

        // Determine number of nodes in each local node set
        Matrix< DDUMat > tNumLocalNodesOnNodeSets( tNodeSetNames.size(), 1 );

        for ( uint tNodeSetIndex = 0; tNodeSetIndex < tNodeSetNames.size(); tNodeSetIndex++ )
        {
            // FIXME: should be replaced by call to MTK to just get number of nodes in node set (local or global?)
            Matrix< IndexMat > tNodeIndices =
                    mMesh->get_set_entity_loc_inds(
                            EntityRank::NODE,
                            tNodeSetNames( tNodeSetIndex ) );

            tNumLocalNodesOnNodeSets( tNodeSetIndex ) = tNodeIndices.numel();
        }

        // Sum number of nodes in node sets across all procs
        Matrix< DDUMat > tNumGlobalNodesOnNodesSets = sum_all_matrix( tNumLocalNodesOnNodeSets );

        // Determine non-empty node sets
        for ( uint tNodeSetIndex = 0; tNodeSetIndex < tNodeSetNames.size(); tNodeSetIndex++ )
        {
            if ( tNumGlobalNodesOnNodesSets( tNodeSetIndex ) > 0 )
            {
                mNodeSetIndices.push_back( tNodeSetIndex );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Writer_Exodus::get_side_sets()
    {
        // Determine number of non-empty side sets across all procs
        Vector< std::string > tSideSetNames = mMesh->get_set_names( mMesh->get_facet_rank() );

        // Determine number of sides in each local side set
        Matrix< DDUMat > tNumLocalSidesOnSideSets( tSideSetNames.size(), 1 );

        for ( uint tSideSetIndex = 0; tSideSetIndex < tSideSetNames.size(); tSideSetIndex++ )
        {
            // FIXME: should be replaced by call to MTK to just get number of faces in side set (local or global?)
            Matrix< IndexMat > tSideSetElements;
            Matrix< IndexMat > tSideSetOrdinals;

            mMesh->get_sideset_elems_loc_inds_and_ords(
                    tSideSetNames( tSideSetIndex ),
                    tSideSetElements,
                    tSideSetOrdinals );

            tNumLocalSidesOnSideSets( tSideSetIndex ) = tSideSetOrdinals.numel();
        }

        // Sum number of sides in sides sets across all procs
        Matrix< DDUMat > tNumGlobalSidesOnSideSets = sum_all_matrix( tNumLocalSidesOnSideSets );

        // Determine non-empty side sets
        for ( uint tSideSetIndex = 0; tSideSetIndex < tSideSetNames.size(); tSideSetIndex++ )
        {
            if ( tNumGlobalSidesOnSideSets( tSideSetIndex ) > 0 )
            {
                mSideSetIndices.push_back( tSideSetIndex );
            }
        }

        // resize to EXODUS_MAX_NUM_SIDESET
        if ( mSideSetIndices.size() > EXODUS_MAX_NUM_SIDESET )
        {
            mSideSetIndices.resize( EXODUS_MAX_NUM_SIDESET );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Writer_Exodus::get_block_sets()
    {
        Vector< std::string > tBlockNames = mMesh->get_set_names( EntityRank::ELEMENT );

        // Initialize map from mesh indices to exodus indices
        mMtkExodusElementIndexMap.set_size( mNumMtkElements, 1, MORIS_INDEX_MAX );

        // Determine number of sides in each local side set
        Matrix< DDUMat > tNumLocalElemsInBlockSets( tBlockNames.size(), 1 );

        // Initialize exodus index
        mNumUniqueExodusElements = 0;

        // Initialize total number of elements in blocks
        mNumTotalExodusElements = 0;

        for ( uint tBlockIndex = 0; tBlockIndex < tBlockNames.size(); tBlockIndex++ )
        {
            Matrix< IndexMat > tElementIndices = mMesh->get_element_indices_in_block_set( tBlockIndex );

            uint tNumElementsInBlock = tElementIndices.length();

            // create exodus index
            for ( uint tInd = 0; tInd < tNumElementsInBlock; tInd++ )
            {
                // mesh-based element index
                moris_index tElementIndex = tElementIndices( tInd );

                // Check that index is not larger than local number of elements
                MORIS_ASSERT( tElementIndex < (moris_index)mNumMtkElements,
                        "Element index larger or equal than number of elements in MTK mesh: index = %d (%d).\n",
                        tElementIndex,
                        mNumMtkElements );

                // create exodus index
                if ( mMtkExodusElementIndexMap( tElementIndex ) == MORIS_INDEX_MAX )
                {
                    mMtkExodusElementIndexMap( tElementIndex ) = mNumUniqueExodusElements++;
                }
            }
            // add number of elements in this block to total number of elements in exodus mesh
            mNumTotalExodusElements += tNumElementsInBlock;

            // save number of elements in block
            tNumLocalElemsInBlockSets( tBlockIndex ) = tNumElementsInBlock;
        }

        // Sum number of elements in block sets across all procs
        Matrix< DDUMat > tNumGlobalElemsInBlockSets = sum_all_matrix( tNumLocalElemsInBlockSets );

        // Determine non-empty side sets
        for ( uint tBlockIndex = 0; tBlockIndex < tBlockNames.size(); tBlockIndex++ )
        {
            if ( tNumGlobalElemsInBlockSets( tBlockIndex ) > 0 )
            {
                mElementBlockIndices.push_back( tBlockIndex );
            }
        }

        // Check that number of elements in blocks is not larger than number of elements in mesh
        MORIS_ASSERT( mNumUniqueExodusElements <= mNumMtkElements,
                "Writer_Exodus::get_block_sets - Number of unique elements in blocks larger than number of elements in mesh: %d vs %d\n",
                mNumUniqueExodusElements,
                mNumMtkElements );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Writer_Exodus::write_nodes()
    {
        // spatial dimension
        int  tSpatialDim = mMesh->get_spatial_dim();
        bool tYDim       = tSpatialDim >= 2;
        bool tZDim       = tSpatialDim >= 3;

        // Set up coordinate and node map arrays based on the number of vertices
        Matrix< IdMat > tNodeMap( mNumNodes, 1, 0 );

        // when using ad-hoc index map, get processor offset and add 1 as exodus used 1 based indices
        uint tProcOffset = 1;
        if ( !mMtkIndexMap )
        {
            tProcOffset = get_processor_offset( mNumNodes ) + 1;
        }

        // Coordinate arrays
        Matrix< DDRMat > tXCoordinates( mNumNodes, 1, 0.0 );
        Matrix< DDRMat > tYCoordinates( mNumNodes, 1, 0.0 );
        Matrix< DDRMat > tZCoordinates( mNumNodes, 1, 0.0 );

        for ( uint tNodeIndex = 0; tNodeIndex < mNumNodes; tNodeIndex++ )
        {
            // Get coordinates
            Matrix< DDRMat > tNodeCoordinates = mMesh->get_node_coordinate( tNodeIndex );

            // Place in coordinate arrays
            tXCoordinates( tNodeIndex ) = tNodeCoordinates( 0 );
            tYCoordinates( tNodeIndex ) = tNodeCoordinates( 1 * tYDim ) * tYDim;
            tZCoordinates( tNodeIndex ) = tNodeCoordinates( 2 * tZDim ) * tZDim;

            // Get global ids for id map using either MTK or ad-hoc node ID map
            if ( mMtkIndexMap )
            {
                tNodeMap( tNodeIndex ) = mMesh->get_glb_entity_id_from_entity_loc_index( tNodeIndex, EntityRank::NODE );
            }
            else
            {
                tNodeMap( tNodeIndex ) = tProcOffset + tNodeIndex;
            }
        }

        // Write coordinate names
        const char* tmp[ 3 ] = { "x", "y", "z" };
        char*       coord_names[ 3 ];
        memcpy( coord_names, tmp, sizeof( coord_names ) );

        ex_put_coord_names(
                mExoID,
                coord_names );

        // Write coordinates
        ex_put_coord(
                mExoID,
                tXCoordinates.data(),
                tYCoordinates.data(),
                tZCoordinates.data() );

        // Write node id map
        ex_put_id_map(
                mExoID,
                EX_NODE_MAP,
                tNodeMap.data() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Writer_Exodus::write_node_sets()
    {
        // Get the number of node sets and their names
        Vector< std::string > tNodeSetNames = mMesh->get_set_names( EntityRank::NODE );

        // Loop through all non-empty node sets
        for ( uint tIndex = 0; tIndex < mNodeSetIndices.size(); tIndex++ )
        {
            // get global index
            uint tNodeSetIndex = mNodeSetIndices( tIndex );

            // get node indices of side set
            Matrix< IndexMat > tNodeIndices =
                    mMesh->get_set_entity_loc_inds(
                            EntityRank::NODE,
                            tNodeSetNames( tNodeSetIndex ) );

            ex_put_set_param(
                    mExoID,
                    EX_NODE_SET,
                    tIndex + 1,
                    tNodeIndices.numel(),
                    0 );

            ex_put_set(
                    mExoID,
                    EX_NODE_SET,
                    tIndex + 1,
                    tNodeIndices.data(),
                    nullptr );

            ex_put_name(
                    mExoID,
                    EX_NODE_SET,
                    tIndex + 1,
                    tNodeSetNames( tNodeSetIndex ).c_str() );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Writer_Exodus::write_blocks()
    {
        // create ad-hoc element IDs FIXME: should be handled by MTK
        uint tProcOffset = get_processor_offset( mNumUniqueExodusElements ) + 1;

        // ids of exodus elements
        Matrix< IdMat > tExodusElementIds;

        // generate ad-hoc element ID map when not using MTK map
        if ( !mMtkIndexMap )
        {
            // set size of exodus element map and fill it
            tExodusElementIds.set_size( mNumUniqueExodusElements, 1 );
            tExodusElementIds.fill( MORIS_ID_MAX );

            for ( uint tInd = 0; tInd < mNumUniqueExodusElements; ++tInd )
            {
                tExodusElementIds( tInd ) = tProcOffset + tInd;
            }
        }

        // position of exodus element in file (only last position stored)
        mExodusElementIndexOrderMap.set_size( mNumUniqueExodusElements, 1, MORIS_ID_MAX );

        // ids for each element in blocks; note: an element might be in multiple blocks
        Matrix< IdMat > tExodusTotalElementIds( mNumTotalExodusElements, 1 );

        // All of the block names
        Vector< std::string > tBlockNames = mMesh->get_set_names( EntityRank::ELEMENT );

        // Initialize element counter: order in which elements are written to exodus file
        uint tElemCounter = 0;

        // Loop through all non-empty blocks
        for ( uint iBlockIndexInExoMesh = 0; iBlockIndexInExoMesh < mElementBlockIndices.size(); iBlockIndexInExoMesh++ )
        {
            // Get global index
            uint tBlockIndexInInputMesh = mElementBlockIndices( iBlockIndexInExoMesh );

            // Get local indices and IDs of elements in current block
            Matrix< IndexMat > tElementIndices = mMesh->get_element_indices_in_block_set( tBlockIndexInInputMesh );
            Matrix< IdMat >    tElementIDs     = mMesh->get_element_ids_in_block_set( tBlockIndexInInputMesh );

            MORIS_ASSERT( tElementIndices.length() == tElementIDs.length(),
                    "Writer_Exodus::write_blocks - Vector of element IDs must match vector of element indices for Exodus file writing." );

            // Number of elements in the block
            uint tNumElementsInBlock = tElementIndices.length();

            // Add name to map
            mBlockNamesMap[ tBlockNames( tBlockIndexInInputMesh ) ] = iBlockIndexInExoMesh;

            // Get the CellTopology of this block
            enum CellTopology tMorisBlockTopology = mMesh->get_blockset_topology( tBlockNames( tBlockIndexInInputMesh ) );

            const char* tExodusBlockTopology = this->get_exodus_block_topology( tMorisBlockTopology );

            // Get the number of nodes/edges/faces/attributes per element
            int tNumNodesPerElement      = this->get_nodes_per_element( tMorisBlockTopology );
            int tNumEdgesPerElement      = 0;
            int tNumFacesPerElement      = 0;
            int tNumAttributesPerElement = 0;

            // Make a block and name it
            ex_put_block(
                    mExoID,
                    EX_ELEM_BLOCK,
                    iBlockIndexInExoMesh + 1,
                    tExodusBlockTopology,
                    tNumElementsInBlock,
                    tNumNodesPerElement,
                    tNumEdgesPerElement,
                    tNumFacesPerElement,
                    tNumAttributesPerElement );

            ex_put_name(
                    mExoID,
                    EX_ELEM_BLOCK,
                    iBlockIndexInExoMesh + 1,
                    tBlockNames( tBlockIndexInInputMesh ).c_str() );

            // Construct matrix of node indices per element
            Matrix< IndexMat > tConnectivityArray( tNumNodesPerElement * tNumElementsInBlock, 1, 0 );

            // Loop through the elements in this block
            uint tConnectivityIndex = 0;

            for ( uint tElementNumber = 0; tElementNumber < tNumElementsInBlock; tElementNumber++ )
            {
                // Local element index
                moris_index tElementIndex = tElementIndices( tElementNumber );

                // Get the vertex indices of this element
                Matrix< IndexMat > tNodeIndices = mMesh->get_nodes_connected_to_element_loc_inds( tElementIndex );

                MORIS_ASSERT( tNodeIndices.numel() >= (uint)tNumNodesPerElement,
                        "Writer_Exodus::write_blocks - number of nodes per element too small for element type." );

                // Build connectivity vector, add 1 since exodus uses 1-based indices
                for ( int tNodeNum = 0; tNodeNum < tNumNodesPerElement; tNodeNum++ )
                {
                    tConnectivityArray( tConnectivityIndex ) = tNodeIndices( tNodeNum ) + 1;
                    tConnectivityIndex++;
                }

                // Check that mesh-based index an exodus index has been assigned
                MORIS_ASSERT( mMtkExodusElementIndexMap( tElementIndex ) != MORIS_INDEX_MAX,
                        "Writer_Exodus::write_blocks - Mesh based element index has not been assigned an exodus index.\n" );

                // Get exodus index
                moris_index tExodusIndex = mMtkExodusElementIndexMap( tElementIndex );

                // Assign tElemCounter to exodus index
                mExodusElementIndexOrderMap( tExodusIndex ) = tElemCounter;

                // Assign exodus element id using either MTK or ad-hoc map
                if ( mMtkIndexMap )
                {
                    tExodusTotalElementIds( tElemCounter++ ) = tElementIDs( tElementNumber );
                }
                else
                {
                    tExodusTotalElementIds( tElemCounter++ ) = tExodusElementIds( tExodusIndex );
                }
            }

            // Write connectivity
            ex_put_conn(
                    mExoID,
                    EX_ELEM_BLOCK,
                    iBlockIndexInExoMesh + 1,
                    tConnectivityArray.data(),
                    nullptr,
                    nullptr );
        }

        // Write the element map
        ex_put_id_map(
                mExoID,
                EX_ELEM_MAP,
                tExodusTotalElementIds.data() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Writer_Exodus::write_side_sets()
    {
        // Get side set names
        Vector< std::string > tSideSetNames = mMesh->get_set_names( mMesh->get_facet_rank() );

        // get the number of side sets in the exo mesh
        uint tNumSideSets = mSideSetIndices.size();

        // set size for the punch card marking which facets are used in exodus mesh
        mFacetUsedInExodus.resize( tNumSideSets );

        // Loop through all non-empty side sets
        for ( uint iSideSetInExoMesh = 0; iSideSetInExoMesh < tNumSideSets; iSideSetInExoMesh++ )
        {
            // Get the index of the current set in the input mesh
            uint tSideSetIndexInInputMesh = mSideSetIndices( iSideSetInExoMesh );

            // get the name/label of the set
            const std::string& tSetLabel = tSideSetNames( tSideSetIndexInInputMesh );

            // Add name to map
            mSideSetNamesMap[ tSetLabel ] = iSideSetInExoMesh;

            // Get the side set element ids
            Matrix< IndexMat > tIgElemIndices;
            Matrix< IndexMat > tIgElemSideOrdinals;
            mMesh->get_sideset_elems_loc_inds_and_ords(
                    tSetLabel,
                    tIgElemIndices,
                    tIgElemSideOrdinals );
            uint tNumFacetsInSideSet = tIgElemSideOrdinals.numel();

            // set size for the punch card marking which facets are used in exodus mesh
            mFacetUsedInExodus( iSideSetInExoMesh ).resize( tNumFacetsInSideSet, true );

            // Side counter
            uint tSideCounter = 0;

            // Loop over all sides in set
            for ( uint iIgElemInSideSet = 0; iIgElemInSideSet < tNumFacetsInSideSet; iIgElemInSideSet++ )
            {
                // Get mesh-based index
                moris_index tIgElemIndex = tIgElemIndices( iIgElemInSideSet );

                // Check that mesh-based index an exodus index has been assigned;
                // if not it is not part of a block on this processor and will be skipped
                if ( mMtkExodusElementIndexMap( tIgElemIndex ) != MORIS_INDEX_MAX )
                {
                    // Get index of IG elem in exodus mesh
                    moris_index tIgElemIndexInExo = mMtkExodusElementIndexMap( tIgElemIndex );

                    // Check that position has been assigned
                    MORIS_ASSERT(
                            mExodusElementIndexOrderMap( tIgElemIndexInExo ) != MORIS_INDEX_MAX,
                            "Writer_Exodus::write_side_sets() - Exodus element index has not been assigned a position in exodus file." );

                    // Get position in exodus file this element has been written to
                    moris_index tExoPos = mExodusElementIndexOrderMap( tIgElemIndexInExo );

                    tIgElemIndices( tSideCounter ) = tExoPos + 1;    // add 1 to side ordinal to match exodus definition
                    tIgElemSideOrdinals( tSideCounter )++;           // add 1 to side ordinal to match exodus definition
                    tSideCounter++;
                }
                // if facet is not part of any block, mark it so facets don't get constructed on non-existent elements in exo mesh
                else
                {
                    mFacetUsedInExodus( iSideSetInExoMesh )( iIgElemInSideSet ) = false;
                }
            }

            // Write the side set with actual number of sides in set
            ex_put_set_param(
                    mExoID,
                    EX_SIDE_SET,
                    iSideSetInExoMesh + 1,
                    tSideCounter,
                    0 );

            ex_put_set(
                    mExoID,
                    EX_SIDE_SET,
                    iSideSetInExoMesh + 1,
                    tIgElemIndices.data(),
                    tIgElemSideOrdinals.data() );

            ex_put_name(
                    mExoID,
                    EX_SIDE_SET,
                    iSideSetInExoMesh + 1,
                    tSetLabel.c_str() );

        }    // end for: each side set

    }    // end function: Writer_Exodus::write_side_sets()

    //--------------------------------------------------------------------------------------------------------------
    // Static (for now)
    //--------------------------------------------------------------------------------------------------------------

    const char*
    Writer_Exodus::get_exodus_block_topology( CellTopology aCellTopology )
    {
        switch ( aCellTopology )
        {
            case CellTopology::TRI3:
                return "TRI3";
            case CellTopology::TRI6:
                return "TRI6";
            case CellTopology::QUAD4:
                return "QUAD4";
            case CellTopology::QUAD8:
                return "QUAD8";
            case CellTopology::QUAD9:
                return "QUAD9";
            case CellTopology::QUAD16:
                return "QUAD4";    // use 4-node QUAD as ParaView does not support 16-node QUAD
            case CellTopology::TET4:
                return "TET4";
            case CellTopology::TET10:
                return "TET10";
            case CellTopology::HEX8:
                return "HEX8";
            case CellTopology::HEX20:
                return "HEX20";
            case CellTopology::HEX27:
                return "HEX27";
            case CellTopology::HEX64:
                return "HEX8";    // use 8-node HEX as ParaView does not support 64-node HEX
            case CellTopology::PRISM6:
                return "PRISM6";
            default:
                MORIS_ERROR( false, "This element is invalid or it hasn't been implemented yet!" );
                return "";
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    int Writer_Exodus::get_nodes_per_element( CellTopology aCellTopology )
    {
        switch ( aCellTopology )
        {
            case CellTopology::TRI3:
                return 3;
            case CellTopology::TRI6:
                return 6;
            case CellTopology::QUAD4:
                return 4;
            case CellTopology::QUAD8:
                return 8;
            case CellTopology::QUAD9:
                return 9;
            case CellTopology::QUAD16:
                return 4;    // use 4-node QUAD as ParaView does not support 16-node QUAD
            case CellTopology::TET4:
                return 4;
            case CellTopology::TET10:
                return 10;
            case CellTopology::HEX8:
                return 8;
            case CellTopology::HEX20:
                return 20;
            case CellTopology::HEX27:
                return 27;
            case CellTopology::HEX64:
                return 8;    // use 8-node HEX as ParaView does not support 64-node HEX
            case CellTopology::PRISM6:
                return 6;
            default:
                MORIS_ERROR( false, "This element is invalid or it hasn't been implemented yet!" );
                return 0;
        }
    }

    //--------------------------------------------------------------------------

    void Writer_Exodus::create_directory( const std::string& aDirectoryName )
    {
        namespace fs = std::filesystem;

        // skip if directory name is empty
        if ( aDirectoryName.empty() )
        {
            return;
        }

        // Define folder path
        fs::path tDirectoryPath( aDirectoryName );

        // Check if path exists
        if ( !fs::exists( tDirectoryPath ) )
        {

            // Create all necessary directories
            MORIS_ERROR( fs::create_directories( tDirectoryPath ),
                    "create_director - failed to create %s",
                    aDirectoryName.c_str() );
        }
    }

    //--------------------------------------------------------------------------

    bool Writer_Exodus::isFileNameOnly( const std::string& aFileName )
    {
        if ( aFileName.rfind( "./", 0 ) == 0 )
        {
            return aFileName.find( '/', 2 ) == std::string::npos;
        }
        return aFileName.find( '/' ) == std::string::npos;
    }

    //--------------------------------------------------------------------------

}    // namespace moris::mtk
