/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Reader_Exodus.cpp
 *
 */

#include <exodusII.h>
#include "cl_MTK_Reader_Exodus.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"

// TODO
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Sets_Info.hpp"

#include <iostream>
#include <utility>

namespace moris::mtk
{

    Reader_Exodus::Reader_Exodus()
    {
        mExoID = 0;
    }

    Reader_Exodus::~Reader_Exodus() = default;

    void Reader_Exodus::set_error_options( bool abort, bool debug, bool verbose )
    {
        ex_opts( abort * EX_ABORT | debug * EX_DEBUG | verbose * EX_VERBOSE );
    }

    void Reader_Exodus::open_file( const std::string& aExodusFileName, float aVersion )
    {
        int tCPUWordSize = 8, tIOWordSize = 8;    // TODO
        mExoID = ex_open( aExodusFileName.c_str(), EX_READ, &tCPUWordSize, &tIOWordSize, &aVersion );
    }

    void Reader_Exodus::close_file()
    {
        ex_close( mExoID );
    }

    void Reader_Exodus::read_file( const std::string& aFileName )
    {
        Reader_Exodus::open_file( aFileName );

        // Get init from Exodus
        char        tDatabaseName[ 100 ] = "";    // FIXME
        moris::uint tNumNodes, tNumElements, tNumBlocks, tNumNodeSets, tNumSideSets;
        ex_get_init( mExoID, tDatabaseName, &mNumSpatialDimensions, &tNumNodes, &tNumElements, &tNumBlocks, &tNumNodeSets, &tNumSideSets );

        // Create and population mesh data input
        moris::mtk::MtkMeshData mMeshDataInput( tNumBlocks );
        mMeshDataInput.SpatialDim             = &mNumSpatialDimensions;
        mMeshDataInput.MarkNoBlockForIO       = false;
        mMeshDataInput.CreateAllEdgesAndFaces = true;
        mMeshDataInput.Verbose                = true;

        // Read the coordinates from the exodus file
        moris::Matrix< moris::DDRMat > tXCoordinates( tNumNodes, 1, 0.0 );
        moris::Matrix< moris::DDRMat > tYCoordinates( tNumNodes, 1, 0.0 );
        moris::Matrix< moris::DDRMat > tZCoordinates( tNumNodes, 1, 0.0 );
        ex_get_coord( mExoID, tXCoordinates.data(), tYCoordinates.data(), tZCoordinates.data() );

        // Combine the coordinates into one matrix for input
        mNumSpatialDimensions = 3;    // FIXME
        mNodeCoordinates.resize( tNumNodes, mNumSpatialDimensions );
        mNodeCoordinates.set_column( 0, tXCoordinates );
        mNodeCoordinates.set_column( 1, tYCoordinates );
        mNodeCoordinates.set_column( 2, tZCoordinates );
        mMeshDataInput.NodeCoords = &mNodeCoordinates;

        // Get the node map
        mNodeMap.resize( tNumNodes, 1 );
        ex_get_id_map( mExoID, EX_NODE_MAP, mNodeMap.data() );
        mMeshDataInput.LocaltoGlobalNodeMap = &mNodeMap;

        // Node ownership
        mNodeOwner.set_size( 1, tNumNodes, moris::par_rank() );
        mMeshDataInput.NodeProcOwner = &mNodeOwner;

        // Prepare for get_block
        mBlockSetInfo.resize( tNumBlocks, moris::mtk::MtkBlockSetInfo() );
        mBlockDescription.resize( tNumBlocks, "" );
        mElemConn.resize( tNumBlocks, moris::Matrix< moris::IndexMat >( 1, 1 ) );
        mLocaltoGlobalElemMap.resize( tNumBlocks, moris::Matrix< moris::IdMat >( 1, 1 ) );
        moris::uint tNumElementsInBlock, tNumNodesPerElement, tNumEdgesPerElement, tNumFacesPerElement, tNumAttributesPerElement;
        char        tReadName[ 100 ] = "";

        // Loop through all of the blocks
        for ( moris::uint tBlockIndex = 0; tBlockIndex < tNumBlocks; tBlockIndex++ )
        {
            // Get the block parameters
            ex_get_block( mExoID, EX_ELEM_BLOCK, tBlockIndex + 1, tReadName, &tNumElementsInBlock, &tNumNodesPerElement, &tNumEdgesPerElement, &tNumFacesPerElement, &tNumAttributesPerElement );

            // Cell topology
            mBlockSetInfo( tBlockIndex ).mBlockSetTopo = this->get_cell_topology( tReadName );

            // Block Name
            ex_get_name( mExoID, EX_ELEM_BLOCK, tBlockIndex + 1, tReadName );
            mBlockDescription( tBlockIndex ).assign( tReadName, 100 );
            mBlockSetInfo( tBlockIndex ).mBlockSetName = mBlockDescription( tBlockIndex );

            // Get connectivity
            moris::Matrix< moris::IndexMat > tBlockNodesVector( tNumElementsInBlock * tNumNodesPerElement, 1, 0 );
            ex_get_conn( mExoID, EX_ELEM_BLOCK, tBlockIndex + 1, tBlockNodesVector.data(), nullptr, nullptr );

            // Reorganize connectivity
            mElemConn( tBlockIndex ).resize( tNumElementsInBlock, tNumNodesPerElement );
            moris::uint tVectorIndex = 0;
            for ( moris::uint tElementIndex = 0; tElementIndex < tNumElementsInBlock; tElementIndex++ )
            {
                for ( moris::uint tNodeIndex = 0; tNodeIndex < tNumNodesPerElement; tNodeIndex++ )
                {
                    mElemConn( tBlockIndex )( tElementIndex, tNodeIndex ) = tBlockNodesVector( tVectorIndex );
                    tVectorIndex++;
                }
            }
            mMeshDataInput.ElemConn( tBlockIndex ) = &mElemConn( tBlockIndex );

            // Element map
            mLocaltoGlobalElemMap( tBlockIndex ).set_size( tNumElementsInBlock, 1 );
            ex_get_id_map( mExoID, EX_ELEM_MAP, mLocaltoGlobalElemMap( tBlockIndex ).data() );    // FIXME
            mMeshDataInput.LocaltoGlobalElemMap( tBlockIndex ) = &mLocaltoGlobalElemMap( tBlockIndex );

            // Cell ids
            mBlockSetInfo( tBlockIndex ).mCellIdsInSet = mMeshDataInput.LocaltoGlobalElemMap( tBlockIndex );

            // Add sets to input data
            mMtkMeshSets.add_block_set( &mBlockSetInfo( tBlockIndex ) );
        }
        mMeshDataInput.SetsInfo = &mMtkMeshSets;

        // Side sets
        moris::uint tNumElementsInSideSet;
        mSideSetInfo.resize( tNumSideSets, moris::mtk::MtkSideSetInfo() );
        mSideSetElements.resize( tNumSideSets, moris::Matrix< moris::IndexMat >( 1, 1 ) );
        mSideSetNames.resize( tNumSideSets, "" );
        moris::Matrix< moris::IndexMat > tSideSetOrdinals;

        for ( moris::uint tSideSetIndex = 0; tSideSetIndex < tNumSideSets; tSideSetIndex++ )
        {
            // Get number of entries in this side set and resize matrices
            ex_get_set_param( mExoID, EX_SIDE_SET, tSideSetIndex + 1, &tNumElementsInSideSet, nullptr );
            mSideSetElements( tSideSetIndex ).set_size( tNumElementsInSideSet, 2 );
            tSideSetOrdinals.set_size( tNumElementsInSideSet, 1 );

            // Get data
            ex_get_set( mExoID, EX_SIDE_SET, tSideSetIndex + 1, mSideSetElements( tSideSetIndex ).data(), tSideSetOrdinals.data() );
            ex_get_name( mExoID, EX_SIDE_SET, tSideSetIndex + 1, tReadName );

            // Adjust for moris
            mSideSetNames( tSideSetIndex ).assign( tReadName );
            for ( moris::uint tElementNum = 0; tElementNum < tNumElementsInSideSet; tElementNum++ )
            {
                mSideSetElements( tSideSetIndex )( tElementNum, 1 )--;
                tSideSetOrdinals( tElementNum )--;
            }
            mSideSetElements( tSideSetIndex ).set_column( 1, tSideSetOrdinals );

            // Assign data
            mSideSetInfo( tSideSetIndex ).mElemIdsAndSideOrds = &mSideSetElements( tSideSetIndex );
            mSideSetInfo( tSideSetIndex ).mSideSetName        = mSideSetNames( tSideSetIndex );
            mMtkMeshSets.add_side_set( &mSideSetInfo( tSideSetIndex ) );
        }

        mMesh = moris::mtk::create_integration_mesh( moris::mtk::MeshType::STK, mMeshDataInput );

        Reader_Exodus::close_file();
    }

    // ---------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------

    moris::mtk::CellTopology Reader_Exodus::get_cell_topology( char* aElementType )
    {
        std::string tElementType( aElementType );
        if ( tElementType == "tri3" )
            return moris::mtk::CellTopology::TRI3;
        else if ( tElementType == "quad4" )
            return moris::mtk::CellTopology::QUAD4;
        else if ( tElementType == "tet4" )
            return moris::mtk::CellTopology::TET4;
        else if ( tElementType == "tet10" )
            return moris::mtk::CellTopology::TET10;
        else if ( tElementType == "hex8" )
            return moris::mtk::CellTopology::HEX8;
        else if ( tElementType == "prism6" )
            return moris::mtk::CellTopology::PRISM6;
        else
            return moris::mtk::CellTopology::UNDEFINED;
    }
}    // namespace moris::mtk
