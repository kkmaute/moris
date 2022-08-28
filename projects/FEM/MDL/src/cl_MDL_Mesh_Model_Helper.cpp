/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MDL_Mesh_Model_Helper.cpp
 *
 */

#include "cl_Stopwatch.hpp" //CHR/src

#include "fn_unique.hpp" //CHR/src

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"                    //MTK/src

#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src
#include "cl_FEM_Enums.hpp"               //FEM/INT/src

#include "cl_FEM_Element_Bulk.hpp"               //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Element_Factory.hpp"
#include "cl_FEM_Set.hpp"

#include "cl_MDL_Mesh_Model_Helper.hpp"

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

    Mesh_Model_Helper::Mesh_Model_Helper(       mtk::Mesh_Manager * aMesh,
                                          const moris_index         aMeshPairIndex ) : mMesh( aMesh )
    {
        // get mesh

        // get colors for blocks
        // colormap build on first entry of list
        // create fem::blocks -- input
        //

        // get blocks
        // get vertices on blocks
        // make them unique
        // build nodes based on unique list
        //

        // list of colors with corresponding blocks

        mMesh->get_mesh_pair( aMeshPairIndex, mInterpolationMesh, mIntegrationMesh );

        //========================================================================
        // Create num vertices on blocks list
        uint tNumBlocks = mIntegrationMesh->get_num_blocks();
        mVerticesOnBlock.set_size( tNumBlocks, 1 , -1 );

        // Loop over the blocks
        for( uint Ik = 0; Ik < tNumBlocks; Ik++ )
        {
            mVerticesOnBlock( Ik ) = mIntegrationMesh->get_block_by_index( Ik )->get_num_vertices_on_set( true );
        }

        //==========================================================================
        // create num vertices on side set list
        uint tNumSideSets = mIntegrationMesh->get_num_side_set();
        mVerticesOnSideSet.set_size( tNumSideSets, 1 , -1 );

//        std::cout<<tNumSideSets<<std::endl;

        // Loop over the side sets
        for( uint Ik = 0; Ik < tNumSideSets; Ik++ )
        {
            mVerticesOnSideSet( Ik ) = mIntegrationMesh->get_side_set_by_index( Ik )->get_num_vertices_on_set( true);
        }

        MORIS_ASSERT( mVerticesOnBlock.min() != -1, "negative number of vertices on block");
        MORIS_ASSERT( mVerticesOnSideSet.min() != -1, "negative number of vertices on side set");
    }

//------------------------------------------------------------------------------

    Mesh_Model_Helper::~Mesh_Model_Helper()
    {}

//------------------------------------------------------------------------------

    void Mesh_Model_Helper::set_block_by_color( const uint             tColor,
                                                const Matrix< DDUMat > tBlockIndex )
    {
        if ( mColorListBlock.size() < tColor )
        {
            mColorListBlock.resize( tColor + 1 );
        }

        MORIS_ASSERT( mColorListBlock( tColor ).numel() == 0, "Mesh_Model_Helper::set_block_by_color(): Blocks for this color were set earlier");

        mColorListBlock( tColor ) = tBlockIndex;
    }

//------------------------------------------------------------------------------

    void Mesh_Model_Helper::set_side_set_by_color( const uint             tColor,
                                                   const Matrix< DDUMat > tSideSetIndex )
    {
        if ( mColorListSideSet.size() < tColor )
        {
            mColorListSideSet.resize( tColor + 1 );
        }

        MORIS_ASSERT( mColorListSideSet( tColor ).numel() == 0, "Mesh_Model_Helper::set_block_by_color(): Blocks for this color were set earlier");

        mColorListSideSet( tColor ) = tSideSetIndex;
    }

//------------------------------------------------------------------------------

    void Mesh_Model_Helper::compute_unique_node_lists()
    {
        // start timer
        tic tTimer;

        // Set Cell containing vertex ind to color size
        VertexIndOnColor.resize( mColorListBlock.size() );

        Matrix< DDSMat > tNumMaxVertPerColor( mColorListBlock.size(), 1 , -1 );

        this->compute_max_num_vertices_on_color( tNumMaxVertPerColor );

        this->compute_unique_vertex_list_on_color( tNumMaxVertPerColor );

        uint tNumNodes = 0;
        for( uint Ik = 0; Ik < VertexIndOnColor.size(); Ik++ )
        {
            tNumNodes = tNumNodes + VertexIndOnColor( Ik ).numel();
        }

        mNodes.resize( tNumNodes, nullptr );
        mNodeToVertexIndMap.set_size( tNumNodes, 1, -1 );

        //print(VertexIndOnColor, "VertexIndOnColor");

        //==================================================================================

        mVertexColorToNodeIndMap.resize(mColorListBlock.size());
        for( uint Ii = 0; Ii<mColorListBlock.size(); ++Ii )
        {
            //FIXME
            mVertexColorToNodeIndMap( Ii ).set_size( 1000000, 1, -1);
        }

        uint tCounter3 = 0;
        for( uint Ii = 0; Ii<VertexIndOnColor.size(); ++Ii )
        {
            for( uint Ia = 0; Ia<VertexIndOnColor( Ii ).numel(); ++Ia )
            {
                mNodes( tCounter3 ) = new fem::Node( &mIntegrationMesh->get_mtk_vertex( VertexIndOnColor( Ii )( Ia ) ) );

                mNodes( tCounter3 )->set_index( tCounter3 );

                mVertexColorToNodeIndMap( Ii )( VertexIndOnColor( Ii )( Ia ) = tCounter3 );

                mNodeToVertexIndMap( tCounter3++, 0 ) = VertexIndOnColor( Ii )( Ia );
            }
        }

        this->assign_node_ids();

        if( par_rank() == 0)
        {
            // stop timer
            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

            // print output
            std::fprintf( stdout,"Model: created %u FEM nodes in %5.3f seconds.\n\n",
                    ( unsigned int ) mNodes.size(),
                    ( double ) tElapsedTime / 1000 );
        }
    }

//------------------------------------------------------------------------------

    void Mesh_Model_Helper::compute_max_num_vertices_on_color( Matrix< DDSMat > & aNumMaxVertPerColor )
    {
        // FIXME Check if block and side sets have the same colors and ordering
        // Loop over the block colors
        for( uint Ik = 0; Ik < mColorListBlock.size(); Ik++ )
        {
            uint tCounter = 0;

            // Loop over the blocks
            for( uint Ij = 0; Ij < mColorListBlock( Ik ).numel(); Ij++ )
            {
                tCounter = tCounter + mVerticesOnBlock( mColorListBlock( Ik )( Ij, 0 ), 0 );
            }

            // Loop over the side sets
            for( uint Ij = 0; Ij < mColorListSideSet( Ik ).numel(); Ij++ )
            {
                tCounter = tCounter + mVerticesOnSideSet( mColorListSideSet( Ik )( Ij, 0 ), 0 );
            }

            aNumMaxVertPerColor( Ik, 0 ) = tCounter;
        }
    }

//------------------------------------------------------------------------------

    void Mesh_Model_Helper::compute_unique_vertex_list_on_color( const Matrix< DDSMat > & aNumMaxVertPerColor )
    {
        // FIXME Check if block and side sets have the same colors and ordering
        // Loop over the block colors
        for( uint Ik = 0; Ik < mColorListBlock.size(); Ik++ )
        {
            Matrix< DDSMat >tVertInds( 1,aNumMaxVertPerColor( Ik ,0 ) );
//            Matrix< DDSMat >tVertInds( aNumMaxVertPerColor( Ik ,0 ), 1 );
            uint tCounter = 0;

            // Loop over the blocks
            for( uint Ij = 0; Ij < mColorListBlock( Ik ).length(); Ij++ )
            {
                //FIXME rewrite for more readability
                tVertInds( { 0, 0 }, { tCounter, tCounter + mVerticesOnBlock( mColorListBlock( Ik )( Ij, 0 ), 0 ) - 1 } ) =
                        mIntegrationMesh->get_block_by_index( mColorListBlock( Ik )( Ij, 0 ) )->get_ig_vertices_inds_on_block( true ).matrix_data();

                tCounter = tCounter + mVerticesOnBlock( mColorListBlock( Ik )( Ij, 0 ), 0 );
            }

            // Loop over the side sets
            for( uint Ij = 0; Ij < mColorListSideSet( Ik ).length(); Ij++ )
            {
                //FIXME rewrite for more readability
                tVertInds( { 0, 0 }, { tCounter, tCounter + mVerticesOnSideSet( mColorListSideSet( Ik )( Ij, 0 ), 0 ) - 1 } ) =
                        mIntegrationMesh->get_block_by_index( mColorListSideSet( Ik )( Ij, 0 ) )->get_ig_vertices_inds_on_block( true ).matrix_data();

                tCounter = tCounter + mVerticesOnSideSet( mColorListSideSet( Ik )( Ij, 0 ), 0 );
            }

            unique( tVertInds, VertexIndOnColor( Ik ) );
        }
    }

//------------------------------------------------------------------------------

    void Mesh_Model_Helper::assign_node_ids()
    {
//        // Get list containing the number of owned nodes of each processor
//        Matrix< DDUMat > tNumOwnedNodeList;
//        comm_gather_and_broadcast( aNumOwnedAdofs, tNumOwnedNodeList );
//
//        Matrix< DDUMat > tOwnedANodesOffsetList( tNumOwnedNodeList.length(), 1, 0 );
//
//        // Loop over all entries to create the offsets. Starting with 1
//        for ( moris::uint Ij = 1; Ij < tOwnedANodesOffsetList.length(); Ij++ )
//        {
//            // Add the number of owned adofs of the previous processor to the offset of the previous processor
//            tOwnedANodesOffsetList( Ij, 0 ) = tOwnedANodesOffsetList( Ij-1, 0 ) + tNumOwnedNodeList( Ij-1, 0 );
//        }
//
//        return tOwnedANodesOffsetList( par_rank(), 0);

        for( uint Ii = 0; Ii<mNodes.size(); ++Ii )
        {
            mNodes( Ii )->set_id( Ii );
        }
    }

//------------------------------------------------------------------------------

    void Mesh_Model_Helper::create_node_set( const moris::uint               aColor,
                                             const Matrix< moris::IndexMat > aNodeSetIndex)
    {
        Cell< fem::Node_Base* > tNodeSet( aNodeSetIndex.numel(), nullptr );

        for( uint Ii = 0; Ii<aNodeSetIndex.numel(); ++Ii )
        {
            MORIS_ASSERT( mVertexColorToNodeIndMap( aColor )( aNodeSetIndex( Ii, 0 )) != -1,
                    " Mesh_Model_Helper::node_set(), Requested node does not exist");

            moris::uint tNodeIndex = mVertexColorToNodeIndMap( aColor )( aNodeSetIndex( Ii, 0 ));

            tNodeSet( Ii ) = mNodes( tNodeIndex );
        }

        mNodeSets.push_back( tNodeSet );
    }
//------------------------------------------------------------------------------

    void Mesh_Model_Helper::create_elements()
    {
//        // start timer
//        tic tTimer;
//
//        // a factory to create the elements
//        fem::Element_Factory tElementFactory;
//
//        moris::Cell< moris_index > tSideSetIndex( mColorListSideSet( 0 ).numel() );
//        for( luint Ii = 0; Ii<mColorListSideSet( 0 ).numel(); ++Ii)
//        {
//            tSideSetIndex( Ii ) = mColorListSideSet( 0 )( Ii );
//        }
//
//        // get the number of element to create
//        luint tNumberOfEquationObjects = mIntegrationMesh->get_num_elems()
//                                       + mIntegrationMesh->get_sidesets_num_faces( tSideSetIndex );     //Fixme
//
//        // create equation objects
//        mElements.reserve( tNumberOfEquationObjects );
//
//        uint tNumberElementBlocks = 1 + mColorListSideSet( 0 ).numel();
//
//        mElementBlocks.resize( tNumberElementBlocks, nullptr );
//
//        //  Create Blockset Elements ---------------------------------------------------
//        std::cout<<" Create Blockset Elements "<<std::endl;
//        //------------------------------------------------------------------------------
//
//        uint tNumBlocks = mIntegrationMesh->get_num_blocks();

    }

//    void Mesh_Model_Helper::create_nodes()
//    {}

    } /* namespace mdl */
} /* namespace moris */

