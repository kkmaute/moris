/*
 * cl_HMR_MTK.cpp
 *
 *  Created on: May 16, 2018
 *      Author: messe
 */

#include "linalg.hpp"
#include "fn_sort.hpp" //LNA/src
#include "cl_Mesh.hpp" //MTK/src
#include "cl_HMR_MTK.hpp" //HMR/src

#include "HMR_Tools.hpp" //HMR/src
namespace moris
{
    namespace hmr
    {

// ----------------------------------------------------------------------------

    void
    MTK::create_mesh_data(
            const luint                    & aNumberOfElements,
            const uint                     & aNumberOfNodesPerElement,
            const luint                    & aNumberOfNodes )
    {
        // start timer
        tic tTimer;

        // get my rank
        uint tMyRank = par_rank();

        // ask settings for number of dumensions
        mNumberOfDimensions = mParameters->get_number_of_dimensions();

        // initialize topology field
        mElementTopology.set_size( aNumberOfElements, aNumberOfNodesPerElement );

        // initialize local to global map
        mElementLocalToGlobal.set_size( aNumberOfElements, 1 );

        // initialize coordinates
        mNodeCoords.set_size( aNumberOfNodes, mNumberOfDimensions );

        // initialize local to global map for nodes
        mNodeLocalToGlobal.set_size( aNumberOfNodes, 1 );

        // initialize node ownership
        mNodeOwner.set_size( aNumberOfNodes, 1 );

        // initialize data fields
        Mat<real> tEmpty;
        mFieldData.clear();
        mFieldData.resize( 1, { tEmpty } );
        mFieldsInfo.FieldsName.clear();
        mFieldsInfo.FieldsRank.clear();

        // Field 0: Element Level
        mFieldsInfo.FieldsName.push_back( "Element_Level");
        mFieldsInfo.FieldsRank.push_back( EntityRank::ELEMENT );
        mFieldData( 0 ).set_size( aNumberOfElements, 1 );

        // initialize counter
        luint tCount = 0;

        // loop over all elements
        for( auto tElement : mAllElementsOnProc )
        {
            // test if element is active and owned
            if ( tElement->get_owner() == tMyRank
                    && tElement->is_active() )
            {
                // loop over all nodes
                for( uint k=0; k<aNumberOfNodesPerElement; ++k )
                {
                    // save domain index into topology
                    mElementTopology( tCount, k )
                            = tElement->get_basis( k )->get_domain_index() + 1;
                }

                // save element index in map
                mElementLocalToGlobal( tCount )
                    = tElement->get_domain_index() + 1;

                // save element level in field
                mFieldData( 0 )( tCount ) = tElement->get_level();

                // increment counter
                ++tCount;
            }
        }

        // reset counter
        tCount = 0;

        // loopover all nodes
        for( auto tNode : mAllBasisOnProc )
        {
            // test if node is used by any owned element
            if ( tNode->is_used() )
            {
                // get coordinate of node
                const real * tXZY = tNode->get_xyz();

                // copy node coordinates into matrix
                for( uint i=0; i<mNumberOfDimensions; ++i )
                {
                    mNodeCoords( tCount, i ) = tXZY[ i ];
                }

                // copy node Owner
                mNodeOwner( tCount ) = tNode->get_owner();

                // copy node index into map
                mNodeLocalToGlobal( tCount ) = tNode->get_domain_index() + 1;

                // increment counter
                ++tCount;
            }
        }

        // link mesh data object
        mMeshData.SpatialDim           = & mNumberOfDimensions;
        mMeshData.ElemConn             = & mElementTopology;
        mMeshData.NodeCoords           = & mNodeCoords;
        mMeshData.EntProcOwner         = & mNodeOwner;
        mMeshData.LocaltoGlobalElemMap = & mElementLocalToGlobal;
        mMeshData.LocaltoGlobalNodeMap = & mNodeLocalToGlobal;
        mMeshData.FieldsInfo           = & mFieldsInfo;
        mFieldsInfo.FieldsData         = & mFieldData;


        if ( mVerboseFlag )
        {
            // stop timer
            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

            // print output
            std::fprintf( stdout,"%s Created MTK output object.\n               Mesh has %lu elements and %lu nodes.\n               Creation took %5.3f seconds.\n\n",
                    proc_string().c_str(),
                    ( long unsigned int ) aNumberOfElements,
                    ( long unsigned int ) aNumberOfNodes,
                    ( double ) tElapsedTime / 1000);
        }

    }

// ----------------------------------------------------------------------------

    void
    MTK::add_node_data(
            const std::string & aLabel,
            const Mat< real > & aData )
    {
        mFieldsInfo.FieldsName.push_back( aLabel );
        mFieldsInfo.FieldsRank.push_back( EntityRank::NODE );

        // get entry
        uint tFieldNr = mFieldData.size();

        // create empty matrix
        Mat< real > tEmpty;
        mFieldData.push_back( tEmpty );

        // initialize memory
        mFieldData( tFieldNr ).set_size( mNodeLocalToGlobal.length(), 1 );

        // initialize counter
        luint tCount = 0;

        // get number of all nodes
        luint tNumberOfNodes = mAllBasisOnProc.size();

        // loop over all nodes
        for( luint k=0; k<tNumberOfNodes; ++k )
        {
            if ( mAllBasisOnProc( k )->is_used() )
            {
                mFieldData( tFieldNr )( tCount++ ) = aData( k );
            }
        }
    }

// ----------------------------------------------------------------------------

    void
    MTK::add_element_data(
            const std::string & aLabel,
            const Mat< real > & aData )
    {
        mFieldsInfo.FieldsName.push_back( aLabel );
        mFieldsInfo.FieldsRank.push_back( EntityRank::ELEMENT );

        // get entry
        uint tFieldNr = mFieldData.size();

        // create empty matrix
        Mat< real > tEmpty;
        mFieldData.push_back( tEmpty );

        // initialize memory
        mFieldData( tFieldNr ).set_size( mElementLocalToGlobal.length(), 1 );

        // initialize counter
        luint tCount = 0;

        uint tMyRank = par_rank();

        // get number of all elements
        luint tNumberOfElements = mAllElementsOnProc.size();

        // loop over all elements
        for( luint e=0; e<tNumberOfElements; ++e  )
        {
            if ( mAllElementsOnProc( e )->get_owner() == tMyRank )
            {
                mFieldData( tFieldNr )( tCount++ ) = aData( e );
            }
        }
    }

// ----------------------------------------------------------------------------

    void
    MTK::save_to_file( const std::string & aFilePath )
    {

        tic tTimer;

        // create database object
        moris::mesh tMesh( MeshType::MTK, mMeshData );

        // copy file path, since tMesh does not like const input
        std::string tFilePath = aFilePath;

        // save file
        tMesh.create_output_mesh( tFilePath );

        if ( mVerboseFlag )
        {
            // stop timer
            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

            // print output
            std::fprintf( stdout,"%s Wrote MTK mesh to file.\n               Writing took %5.3f seconds.\n\n",
                    proc_string().c_str(),
                    ( double ) tElapsedTime / 1000 );
        }

    }

// ----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */
