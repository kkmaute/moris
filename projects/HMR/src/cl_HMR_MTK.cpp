/*
 * cl_HMR_MTK.cpp
 *
 *  Created on: May 16, 2018
 *      Author: messe
 */

#include "linalg.hpp"
#include "fn_trans.hpp"
#include "fn_sort.hpp" //LNA/src
#include "cl_Mesh.hpp" //MTK/src
#include "cl_HMR_MTK.hpp" //HMR/src

#include "HMR_Tools.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp"  //HMR/src

namespace moris
{
    namespace hmr
    {

// ----------------------------------------------------------------------------

    MTK::MTK( Lagrange_Mesh_Base * aMesh ) :
            mParameters( aMesh->get_parameters() ),
            mMesh( aMesh ),
            mNumberOfDimensions( aMesh->get_parameters()->get_number_of_dimensions() )
    {

    }

// ----------------------------------------------------------------------------
    void
    MTK::create_mesh_data()
    {
        // start timer
        tic tTimer;

        // activate this pattern on background mesh
        mMesh->select_activation_pattern();

        // get number of elements
        auto tNumberOfElements        = mMesh->get_number_of_elements();

        // get number of nodes
        auto tNumberOfNodes           = mMesh->get_number_of_nodes_on_proc();

        // get number of nodes per element
        auto tNumberOfNodesPerElement = mMesh->get_number_of_basis_per_element();

        // initialize topology field
        mElementTopology.set_size( tNumberOfElements, tNumberOfNodesPerElement );

        // initialize local to global map
        mElementLocalToGlobal.set_size( tNumberOfElements, 1 );

        // initialize coordinates
        mNodeCoords.set_size( tNumberOfNodes, mNumberOfDimensions );

        // initialize local to global map for nodes
        mNodeLocalToGlobal.set_size( tNumberOfNodes, 1 );

        // initialize node ownership
        mNodeOwner.set_size( tNumberOfNodes, 1 );


        // get nunber of fields
        uint tNumberOfFields = mMesh->get_number_of_fields();

        // reset field info
        mFieldsInfo.FieldsName.clear();
        mFieldsInfo.FieldsRank.clear();

        // first field is always element level
        mFieldsInfo.FieldsName.push_back( mMesh->get_field_label( 0 ) );
        mFieldsInfo.FieldsRank.push_back( EntityRank::ELEMENT );

        // add nodal fields
        for( uint f=1; f<tNumberOfFields; ++f )
        {
            mFieldsInfo.FieldsName.push_back( mMesh->get_field_label( f ) );
            mFieldsInfo.FieldsRank.push_back( EntityRank::NODE );
        }

        // create new matrix with element data
        //Mat< real > tElementLevels( tNumberOfElements, 1 );
        Mat< real > & tElementLevels = mMesh->get_field_data( 0 );
        tElementLevels.set_size( tNumberOfElements, 1 );

        // loop over all elements
        for( uint e=0; e<tNumberOfElements; ++e )
        {
            // get pointer to element
            auto tElement = mMesh->get_element( e );

            MORIS_ASSERT( tElement->is_active(), "Tried to create inactive element" );
            // get node IDs
            auto tNodeIDs = tElement->get_vertex_ids();

            // cast copy node IDs to topology matrix
            for( uint k=0; k<tNumberOfNodesPerElement; ++k )
            {
                mElementTopology( e, k ) = tNodeIDs( k ) + 1;
            }

            // save element index in map
            mElementLocalToGlobal( e ) = tElement->get_domain_index() + 1;

            // save level of element
            tElementLevels( e ) = tElement->get_level();
        }

        Mat< real > & tVertexIDs = mMesh->get_field_data( 1 );
        tVertexIDs.set_size( tNumberOfNodes, 1 );

        for( uint k=0; k<tNumberOfNodes; ++k )
        {
            auto tNode = mMesh->get_node_by_index( k );

            // get coordinates of node
            Mat< real > tCoords = trans( tNode->get_coords() );

            // copy coords to output matrix
            mNodeCoords.rows( k, k ) = tCoords.rows( 0, 0 );

            // copy node Owner
            mNodeOwner( k ) = tNode->get_owner();

            // copy node index into map
            mNodeLocalToGlobal( k ) = tNode->get_domain_index() + 1;

            // save vertex id
            tVertexIDs( k ) = tNode->get_id();
        }

        // link mesh data object
        mMeshData.SpatialDim           = & mNumberOfDimensions;
        mMeshData.ElemConn             = & mElementTopology;
        mMeshData.NodeCoords           = & mNodeCoords;
        mMeshData.EntProcOwner         = & mNodeOwner;
        mMeshData.LocaltoGlobalElemMap = & mElementLocalToGlobal;
        mMeshData.LocaltoGlobalNodeMap = & mNodeLocalToGlobal;
        mMeshData.FieldsInfo           = & mFieldsInfo;
        mFieldsInfo.FieldsData         =   mMesh->get_field_data();

        if ( mParameters->is_verbose() )
        {
            // stop timer
            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

            // print output
            std::fprintf( stdout,"%s Created MTK output object.\n               Mesh has %lu elements and %lu nodes.\n               Creation took %5.3f seconds.\n\n",
                    proc_string().c_str(),
                    ( long unsigned int ) tNumberOfElements,
                    ( long unsigned int ) tNumberOfNodes,
                    ( double ) tElapsedTime / 1000);
        }

    }

// ----------------------------------------------------------------------------

/*    void
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
    } */

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

        if ( mParameters->is_verbose() )
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
