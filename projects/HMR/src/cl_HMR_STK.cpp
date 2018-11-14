/*
 * cl_HMR_MTK.cpp
 *
 *  Created on: May 16, 2018
 *      Author: messe
 */

#include "fn_trans.hpp"
#include "fn_sort.hpp" //LINALG/src
#include "fn_print.hpp" //LINALG/src
#include "cl_MTK_Mesh.hpp" //MTK/src
#include "cl_HMR_STK.hpp" //HMR/src

#include "HMR_Tools.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp"  //HMR/src
#include "stk_impl/cl_MTK_Mesh_STK.hpp"

namespace moris
{
    namespace hmr
    {

// ----------------------------------------------------------------------------

    STK::STK( Lagrange_Mesh_Base * aMesh ) :
            mParameters( aMesh->get_parameters() ),
            mMesh( aMesh ),
            mNumberOfDimensions( aMesh->get_parameters()->get_number_of_dimensions() )
    {

    }

// ----------------------------------------------------------------------------
    void
    STK::create_mesh_data( const double aTimeStep )
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

        // get nunmber of fields
        uint tNumberOfFields = mMesh->get_number_of_fields();

        // reset field info
        mFieldsInfo.FieldsName.clear();
        mFieldsInfo.FieldsRank.clear();

        // first field is always element level
        mFieldsInfo.FieldsName.push_back( mMesh->get_field_label( 0 ) );
        mFieldsInfo.FieldsRank.push_back( EntityRank::ELEMENT );

        // second field is alwase element owner
        mFieldsInfo.FieldsName.push_back( mMesh->get_field_label( 1 ) );
        mFieldsInfo.FieldsRank.push_back( EntityRank::ELEMENT );

        // add nodal fields
        for( uint f=2; f<tNumberOfFields; ++f )
        {
            mFieldsInfo.FieldsName.push_back( mMesh->get_field_label( f ) );
            mFieldsInfo.FieldsRank.push_back( EntityRank::NODE );
        }

        // create new matrix with element data
        Matrix< DDRMat > & tElementLevels = mMesh->get_field_data( 0 );
        Matrix< DDRMat > & tElementOwners = mMesh->get_field_data( 1 );
        tElementLevels.set_size( tNumberOfElements, 1 );
        tElementOwners.set_size( tNumberOfElements, 1 );

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
                mElementTopology( e, k ) = tNodeIDs( k );
            }

            // save element index in map
            mElementLocalToGlobal( e ) = tElement->get_id();

            // save level of element
            tElementLevels( e ) = tElement->get_level();

            // save owners of element
            tElementOwners( e ) = tElement->get_owner();
        }

        // field 1 is always vertex ids
        Matrix< DDRMat > & tVertexIDs = mMesh->get_field_data( 2 );
        tVertexIDs.set_size( tNumberOfNodes, 1 );

        for( uint k=0; k<tNumberOfNodes; ++k )
        {
            auto tNode = mMesh->get_node_by_index( k );

            Matrix< DDRMat > tNodeCoords = tNode->get_coords();

            // copy coords to output matrix
            mNodeCoords.set_row( k , tNodeCoords );

            // copy node Owner
            mNodeOwner( k ) = tNode->get_owner();

            // copy node index into map
            mNodeLocalToGlobal( k ) = tNode->get_id();

            // save vertex id
            tVertexIDs( k ) = tNode->get_id();
        }

        // link mesh data object
        mMeshData.SpatialDim              = & mNumberOfDimensions;
        mMeshData.ElemConn(0)             = & mElementTopology;
        mMeshData.NodeCoords              = & mNodeCoords;
        mMeshData.EntProcOwner            = & mNodeOwner;
        mMeshData.LocaltoGlobalElemMap(0) = & mElementLocalToGlobal;
        mMeshData.LocaltoGlobalNodeMap    = & mNodeLocalToGlobal;
        mMeshData.FieldsInfo              = & mFieldsInfo;
        mFieldsInfo.FieldsData            = & mMesh->get_field_data();
        mMeshData.SetsInfo                = & mSetsInfo;

        /* if( par_rank() == 1 )
        {
            print( mElementTopology, "topo" );
            print( mElementLocalToGlobal, "Elements" );
            print( mNodeLocalToGlobal, "Nodes" );
            print( mNodeCoords, "Coords" );
            print( mNodeOwner, "owner" );

            moris::Cell< Matrix< DDRMat > > & tData = *mFieldsInfo.FieldsData;

            for( uint k=0; k<mFieldsInfo.FieldsData->size(); ++k )
            {
                std::string tRank;
                if( mFieldsInfo.FieldsRank( k ) == EntityRank::ELEMENT )
                {
                    tRank = " e ";
                }
                else if( mFieldsInfo.FieldsRank( k ) == EntityRank::NODE )
                {
                    tRank = " n ";
                }

                std::cout << k << " " << mFieldsInfo.FieldsName( k ) << tRank
                        << tData( k ).length() << std::endl;
            }
        } */

        // set timestep of mesh data object
        mMeshData.TimeStamp = aTimeStep;
        mMeshData.AutoAuraOptionInSTK = false;

        // get number of sets
        uint tNumberOfSideSets = mMesh->get_number_of_side_sets();

        // clear info table
        mSetsInfo.SideSetsInfo.clear();

        for( uint k=0; k<tNumberOfSideSets; ++k )
        {

            // get info
            mtk::MtkSideSetInfo & tInfo =  mMesh->get_side_set_info( k );

            // push back as pointer
            mSetsInfo.SideSetsInfo.push_back( &tInfo );
        }

        // special function for old mesh
        this->flag_old_and_new_elements();

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

    void
    STK::save_to_file( const std::string & aFilePath )
    {
        tic tTimer;

        // create database object
        moris::mtk::Mesh_STK tMesh( mMeshData );

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

        void
        STK::flag_old_and_new_elements()
        {
            uint tInputPattern  = mParameters->get_input_pattern();
            uint tOutputPattern = mParameters->get_output_pattern();

            if( mMesh->get_activation_pattern() ==  tInputPattern )
            {
                // get number of elements on mesh
                uint tNumberOfElements = mMesh->get_number_of_elements();

                uint tFieldIndex = mMesh->create_field_data( "Refinement"  );

                mFieldsInfo.FieldsName.push_back(  mMesh->get_field_label( tFieldIndex ) );
                mFieldsInfo.FieldsRank.push_back( EntityRank::ELEMENT );

                // link to field
                Matrix< DDRMat > & tData = mMesh->get_field_data()( tFieldIndex );

                tData.set_size( tNumberOfElements, 1 );

                // loop over all elements
                for( uint e=0; e<tNumberOfElements; ++e )
                {
                    // get background element
                    Background_Element_Base * tElement = mMesh->get_element( e )->get_background_element();

                    if( tElement->is_active( tOutputPattern ) )
                    {
                        tData( e ) = 0.0;
                    }
                    else if(  tElement->is_refined( tOutputPattern ) )
                    {
                        tData( e ) = 1.0;
                    }
                    else
                    {
                        tData( e ) = -1.0;
                    }
                }
            }
        }

    } /* namespace hmr */
} /* namespace moris */
