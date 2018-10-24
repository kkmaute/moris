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
        mFieldsInfo.clear_fields();

        // Initialize scalar field data
        mFields = moris::Cell<mtk::Scalar_Field_Info<DDRMat>>(tNumberOfFields+2);

        // first field is always element level
        mFields(0).set_field_name(mMesh->get_field_label( 0 ));
        mFields(0).set_field_entity_rank(EntityRank::ELEMENT);
        mFields(0).add_field_data(&mElementLocalToGlobal, & mMesh->get_field_data()(0));
        mFieldsInfo.mRealScalarFields.push_back(&mFields(0));

        // Second field is always element owner
        mFields(1).set_field_name(mMesh->get_field_label( 1 ));
        mFields(1).set_field_entity_rank(EntityRank::ELEMENT);
        mFields(1).add_field_data(&mElementLocalToGlobal,& mMesh->get_field_data()(1));
        mFieldsInfo.mRealScalarFields.push_back(&mFields(1));

        // add nodal fields
        for( uint f=2; f<tNumberOfFields; ++f )
        {

            mFields(f).set_field_name(mMesh->get_field_label( f ));
            mFields(f).set_field_entity_rank(EntityRank::NODE);
            mFields(f).add_field_data(&mNodeLocalToGlobal,&mMesh->get_field_data()(f));
            mFieldsInfo.mRealScalarFields.push_back(&mFields(f));
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

    } /* namespace hmr */
} /* namespace moris */
