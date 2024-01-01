/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_STK.cpp
 *
 */

#include "cl_HMR_STK.hpp" //HMR/src

#include "cl_HMR_Lagrange_Mesh_Base.hpp"  //HMR/src
#include "HMR_Tools.hpp" //HMR/src
#include "fn_trans.hpp"
#include "fn_sort.hpp" //LINALG/src
#include "cl_MTK_Mesh.hpp" //MTK/src
#include "HMR_Tools.hpp" //HMR/src
#include "stk_impl/cl_MTK_Mesh_Core_STK.hpp"

namespace moris::hmr
{

    // ----------------------------------------------------------------------------

    STK::STK( Lagrange_Mesh_Base * aMesh ) : mParameters( aMesh->get_parameters() ),
                                             mMesh( aMesh ),
                                             mNumberOfDimensions( aMesh->get_parameters()->get_number_of_dimensions() )
    {
    }

    // ----------------------------------------------------------------------------
    void STK::create_mesh_data( const double aTimeStep )
    {
        // start timer
        tic tTimer;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: get number of nodes, elements etc
        // - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // activate this pattern on background mesh
        mMesh->select_activation_pattern();

        // get number of elements
        uint tNumberOfElements        = mMesh->get_number_of_elements();

        // get number of nodes
        uint tNumberOfNodes           = mMesh->get_number_of_nodes_on_proc();

        // get number of nodes per element
        uint tNumberOfNodesPerElement = mMesh->get_number_of_bases_per_element();

        // - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: Allocate Matrices
        // - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // first field is element level
        Matrix< DDRMat> & tElementLevels = mMesh->get_real_scalar_field_data( 0 );
        tElementLevels.set_size( tNumberOfElements, 1 );

        // second field is element owner
        Matrix< DDRMat> & tElementOwners= mMesh->get_real_scalar_field_data( 1 );
        tElementOwners.set_size( tNumberOfElements, 1 );

        // third field is node Ids
        Matrix< DDRMat> & tNodeIDs = mMesh->get_real_scalar_field_data( 2 );
        tNodeIDs.set_size( tNumberOfNodes, 1 );

        // forth field is node Inds
    //        Matrix< DDRMat> & tElementIndices = mMesh->get_real_scalar_field_data( 3 );
    //        tElementIndices.set_size( tNumberOfElements, 1 );

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

        // initialize node ownership
        mNodeSharing = Matrix<IdMat> ( tNumberOfNodes, mMesh->get_background_mesh()->get_number_of_proc_neighbors() );
        mNodeSharing.fill(MORIS_ID_MAX);

        // - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: Populate Matrices
        // - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // loop over all elements
        for( uint e=0; e<tNumberOfElements; ++e )
        {
            // get pointer to element
            Element * tElement = mMesh->get_element( e );

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

            // save index of element
    //            tElementIndices( e ) = tElement->get_index();
        }

        // loop over all nodes
        for( uint k=0; k<tNumberOfNodes; ++k )
        {
            auto tNode = mMesh->get_node_by_index( k );

            Matrix< DDRMat > tNodeCoords = tNode->get_coords();

            // copy coords to output matrix
            mNodeCoords.set_row( k , tNodeCoords );

            // copy node Owner
            mNodeOwner( k ) = tNode->get_owner();

            // Add node sharing
           if(tNode->has_node_sharing())
           {
               // Get node sharing
               Matrix<IdMat> tNodeSharing = tNode->get_node_sharing( );
               // Iterate through sharing
               for(uint iShare = 0; iShare<tNodeSharing.numel(); iShare++)
               {
                   mNodeSharing(k,iShare) = tNodeSharing(iShare);
               }
           }

            // copy node index into map
            mNodeLocalToGlobal( k ) = tNode->get_id();

            // save vertex id
            tNodeIDs( k ) = tNode->get_id();
        }

        // special function for old mesh
        this->flag_old_and_new_elements();

        // - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Link MTK Object
        // - - - - - - - - - - - - - - - - - - - - - - - - - - -

        uint tNumberOfRealScalarFields = mMesh->get_number_of_real_scalar_fields();

        if( mMesh->get_activation_pattern() == mParameters->get_lagrange_input_pattern() )
        {
            // Initialize scalar field data
            mRealScalarFields = moris::Vector< mtk::Scalar_Field_Info<DDRMat> >( tNumberOfRealScalarFields+1 );
        }
        else
        {
            // Initialize scalar field data
            mRealScalarFields = moris::Vector< mtk::Scalar_Field_Info<DDRMat> >( tNumberOfRealScalarFields );
        }

        // add real scalar fields
        for( uint f=0; f<tNumberOfRealScalarFields; ++f )
        {
            mRealScalarFields( f ).set_field_name(        mMesh->get_real_scalar_field_label( f ) );
            mRealScalarFields( f ).set_field_entity_rank( mMesh->get_real_scalar_field_rank ( f ) );

            switch( mMesh->get_real_scalar_field_rank( f ) )
            {
                case( mtk::EntityRank::NODE ) :
                {
                    mRealScalarFields( f ).add_field_data( & mNodeLocalToGlobal,
                                                           & mMesh->get_real_scalar_field_data( f ) );
                    break;
                }
                case( mtk::EntityRank::ELEMENT ) :
                {
                    mRealScalarFields( f ).add_field_data( & mElementLocalToGlobal,
                                                           & mMesh->get_real_scalar_field_data( f ) );
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "Invalid entity of real scalar field" );
                    break;
                }
            }

            mFieldsInfo.mRealScalarFields.push_back( &mRealScalarFields( f ) );
        }

        // link mesh data object
        mMeshData.SpatialDim                = & mNumberOfDimensions;
        mMeshData.ElemConn( 0 )             = & mElementTopology;
        mMeshData.NodeCoords                = & mNodeCoords;
        mMeshData.NodeProcsShared           = & mNodeSharing;
        mMeshData.LocaltoGlobalElemMap( 0 ) = & mElementLocalToGlobal;
        mMeshData.LocaltoGlobalNodeMap      = & mNodeLocalToGlobal;
        mMeshData.FieldsInfo                = & mFieldsInfo;
        mMeshData.SetsInfo                  = & mSetsInfo;

        // set timestep of mesh data object
        mMeshData.TimeStamp = aTimeStep;
        mMeshData.AutoAuraOptionInSTK = false;
        mMeshData.CreateAllEdgesAndFaces = true;

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

        // stop timer
        real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

        // print output
    //        MORIS_LOG_INFO( "%s Created MTK output object.\n               Mesh has %lu elements and %lu nodes.\n               Creation took %5.3f seconds.\n\n",
    //                proc_string().c_str(),
    //                ( long unsigned int ) tNumberOfElements,
    //                ( long unsigned int ) tNumberOfNodes,
    //                ( double ) tElapsedTime / 1000);
        MORIS_LOG_INFO( "%s Created MTK output object.",
                proc_string().c_str());
        MORIS_LOG_INFO( "Mesh has %lu elements and %lu nodes.",
                ( long unsigned int ) tNumberOfElements,
                ( long unsigned int ) tNumberOfNodes);
        MORIS_LOG_INFO(  "Creation took %5.3f seconds.",
                ( double ) tElapsedTime / 1000);
        MORIS_LOG_INFO( " " );

    }

    // ----------------------------------------------------------------------------

    void STK::save_to_file( const std::string & aFilePath )
    {
        tic tTimer;

    //        mMeshData.print_details();

        // create database object
        moris::mtk::Mesh_Core_STK tMesh( mMeshData );

        // copy file path, since tMesh does not like const input
        std::string tFilePath = aFilePath;

        // flag to add element cmap to exodus file
        bool tAddElemCmap = true;

        // save file
        tMesh.create_output_mesh( tFilePath, tAddElemCmap );

        // stop timer
        real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

        // print output
    //        MORIS_LOG_INFO( "%s Wrote MTK mesh to file %s.\n               Writing took %5.3f seconds.\n\n",
    //                proc_string().c_str(),
    //                aFilePath.c_str(),
    //                ( double ) tElapsedTime / 1000 );
        MORIS_LOG_INFO( "%s Wrote MTK mesh to file %s.",
                proc_string().c_str(),
                aFilePath.c_str());
        MORIS_LOG_INFO( "Writing took %5.3f seconds.",
                ( double ) tElapsedTime / 1000 );
        MORIS_LOG_INFO( " " );
    }

    // ----------------------------------------------------------------------------

    void STK::flag_old_and_new_elements()
    {
        uint tInputPattern  = mParameters->get_lagrange_input_pattern();
        uint tOutputPattern = mParameters->get_lagrange_output_pattern();

        if( mMesh->get_activation_pattern() ==  tInputPattern )
        {
            // get number of elements on mesh
            uint tNumberOfElements = mMesh->get_number_of_elements();

            uint tFieldIndex = mMesh->create_real_scalar_field_data( "Refinement", mtk::EntityRank::ELEMENT );

            // link to field
            Matrix< DDRMat > & tData = mMesh->get_real_scalar_field_data( tFieldIndex );

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
// ----------------------------------------------------------------------------
} /* namespace moris */

