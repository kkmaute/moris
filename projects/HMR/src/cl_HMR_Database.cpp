/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Database.cpp
 *
 */

#include "cl_HMR_Database.hpp"

#include "cl_HMR_Field.hpp"
#include "cl_HMR_File.hpp"
#include "cl_HMR_Mesh.hpp"
#include "MTK_Tools.hpp"
#include "cl_Tracer.hpp"

#include "op_times.hpp"
#include "fn_dot.hpp"
#include "fn_eye.hpp"
#include "fn_unique.hpp"

namespace moris::hmr
{
    // -----------------------------------------------------------------------------

    Database::Database( Parameters* aParameters )
            : mParameters( aParameters )
    {
        // Log Trace this function
        Tracer tTracer( "HMR", "Database", "Create" );

        // lock parameters ( number of elements per direction etc )
        aParameters->lock();

        // create factory object
        Factory tFactory( mParameters );

        // create background mesh object and initialize it with the coarsest level
        mBackgroundMesh = tFactory.create_background_mesh();

        // update element table
        mBackgroundMesh->collect_active_elements();

        // reset other patterns
        for ( uint iPattern = 0; iPattern < gNumberOfPatterns; ++iPattern )
        {
            mBackgroundMesh->reset_pattern( iPattern );
        }

        mBackgroundMesh->set_activation_pattern( 0 );

        // initialize mesh objects
        this->create_meshes();

        mAdditionalLagrangeMeshes.resize( gNumberOfPatterns );
    }

    // -----------------------------------------------------------------------------

    Database::Database( const std::string& aPath )
    {
        // create file object
        File tHDF5;

        // open file on disk
        tHDF5.open( aPath );

        // create new parameter pointer
        mParameters = new Parameters;

        // load settings
        tHDF5.load_settings( mParameters );

        // close file
        tHDF5.close();

        // create factory
        Factory tFactory( mParameters );

        // create background mesh object
        mBackgroundMesh = tFactory.create_background_mesh();

        // reset all patterns
        for ( uint k = 0; k < gNumberOfPatterns; ++k )
        {
            mBackgroundMesh->reset_pattern( k );
        }

        this->load_pattern_from_hdf5_file( aPath );

        // initialize mesh objects
        this->create_meshes();

        mAdditionalLagrangeMeshes.resize( gNumberOfPatterns );

        // activate input pattern
        this->set_activation_pattern( Parameters::mLagrangeInputPattern );
    }

    // -----------------------------------------------------------------------------

    Database::~Database()
    {
        // delete B-Spline and Lagrange meshes
        this->delete_meshes();

        // delete Background Mesh
        delete mBackgroundMesh;
        mBackgroundMesh = nullptr;

        // delete parameters
        if ( mDeleteParametersOnDestruction )
        {
            delete mParameters;
            mParameters = nullptr;
        }
    }

    // -----------------------------------------------------------------------------
    void
    Database::delete_additional_meshes( uint aPattern )
    {
        for ( auto tMesh : mAdditionalLagrangeMeshes( aPattern ) )
        {
            // delete this mesh
            delete tMesh;
        }
        mAdditionalLagrangeMeshes( aPattern ).clear();
    }

    // -----------------------------------------------------------------------------

    void
    Database::load_pattern_from_hdf5_file(
            const std::string& aPath )
    {
        // create file object
        File tHDF5;

        // open file on disk
        tHDF5.open( aPath );

        // load input pattern into file
        tHDF5.load_refinement_pattern( mBackgroundMesh );

        // close hdf5 file
        tHDF5.close();
    }

    // -----------------------------------------------------------------------------

    void
    Database::load_refinement_pattern(
            Matrix< DDLUMat >&                aElementCounterPerLevelAndPattern,
            Vector< Matrix< DDLUMat > >& aElementPerPattern,
            Matrix< DDUMat >&                 aPatternListUniqueMat )
    {
        uint tNumPattern = aElementPerPattern.size();

        // get number of levels
        uint tNumberOfLevels = aElementCounterPerLevelAndPattern.n_rows();

        for ( uint Ik = 0; Ik < tNumPattern; Ik++ )
        {
            // reset counter
            luint tCount = 0;

            // select B-Spline pattern
            mBackgroundMesh->set_activation_pattern( aPatternListUniqueMat( Ik ) );

            // loop over all levels
            for ( uint l = 0; l < tNumberOfLevels; ++l )
            {
                // cell which contains elements
                Vector< Background_Element_Base* > tElements;

                // collect elements from this level
                mBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tElements );

                // create a map with ids
                map< moris_id, luint > tMap;

                luint j = 0;
                for ( Background_Element_Base* tElement : tElements )
                {
                    tMap[ tElement->get_hmr_id() ] = j++;
                }

                luint tNumberOfElements = aElementCounterPerLevelAndPattern( l, Ik );

                for ( luint k = 0; k < tNumberOfElements; ++k )
                {
                    tElements( tMap.find( aElementPerPattern( Ik )( tCount++ ) ) )->put_on_refinement_queue();
                }

                // refine mesh
                mBackgroundMesh->perform_refinement( aPatternListUniqueMat( Ik ) );
            }
        }

        mBackgroundMesh->update_database();
    }

    // -----------------------------------------------------------------------------

    /**
     * sets the flag that the parameter object must be deleted
     * by the destructor
     */
    void
    Database::set_parameter_owning_flag()
    {
        mDeleteParametersOnDestruction = true;
    }

    // -----------------------------------------------------------------------------

    /**
     * sets the flag that the parameter object must not be deleted
     * by the destructor
     */
    void
    Database::unset_parameter_owning_flag()
    {
        mDeleteParametersOnDestruction = false;
    }

    // -----------------------------------------------------------------------------

    //        void Database::create_meshes()
    //        {
    //            // delete existing meshes
    //            this->delete_meshes();
    //
    //            // create factory object
    //            Factory tFactory;
    //
    //            // create BSpline meshes
    //            uint tNumberOfBSplineMeshes = mParameters->get_number_of_bspline_meshes();
    //
    //            // assign memory for B-Spline meshes
    //            mBSplineMeshes.resize ( tNumberOfBSplineMeshes, nullptr );
    //
    //            for( uint k=0; k<tNumberOfBSplineMeshes; ++k )
    //            {
    //                mBSplineMeshes( k ) = tFactory.create_bspline_mesh(
    //                                                                    mBackgroundMesh,
    //                                                                    mParameters->get_bspline_pattern( k ),
    //                                                                    mParameters->get_bspline_order( k ) );
    //
    //                mBSplineMeshes( k )->set_index( k );
    //            }
    //
    //            // create Lagrange meshes
    //            uint tNumberOfLagrangeMeshes = mParameters->get_number_of_lagrange_meshes();
    //
    //            // assign memory for Lagrange meshes
    //            mLagrangeMeshes.resize ( tNumberOfLagrangeMeshes, nullptr );
    //
    //            for( uint k=0; k<tNumberOfLagrangeMeshes; ++k )
    //            {
    //                mLagrangeMeshes( k ) = tFactory.create_lagrange_mesh(
    //                                                                      mBackgroundMesh,
    //                                                                      mBSplineMeshes,
    //                                                                      mParameters->get_lagrange_pattern( k ),
    //                                                                      mParameters->get_lagrange_order( k ) );
    //
    //                mLagrangeMeshes( k )->set_index( k );
    //
    //                // link to sideset if this is an output mesh
    //                if ( mLagrangeMeshes( k )->get_activation_pattern() == mParameters->get_lagrange_output_pattern() )
    //                {
    //                    mLagrangeMeshes( k )->set_create_side_sets( true );
    //                }
    //            }
    //        }

    // -----------------------------------------------------------------------------

    void
    Database::create_meshes()
    {
        // log & trace this function
        Tracer tTracer( "HMR", "Database", "Create Meshes" );

        // delete existing meshes
        this->delete_meshes();

        // create factory object
        Factory tFactory( mParameters );

        // create B-Spline meshes
        uint tNumberOfBSplineMeshes = mParameters->get_number_of_bspline_meshes();

        // assign memory for B-Spline meshes
        mBSplineMeshes.resize( tNumberOfBSplineMeshes, nullptr );

        // create all B-Spline meshes requested
        for ( uint iBspMesh = 0; iBspMesh < tNumberOfBSplineMeshes; ++iBspMesh )
        {
            // get the mesh pattern and polynomial order associated with the current B-spline mesh
            uint tPatternForBspMesh   = mParameters->get_bspline_pattern( iBspMesh );
            uint tPolyOrderForBspMesh = mParameters->get_bspline_order( iBspMesh );

            // create and initialize the HMR B-spline meshes
            mBSplineMeshes( iBspMesh ) = tFactory.create_bspline_mesh(
                    mBackgroundMesh,
                    tPatternForBspMesh,
                    tPolyOrderForBspMesh,
                    iBspMesh );
        }

        // get number of Lagrange meshes to be created
        uint tNumberOfLagrangeMeshes = mParameters->get_number_of_lagrange_meshes();

        // assign memory for Lagrange meshes
        mLagrangeMeshes.resize( tNumberOfLagrangeMeshes, nullptr );

        // create all Lagrange meshes requested
        for ( uint iLagMesh = 0; iLagMesh < tNumberOfLagrangeMeshes; ++iLagMesh )
        {
            // get B-Spline index corresponding to current Lagrange mesh
            Vector< uint > tBsplineMeshIndices = mParameters->get_lagrange_to_bspline_mesh( iLagMesh );

            // create a cell containing B-Spline meshes associated with current Lagrange mesh
            Vector< BSpline_Mesh_Base* > tBsplineMeshes( tBsplineMeshIndices.size() );

            // pick out B-Spline meshes associated with current Lagrange mesh from HMR-global
            // list of B-Spline meshes and fill container with them
            for ( uint iBspMeshOnLagMesh = 0; iBspMeshOnLagMesh < tBsplineMeshIndices.size(); ++iBspMeshOnLagMesh )
            {
                // assign existing b-spline mesh to list of b-spline meshes
                if ( tBsplineMeshIndices( iBspMeshOnLagMesh ) != UINT_MAX )
                {
                    tBsplineMeshes( iBspMeshOnLagMesh ) = mBSplineMeshes( tBsplineMeshIndices( iBspMeshOnLagMesh ) );
                }
                else
                {
                    tBsplineMeshes( iBspMeshOnLagMesh ) = nullptr;
                }
            }

            // get information specifying what Lagrange mesh to create
            uint tPatternForLagMesh = mParameters->get_lagrange_pattern( iLagMesh );
            uint tOrderForLagMesh   = mParameters->get_lagrange_order( iLagMesh );

            // create Lagrange mesh with links to all B-spline meshes associated to them
            mLagrangeMeshes( iLagMesh ) = tFactory.create_lagrange_mesh(
                    mBackgroundMesh,
                    tBsplineMeshes,
                    tPatternForLagMesh,
                    tOrderForLagMesh,
                    iLagMesh );

            // link to side-set if this is an output mesh
            if ( mParameters->is_output_mesh( iLagMesh ) )
            {
                // FIXME make mOutputSideSets more general. Set with mesh
                mLagrangeMeshes( iLagMesh )->set_side_sets( mOutputSideSets );
            }
        }    // end for: each Lagrange mesh to be created
    }        // end function: hmr::Database::create_meshes()

    // -----------------------------------------------------------------------------

    void
    Database::delete_meshes()
    {
        // delete all pointers
        for ( auto tMesh : mBSplineMeshes )
        {
            // delete this mesh
            delete tMesh;
        }

        mBSplineMeshes.clear();

        // delete all pointers
        for ( auto tMesh : mLagrangeMeshes )
        {
            // delete this mesh
            delete tMesh;
        }

        mLagrangeMeshes.clear();

        // delete all pointers
        for ( auto tMeshes : mAdditionalLagrangeMeshes )
        {
            for ( auto tMesh : tMeshes )
            {
                // delete this mesh
                delete tMesh;
            }
            tMeshes.clear();
        }

        mAdditionalLagrangeMeshes.clear();
    }

    // -----------------------------------------------------------------------------

    void
    Database::update_bspline_meshes()
    {
        // remember active pattern // uint
        auto tActivePattern = mBackgroundMesh->get_activation_pattern();

        // update all B-Spline meshes
        for ( auto tMesh : mBSplineMeshes )
        {
            // synchronize mesh with background mesh
            tMesh->update_mesh();
        }

        // reset pattern
        mBackgroundMesh->set_activation_pattern( tActivePattern );
    }

    // -----------------------------------------------------------------------------

    void
    Database::update_bspline_meshes( uint aPattern )
    {
        // remember active pattern // uint
        auto tActivePattern = mBackgroundMesh->get_activation_pattern();

        // update all B-Spline meshes
        for ( auto tMesh : mBSplineMeshes )
        {
            if ( tMesh->get_activation_pattern() == aPattern )
            {
                // synchronize mesh with background mesh
                tMesh->update_mesh();
            }
        }

        // reset pattern
        mBackgroundMesh->set_activation_pattern( tActivePattern );
    }

    // -----------------------------------------------------------------------------

    void
    Database::update_lagrange_meshes()
    {
        // remember active pattern // uint
        auto tActivePattern = mBackgroundMesh->get_activation_pattern();

        // update all Lagrange meshes and link elements to their
        // B-Spline twins
        for ( auto tMesh : mLagrangeMeshes )
        {
            // synchronize mesh with background mesh
            tMesh->update_mesh();
        }

        // reset pattern
        mBackgroundMesh->set_activation_pattern( tActivePattern );
    }

    // -----------------------------------------------------------------------------

    void
    Database::update_lagrange_meshes( uint aPattern )
    {
        // remember active pattern // uint
        auto tActivePattern = mBackgroundMesh->get_activation_pattern();

        // update all Lagrange meshes and link elements to their
        // B-Spline twins
        for ( auto tMesh : mLagrangeMeshes )
        {
            if ( tMesh->get_activation_pattern() == aPattern )
            {
                // synchronize mesh with background mesh
                tMesh->update_mesh();
            }
        }

        // reset pattern
        mBackgroundMesh->set_activation_pattern( tActivePattern );
    }

    // -----------------------------------------------------------------------------

    void
    Database::finalize()
    {
        // log & trace this operation
        Tracer tTracer( "HMR", "Database", "Finalize" );

        // make sure this function is only called once
        MORIS_ERROR(
                !mFinalizedCalled,
                "hmr::Database::finalize() - Finalize was called earlier. You should only call it once." );

        // remember active pattern
        auto tActivePattern = mBackgroundMesh->get_activation_pattern();

        const Vector< Vector< uint > >& tOutputMeshIndices = mParameters->get_output_mesh();

        // finalize each of the output meshes
        if ( tOutputMeshIndices.size() > 0 )
        {
            for ( uint iOutputMesh = 0; iOutputMesh < tOutputMeshIndices( 0 ).size(); iOutputMesh++ )
            {
                // get the corresponding mesh index
                uint tMeshIndex = tOutputMeshIndices( 0 )( iOutputMesh );

                // activate output pattern
                uint tOutputMeshPattern = mLagrangeMeshes( tMeshIndex )->get_activation_pattern();
                mBackgroundMesh->set_activation_pattern( tOutputMeshPattern );

                // create communication table
                this->create_communication_table();

                if ( mParameters->get_number_of_dimensions() == 3 )
                {
                    mBackgroundMesh->create_faces_and_edges();
                }
                else
                {
                    mBackgroundMesh->create_facets();
                }

                if ( not mParameters->get_background_mesh_file_name().empty() )
                {
                    mBackgroundMesh->save_to_vtk( mParameters->get_background_mesh_file_name() );
                }
            }
        }

        // finalize each of the Lagrange meshes
        for ( Lagrange_Mesh_Base* tMesh : mLagrangeMeshes )
        {
            // fixme: check effect of this flag
            // if ( ( ! mHaveInputTMatrix || mParameters->get_lagrange_input_pattern() != tMesh->get_activation_pattern() ) )
            //{
            tMesh->calculate_node_indices();
            tMesh->calculate_node_sharing();
            tMesh->calculate_t_matrices();
            //}

            // only needed for output mesh
            if ( mParameters->is_output_mesh( tMesh->get_index() ) )
            {
                // create facets
                tMesh->create_facets();

                // With the constructed facet information create facet clustering data
                tMesh->create_facet_clusters();

                // create edges
                if ( mParameters->get_number_of_dimensions() == 3 )
                {
                    tMesh->create_edges();
                }

                if ( not mParameters->get_lagrange_mesh_file_name().empty() )
                {
                    tMesh->save_to_vtk( mParameters->get_lagrange_mesh_file_name() );
                }
            }
        }

        // set flag for input t-matrices
        mHaveInputTMatrix = true;

        // finalize each of the B-spline meshes
        for ( auto tMesh : mBSplineMeshes )
        {
            tMesh->calculate_basis_indices( mCommunicationTable );
#ifdef MORIS_HAVE_DEBUG
            tMesh->calculate_basis_coordinates();

            const std::string tBasisFunctionVtkFileName = mParameters->get_basis_fuction_vtk_file_name();

            if ( tBasisFunctionVtkFileName.size() > 0 )
            {
                MORIS_LOG_INFO( "Save basis function to vtk file: %s", tBasisFunctionVtkFileName.c_str() );

                tMesh->save_to_vtk( tBasisFunctionVtkFileName );
            }
#endif
        }

        // reset active pattern
        if ( mBackgroundMesh->get_activation_pattern() != tActivePattern )
        {
            mBackgroundMesh->set_activation_pattern( tActivePattern );
        }

#ifdef MORIS_HAVE_DEBUG
        this->check_entity_ids();
#endif

        // create sidesets for output pattern
        this->create_side_sets();

        mFinalizedCalled = true;
    }

    // -----------------------------------------------------------------------------

    void
    Database::reset_refined_meshes()
    {
        // remember active pattern
        auto tActivePattern = mBackgroundMesh->get_activation_pattern();

        this->delete_side_sets();

        //            if( mParameters->get_number_of_dimensions() == 3 )
        //            {
        //                mBackgroundMesh->delete_faces();
        //                mBackgroundMesh->delete_edges();
        //            }
        //            else
        //            {
        //                mBackgroundMesh->delete_faces();
        //            }
        //
        //            mBackgroundMesh->reset_neighbors();

        // reset other patterns
        for ( uint k = 0; k < gNumberOfPatterns; ++k )
        {
            mBackgroundMesh->reset_pattern( k );
        }

        // update neighborhood tables
        mBackgroundMesh->update_database();

        // reset active pattern
        if ( mBackgroundMesh->get_activation_pattern() != tActivePattern )
        {
            mBackgroundMesh->set_activation_pattern( tActivePattern );
        }

        this->create_meshes();

        // set flag for input t-matrices
        mHaveInputTMatrix = false;

        mFinalizedCalled = false;
    }

    // -----------------------------------------------------------------------------

    void
    Database::create_communication_table()
    {
        moris_id tParSize = par_size();
        moris_id tMyRank  = par_rank();

        if ( tParSize > 1 )
        {
            // in a first step, we identify all processors this proc wants
            // to talk to

            // this is a Boolean-like matrix
            Matrix< IdMat > tColumn( tParSize, 1, 0 );

            // test owners of B-Splines
            for ( auto tMesh : mBSplineMeshes )
            {
                // get number of active B-Splines
                auto tNumberOfBSplines = tMesh->get_number_of_active_basis_on_proc();

                // loop over all active basis on this mesh
                for ( uint k = 0; k < tNumberOfBSplines; ++k )
                {
                    // set flag for this proc
                    tColumn( tMesh->get_active_basis( k )->get_owner() ) = 1;
                }
            }

            // remove self from row
            tColumn( tMyRank ) = 0;

            // communication table
            Matrix< IdMat > tCommTable;

            // matrices to send
            Vector< Matrix< IdMat > > tSend;

            // matrices to receive
            Vector< Matrix< IdMat > > tRecv;

            if ( tMyRank != 0 )
            {
                // create communication table with one entry
                tCommTable.set_size( 1, 1, 0 );
                tSend.resize( 1, tColumn );
            }
            else
            {
                // create comm matrix
                tCommTable.set_size( tParSize, 1, 0 );

                // communicate with all other procs
                for ( moris_id k = 1; k < tParSize; ++k )
                {
                    tCommTable( k ) = k;
                }

                // nothing to send
                Matrix< IdMat > tEmpty;
                tSend.resize( tParSize, tEmpty );
            }

            // exchange matrices
            communicate_mats( tCommTable, tSend, tRecv );

            // process information on leader proc
            if ( tMyRank == 0 )
            {
                // create communication matrix
                Matrix< IdMat > tCommMatrix( tParSize, tParSize, 0 );

                // process first row
                tRecv( 0 ) = tColumn;

                // loop over all procs and create comm matrix
                for ( moris_id j = 0; j < tParSize; ++j )
                {
                    for ( moris_id i = 0; i < tParSize; ++i )
                    {
                        if ( tRecv( j )( i, 0 ) != 0 )
                        {
                            tCommMatrix( i, j ) = 1;
                            tCommMatrix( j, i ) = 1;
                        }
                    }
                }

                // remove diagonal
                for ( moris_id i = 0; i < tParSize; ++i )
                {
                    tCommMatrix( i, i ) = 0;
                }

                // create sending list
                Matrix< IdMat > tEmpty;
                tSend.resize( tParSize, tEmpty );

                for ( moris_id j = 0; j < tParSize; ++j )
                {
                    // count nonzero entries
                    uint tCount = 0;
                    for ( moris_id i = 0; i < tParSize; ++i )
                    {
                        if ( tCommMatrix( i, j ) != 0 )
                        {
                            ++tCount;
                        }
                    }

                    // assign memory
                    tSend( j ).set_size( tCount, 1, 0 );

                    // reset counter
                    tCount = 0;

                    // write values into matrix
                    for ( moris_id i = 0; i < tParSize; ++i )
                    {
                        if ( tCommMatrix( i, j ) != 0 )
                        {
                            tSend( j )( tCount++ ) = i;
                        }
                    }
                }
            }

            // exchange matrices with other procs
            communicate_mats( tCommTable, tSend, tRecv );

            // write data into communication table
            if ( tMyRank == 0 )
            {
                mCommunicationTable = tSend( 0 );
            }
            else
            {
                mCommunicationTable = tRecv( 0 );
            }
        }
        else    // if run in serial
        {
            // communication table is empty
            mCommunicationTable.set_size( 0, 1 );
        }
    }

    // -----------------------------------------------------------------------------

    Background_Mesh_Base*
    Database::get_background_mesh()
    {
        return mBackgroundMesh;
    }

    // -----------------------------------------------------------------------------
    /**
     * creates a union of two patterns
     */
    void
    Database::unite_patterns(
            uint aSourceA,
            uint aSourceB,
            uint aTarget )
    {
        tic tTimer;

        this->delete_additional_meshes( aTarget );

        mBackgroundMesh->unite_patterns( aSourceA,
                aSourceB,
                aTarget );

        // stop timer
        real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;

        MORIS_LOG_INFO( "%s United patterns %lu and %lu to %lu.",
                proc_string().c_str(),
                (long unsigned int)aSourceA,
                (long unsigned int)aSourceB,
                (long unsigned int)aTarget );

        MORIS_LOG_INFO( "Calculation took %5.3f seconds.",
                (double)tElapsedTime / 1000 );
        MORIS_LOG_INFO( " " );
    }

    // -----------------------------------------------------------------------------

    /**
     * copies a source pattern to a target pattern
     */
    void
    Database::copy_pattern( uint aSource,
            uint                 aTarget )
    {
        tic tTimer;

        mBackgroundMesh->copy_pattern( aSource,
                aTarget );

        // stop timer
        real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;

        // print output
        MORIS_LOG_INFO( "%s Copied pattern %lu to %lu.",
                proc_string().c_str(),
                (long unsigned int)aSource,
                (long unsigned int)aTarget );

        MORIS_LOG_INFO( "Calculation took %5.3f seconds.",
                (double)tElapsedTime / 1000 );
        MORIS_LOG_INFO( " " );
    }

    // -----------------------------------------------------------------------------

    //        void Database::add_extra_refinement_step_for_exodus()
    //        {
    //            // get refined pattern
    //            auto tPattern = mParameters->get_refined_output_pattern();
    //
    //            // create refined pattern
    //            mBackgroundMesh->copy_pattern( mParameters->get_lagrange_output_pattern(),
    //                                           tPattern );
    //
    //            // activate output pattern
    //            mBackgroundMesh->set_activation_pattern( tPattern );
    //
    //            // collect all active elements on background mesh
    //            mBackgroundMesh->collect_active_elements();
    //
    //            // get number of elements
    //            luint tNumberOfElements = mBackgroundMesh->get_number_of_active_elements_on_proc();
    //
    //            // flag all active elements
    //            for( luint e = 0; e < tNumberOfElements; ++e )
    //            {
    //                mBackgroundMesh->get_element( e )->put_on_refinement_queue();
    //            }
    //
    //            // perform refinement
    //            mBackgroundMesh->perform_refinement();
    //        }

    // -----------------------------------------------------------------------------

    void
    Database::perform_refinement(
            const uint aActivePattern,
            const bool aResetPattern )
    {
        // flag for output
        bool tFlag = mHaveRefinedAtLeastOneElement;

        // get pointer to working pattern
        uint tWorkingPattern = mParameters->get_working_pattern();

        uint tMinLevel = 0;

        if ( not aResetPattern )
        {
            // minimum refinement level. Is zero by default
            tMinLevel = mParameters->get_initial_refinement( aActivePattern );
        }

        // this function resets the working pattern
        if ( aResetPattern )
        {
            mBackgroundMesh->reset_pattern( aActivePattern );
        }

        // set activation pattern
        mBackgroundMesh->set_activation_pattern( aActivePattern );

        // get max level on this mesh
        uint tMaxLevel = mBackgroundMesh->get_max_level();

        // initial refinement
        for ( uint l = 0; l < tMinLevel; ++l )
        {
            // container for elements on this level
            Vector< Background_Element_Base* > tElementList;

            // collect all elements on this level ( without aura )
            mBackgroundMesh->collect_elements_on_level( l, tElementList );

            // loop over all elements and flag them for refinement
            for ( Background_Element_Base* tElement : tElementList )
            {
                // test if element is marked on working pattern
                if ( !tElement->is_refined( aActivePattern ) )
                {
                    tElement->put_on_refinement_queue();
                }
            }

            if ( mBackgroundMesh->collect_refinement_queue() )
            {
                mBackgroundMesh->perform_refinement( aActivePattern );
            }
        }

        // loop over all levels
        for ( uint l = tMinLevel; l <= tMaxLevel; ++l )
        {
            // container for elements on this level
            Vector< Background_Element_Base* > tElementList;

            // collect all elements on this level ( without aura )
            mBackgroundMesh->collect_elements_on_level( l, tElementList );

            // loop over all elements and flag them for refinement
            for ( Background_Element_Base* tElement : tElementList )
            {
                // test if element is marked on working pattern
                if ( tElement->is_refined( tWorkingPattern ) )
                {
                    // flag this element
                    tElement->put_on_refinement_queue();
                }
            }
            tFlag = tFlag || mBackgroundMesh->perform_refinement( aActivePattern );
        }

        // #MINREF check for minimum refinement criterion
        while ( mBackgroundMesh->collect_refinement_queue() )
        {
            mBackgroundMesh->perform_refinement( aActivePattern );
        }

        // tidy up working pattern
        mBackgroundMesh->reset_pattern( tWorkingPattern );

        // create new B-Spline Meshes
        this->update_bspline_meshes( aActivePattern );

        // create new Lagrange meshes
        this->update_lagrange_meshes( aActivePattern );

        // remember flag
        mHaveRefinedAtLeastOneElement = tFlag;
    }

    // -----------------------------------------------------------------------------

    // interpolate field values from source Lagrange to target Lagrange mesh
    void
    Database::interpolate_field(
            uint                            aSourcePattern,
            const std::shared_ptr< Field >& aSource,
            uint                            aTargetPattern,
            const std::shared_ptr< Field >& aTarget )
    {
        // make sure that mesh orders match
        MORIS_ERROR( aSource->get_interpolation_order() == aTarget->get_interpolation_order(),
                "Database::interpolate_field: Source and Target Field must have same interpolation order" );

        // make sure that both fields are scalar or of equal dimension
        MORIS_ERROR( aSource->get_number_of_dimensions() == aTarget->get_number_of_dimensions(),
                "Database::interpolate_field: Source and Target Field must have same dimension" );

        // get interpolation order
        uint tOrder = mtk::interpolation_order_to_uint( aSource->get_interpolation_order() );

        // pointer to mesh that is linked to input field
        Lagrange_Mesh_Base* tSourceMesh = nullptr;

        bool tMeshFound = false;
        // find pointer to input mesh
        for ( Lagrange_Mesh_Base* tMesh : mLagrangeMeshes )
        {
            if ( tMesh->get_order() == tOrder && tMesh->get_activation_pattern() == aSourcePattern )
            {
                tSourceMesh = tMesh;
                tMeshFound  = true;
                break;
            }
        }

        if ( !tMeshFound )
        {
            // find pointer to input mesh
            for ( Lagrange_Mesh_Base* tMesh : mAdditionalLagrangeMeshes( aSourcePattern ) )
            {
                if ( tMesh->get_order() == tOrder && tMesh->get_activation_pattern() == aSourcePattern )
                {
                    tSourceMesh = tMesh;
                    break;
                }
            }
        }

        // pointer to mesh that is linked to output field
        Lagrange_Mesh_Base* tTargetMesh = nullptr;

        tMeshFound = false;

        // find pointer to output mesh
        for ( Lagrange_Mesh_Base* tMesh : mLagrangeMeshes )
        {
            if ( tMesh->get_order() == tOrder && tMesh->get_activation_pattern() == aTargetPattern )
            {
                tTargetMesh = tMesh;
                tMeshFound  = true;
                break;
            }
        }

        if ( !tMeshFound )
        {
            // find pointer to output mesh
            for ( Lagrange_Mesh_Base* tMesh : mAdditionalLagrangeMeshes( aTargetPattern ) )
            {
                if ( tMesh->get_order() == tOrder && tMesh->get_activation_pattern() == aTargetPattern )
                {
                    tTargetMesh = tMesh;
                    break;
                }
            }
        }

        tTargetMesh->select_activation_pattern();

        // unflag nodes on target
        tTargetMesh->unflag_all_basis();

        // number of elements on target mesh
        auto tNumberOfElements = tTargetMesh->get_number_of_elements();

        // number of nodes per element
        auto tNumberOfNodesPerElement = tTargetMesh->get_number_of_bases_per_element();

        // create unity matrix
        Matrix< DDRMat > tEye;
        eye( tNumberOfNodesPerElement, tNumberOfNodesPerElement, tEye );

        // get values of source field
        const Matrix< DDRMat >& tSourceData = aSource->get_node_values();

        // get target data
        Matrix< DDRMat >& tTargetData = aTarget->get_node_values();

        // allocate value matrix
        tTargetData.set_size( tTargetMesh->get_number_of_all_basis_on_proc(), aTarget->get_number_of_dimensions() );

        // containers for source and target data
        Matrix< DDRMat > tElementSourceData( tNumberOfNodesPerElement, aSource->get_number_of_dimensions() );

        // Create T-matrix
        Factory        tFactory( mParameters );
        T_Matrix_Base* tTMatrix = tFactory.create_t_matrix( tTargetMesh );

        // loop over all elements
        for ( luint e = 0; e < tNumberOfElements; ++e )
        {
            // get pointer to target element
            auto tTargetElement = tTargetMesh->get_element( e );

            // get background element
            auto tBackgroundElement = tTargetElement->get_background_element();

            // initialize refinement Matrix
            Matrix< DDRMat > tR( tEye );

            while ( !tBackgroundElement->is_active( aSourcePattern ) )
            {
                // right multiply refinement matrix
                tR = tR * tTMatrix->get_refinement_matrix( tBackgroundElement->get_child_index() );

                // jump to parent
                tBackgroundElement = tBackgroundElement->get_parent();
            }

            // get pointer to source element
            auto tSourceElement = tSourceMesh->get_element_by_memory_index( tBackgroundElement->get_memory_index() );

            // fill source data vector
            for ( uint k = 0; k < tNumberOfNodesPerElement; ++k )
            {
                // get pointer to source node
                auto tNode  = tSourceElement->get_basis( k );
                auto tIndex = tNode->get_index();

                // copy data from source mesh
                tElementSourceData.set_row( k, tSourceData.get_row( tIndex ) );
            }

            // copy target data to target mesh
            for ( uint k = 0; k < tNumberOfNodesPerElement; ++k )
            {
                // get pointer to target node
                auto tNode = tTargetElement->get_basis( k );

                // test if data has already been written to target
                if ( !tNode->is_flagged() )
                {
                    // get node index
                    auto tIndex = tNode->get_index();

                    tTargetData.set_row( tIndex, tR.get_row( k ) * tElementSourceData );

                    // flag this node
                    tNode->flag();
                }
            }
        }

        delete ( tTMatrix );
    }

    // -----------------------------------------------------------------------------

    void
    Database::change_field_order(
            const std::shared_ptr< Field >& aSource,
            const std::shared_ptr< Field >& aTarget )
    {
        // pointer to in mesh
        Lagrange_Mesh_Base* tSourceMesh = aSource->get_mesh();

        // pointer to out mesh
        Lagrange_Mesh_Base* tTargetMesh = aTarget->get_mesh();

        // make sure that meshes are compatible
        MORIS_ASSERT( tSourceMesh->get_activation_pattern() == tTargetMesh->get_activation_pattern(),
                "incompatible meshes in change_field_order()" );

        this->set_activation_pattern( tSourceMesh->get_activation_pattern() );

        // unflag all nodes on out mesh
        tTargetMesh->unflag_all_basis();

        // source values
        Matrix< DDRMat >& tSourceValues = aSource->get_node_values();

        // target values
        Matrix< DDRMat >& tTargetValues = aTarget->get_node_values();

        // allocate output memory
        tTargetValues.set_size( tTargetMesh->get_number_of_nodes_on_proc(), 1 );

        // get number of elements
        uint tNumberOfElements = tSourceMesh->get_number_of_elements();

        // create T-matrix object
        Factory        tFactory( mParameters );
        T_Matrix_Base* tTMatrix = tFactory.create_t_matrix( tSourceMesh );

        // Get target mesh order and change order matrix
        uint             tTargetMeshOrder = tTargetMesh->get_order();
        Matrix< DDRMat > tT               = tTMatrix->get_change_order_matrix( tTargetMeshOrder );

        // delete t-Matrix object
        delete tTMatrix;

        uint tNumberOfNodesPerSourceElement = tSourceMesh->get_number_of_bases_per_element();
        uint tNumberOfNodesPerTargetElement = tTargetMesh->get_number_of_bases_per_element();

        Matrix< DDRMat > tLocalSourceValues( tNumberOfNodesPerSourceElement, 1 );
        Matrix< DDRMat > tN( 1, tNumberOfNodesPerSourceElement );

        // loop over all elements
        for ( uint e = 0; e < tNumberOfElements; ++e )
        {
            // get pointer to source element
            const Element* tSourceElement = tSourceMesh->get_element( e );

            for ( uint i = 0; i < tNumberOfNodesPerSourceElement; ++i )
            {
                tLocalSourceValues( i ) = tSourceValues( tSourceElement->get_basis( i )->get_index() );
            }

            // get pointer to target element
            Element* tTargetElement = tTargetMesh->get_element( e );

            // loop over all nodes on target
            for ( uint k = 0; k < tNumberOfNodesPerTargetElement; ++k )
            {
                // get basis
                Basis* tNode = tTargetElement->get_basis( k );

                if ( !tNode->is_flagged() )
                {
                    // copy row from T-Matrix
                    tT.get_row( k, tN );

                    // interpolate values
                    tTargetValues( tNode->get_index() ) = dot( tN, tLocalSourceValues );

                    // flag node
                    tNode->flag();
                }
            }
        }
    }

    // -----------------------------------------------------------------------------

    void
    Database::check_entity_ids()
    {
        if ( par_size() > 1 )
        {
            tic tTimer;

            // loop over all Lagrange meshes
            for ( Lagrange_Mesh_Base* tMesh : mLagrangeMeshes )
            {
                tMesh->select_activation_pattern();

                // check elements
                uint tNumberOfEntities = tMesh->get_number_of_elements();

                moris_id tMaxID = tMesh->get_max_element_id();

                for ( uint k = 0; k < tNumberOfEntities; ++k )
                {
                    moris_id tID = tMesh->get_element( k )->get_id();

                    MORIS_ERROR( 0 < tID && tID <= tMaxID, "Invalid Element ID" );
                }

                if ( mParameters->is_output_mesh( tMesh->get_index() ) )
                {
                    // check facets
                    tNumberOfEntities = tMesh->get_number_of_facets();

                    tMaxID = tMesh->get_max_facet_id();

                    for ( uint k = 0; k < tNumberOfEntities; ++k )
                    {
                        moris_id tID = tMesh->get_facet( k )->get_id();
                        MORIS_ERROR( 0 < tID && tID <= tMaxID, "Database::check_entity_ids(), Invalid Facet ID" );
                    }

                    // check edges
                    if ( mParameters->get_number_of_dimensions() == 3 )
                    {
                        tNumberOfEntities = tMesh->get_number_of_edges();
                        tMaxID            = tMesh->get_max_edge_id();

                        for ( uint k = 0; k < tNumberOfEntities; ++k )
                        {
                            moris_id tID = tMesh->get_edge( k )->get_id();

                            MORIS_ERROR( 0 < tID && tID <= tMaxID, "Database::check_entity_ids(), Invalid Edge ID" );
                        }
                    }
                }

                tNumberOfEntities = tMesh->get_number_of_nodes_on_proc();

                tMaxID = tMesh->get_max_node_id();

                for ( uint k = 0; k < tNumberOfEntities; ++k )
                {
                    moris_id tID = tMesh->get_node_by_index( k )->get_id();

                    MORIS_ERROR( 0 < tID && tID <= tMaxID,
                            "Database::check_entity_ids(), Invalid Node ID %-5i for node with index %-5i",
                            tID,
                            k );
                }
            }

            // loop over all B-Spline meshes
            for ( BSpline_Mesh_Base* tMesh : mBSplineMeshes )
            {
                tMesh->select_activation_pattern();

                // get number of splines
                uint tNumberOfEntities = tMesh->get_number_of_active_basis_on_proc();

                for ( uint k = 0; k < tNumberOfEntities; ++k )
                {
                    if ( tMesh->get_active_basis( k )->is_flagged() )
                    {
                        MORIS_ERROR( tMesh->get_active_basis( k )->get_hmr_index() < gNoEntityID, "Database::check_entity_ids(), Invalid B-Spline ID" );
                    }
                }
            }

            // stop timer
            real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;

            // print output
            MORIS_LOG_INFO( "%s passed entity ID test, testing took %5.3f seconds.",
                    proc_string().c_str(),
                    (double)tElapsedTime / 1000 );
        }
    }

    // -----------------------------------------------------------------------------

    Matrix< DDUMat >
    Database::create_output_pattern_list()
    {
        const Vector< Vector< uint > >& OutputMeshIndex = mParameters->get_output_mesh();

        Matrix< DDUMat > tPatternList( OutputMeshIndex( 0 ).size(), 1, MORIS_UINT_MAX );

        for ( uint Ik = 0; Ik < OutputMeshIndex( 0 ).size(); ++Ik )
        {
            tPatternList( Ik ) = mLagrangeMeshes( OutputMeshIndex( 0 )( Ik ) )->get_activation_pattern();
        }

        MORIS_ASSERT( tPatternList.max() != MORIS_UINT_MAX, "Database::create_side_sets(), Output mesh return false pattern index" );

        // get pointer to a Lagrange Mesh that uses this pattern
        // which one does not matter, since all elements with same pattern
        // have the same IDs
        Matrix< DDUMat > tUniquePatternList;
        unique( tPatternList, tUniquePatternList );

        MORIS_ERROR( tUniquePatternList.numel() == 1, "Database::create_side_sets(), Side sets are only for one pattern implemented yet" );

        return tUniquePatternList;
    }

    // -----------------------------------------------------------------------------

    void
    Database::create_side_sets()
    {
        // report on this operation
        MORIS_LOG_INFO( "Creating side sets in HMR Database" );

        if ( mParameters->get_create_side_sets() )
        {
            uint tNumberOfSets = 2 * mParameters->get_number_of_dimensions();
            Side_Set tEmpty;

            // allocate output side set
            mOutputSideSets.resize( tNumberOfSets, tEmpty );

            // get pattern number
            Matrix< DDUMat > tPatternList = this->create_output_pattern_list();

            Lagrange_Mesh_Base* tMesh = nullptr;

            for ( Lagrange_Mesh_Base* tLMesh : mLagrangeMeshes )    // FIXME
            {
                if ( tLMesh->get_activation_pattern() == tPatternList( 0 ) )
                {
                    tMesh = tLMesh;
                    break;
                }
            }
            mOutputSideSetMap.clear();

            // create side sets for output mesh
            for ( uint iSet = 0; iSet < tNumberOfSets; iSet++ )
            {
                // collect elements from background mesh
                Vector< Background_Element_Base* > tBackElements;

                mBackgroundMesh->collect_side_set_elements(
                        tPatternList( 0 ),
                        iSet + 1,
                        tBackElements );

                // get number of elements
                uint tNumberOfElements = tBackElements.size();

                // get ref to side set
                Side_Set& tSideSet = mOutputSideSets( iSet );

                // create name
                tSideSet.mInfo.mSideSetName = "SideSet_" + std::to_string( iSet + 1 );

                mOutputSideSetMap[ tSideSet.mInfo.mSideSetName ] = iSet;

                // allocate memory for ids
                tSideSet.mElemIdsAndSideOrds.set_size( tNumberOfElements, 2 );

                // allocate memory for indices
                tSideSet.mElemIndices.set_size( tNumberOfElements, 1 );

                // initialize counter
                luint tCount = 0;

                // loop over all Background Elements
                for ( Background_Element_Base* tBackElement : tBackElements )
                {
                    // get pointer to element on Lagrange Mesh
                    Element* tElement = tMesh->get_element_by_memory_index( tBackElement->get_memory_index() );

                    // write element ID
                    tSideSet.mElemIdsAndSideOrds( tCount, 0 ) = tElement->get_id();

                    // write sideset ordinal
                    tSideSet.mElemIdsAndSideOrds( tCount, 1 ) = iSet;

                    // write element index
                    tSideSet.mElemIndices( tCount++ ) = tElement->get_index();
                }
            }

            // link sets with Lagrange meshes
            for ( Lagrange_Mesh_Base* tLMesh : mLagrangeMeshes )
            {
                if ( tLMesh->get_activation_pattern() == tPatternList( 0 ) )
                {
                    tLMesh->set_side_sets( mOutputSideSets );
                }
            }
        }
    }

    // -----------------------------------------------------------------------------

    /**
     * creates the sidesets
     */
    void
    Database::delete_side_sets()
    {
        mOutputSideSets.clear();

        mOutputSideSetMap.clear();
    }

    // ----------------------------------------------------------------------------

    void
    Database::flag_element( luint aIndex )
    {
        // flag element implies that a manual refinement is performed
        // therefore, we set the flag
        mHaveRefinedAtLeastOneElement = true;

        // manually put this element on the queue FIXME this function switches the mRefinementQueueFlag to true.
        // However, in set_refined_flag the mRefinementQueueFlag flag is set to false;
        mBackgroundMesh->get_element( aIndex )->put_on_refinement_queue();

        // also remember this element on the working pattern
        mBackgroundMesh->get_element( aIndex )->set_refined_flag( mParameters->get_working_pattern() );
    }

    // ----------------------------------------------------------------------------

    void
    Database::flag_parent( luint aIndex )
    {
        // get pointer to this element
        Background_Element_Base* tElement = mBackgroundMesh->get_element( aIndex );

        // check level
        if ( tElement->get_level() > 0 )
        {
            mHaveRefinedAtLeastOneElement = true;

            // get parent
            Background_Element_Base* tParent = tElement->get_parent();

            // flag parent
            tParent->put_on_refinement_queue();

            // also remember this element on the working pattern
            tParent->set_refined_flag( mParameters->get_working_pattern() );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Database::create_extra_refinement_buffer_for_level( const uint aLevel )
    {
        // collect elements from level
        Vector< Background_Element_Base* > tElementsOfLevel;

        // uint tActivationPatter = mParameters->get_working_pattern();

        // collect elements on level
        mBackgroundMesh->collect_elements_on_level( aLevel, tElementsOfLevel );

        // count flagged elements
        luint tCount = 0;
        for ( Background_Element_Base* tElement : tElementsOfLevel )
        {
            if ( tElement->is_queued_for_refinement() )
            {
                ++tCount;
            }
        }

        // create container for elements
        Vector< Background_Element_Base* > tElements( tCount, nullptr );

        // reset counter
        tCount = 0;

        for ( Background_Element_Base* tElement : tElementsOfLevel )
        {
            if ( tElement->is_queued_for_refinement() )
            {
                tElements( tCount++ ) = tElement;
            }
        }

        // get buffer
        uint tBuffer = mParameters->get_refinement_buffer();

        if ( aLevel > 0 )
        {
            // loop over all elements and make sure that neighbors exist
            for ( Background_Element_Base* tElement : tElements )
            {
                // grab parent
                Background_Element_Base* tParent = tElement->get_parent();

                // container for neighbors
                Vector< Background_Element_Base* > tNeighbors;

                // get neighbors of parent ( because of the staircase buffer, they must exist )
                tParent->get_neighbors_from_same_level( tBuffer, tNeighbors );

                // loop over all neighbors
                for ( Background_Element_Base* tNeighbor : tNeighbors )
                {
                    // check if neighbor has children
                    if ( !tNeighbor->has_children() )
                    {
                        // remember refinement flag
                        bool tFlag = tNeighbor->is_queued_for_refinement();

                        // tell mesh to create children and keep state
                        mBackgroundMesh->refine_element( tNeighbor, true );

                        if ( tFlag )
                        {
                            tNeighbor->put_on_refinement_queue();
                        }
                    }
                }
            }
        }

        // update neighborhood tables
        mBackgroundMesh->update_database();

        // loop over all elements in queue
        for ( Background_Element_Base* tElement : tElements )
        {
            // get element neighbors
            // container for neighbors
            Vector< Background_Element_Base* > tNeighbors;

            // get neighbors of element
            tElement->get_neighbors_from_same_level( tBuffer, tNeighbors );

            for ( Background_Element_Base* tNeighbor : tNeighbors )
            {
                // flag neighbor
                tNeighbor->put_on_refinement_queue();

                // set flag on working pattern
                // tNeighbor->set_refined_flag( tActivationPatter );
            }
        }
    }

    // ----------------------------------------------------------------------------
}    // namespace moris::hmr
