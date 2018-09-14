/*
 * cl_HMR.cpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */
#include "op_times.hpp" //LNA/src
#include "fn_trans.hpp" //LNA/src
#include "cl_HMR.hpp" //HMR/src
#include "cl_HMR_Interface.hpp" //HMR/src
#include "cl_HMR_STK.hpp" //HMR/src
#include "cl_HMR_File.hpp" //HMR/src
#include "cl_HMR_Refinement_Manager.hpp"  //HMR/src
#include "cl_MDL_Model.hpp"
#include "cl_FEM_IWG_L2.hpp"

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        // default constuctor
        HMR::HMR ( Parameters * aParameters ) :
                mParameters( aParameters )
        {
            // locl parameters ( number of elements per direction etc )
            aParameters->lock();

            // create factory object
            Factory tFactory;

            // create background mesh object
            mBackgroundMesh = tFactory.create_background_mesh( mParameters );

            // update element table
            mBackgroundMesh->collect_active_elements();

            // reset other patterns
            for( uint k=0; k<gNumberOfPatterns; ++ k )
            {
                mBackgroundMesh->reset_pattern( k );
            }

            // initialize mesh objects
            this->create_meshes();

            // initialize T-Matrix objects
            this->init_t_matrices();
        }
// -----------------------------------------------------------------------------

        // alternative constuctor that converts ref to a pointer
        HMR::HMR ( Parameters & aParameters ) :
                                HMR( & aParameters )
        {

        }

// -----------------------------------------------------------------------------

        // alternative constuctor that uses parameter list
        HMR::HMR ( ParameterList & aParameterList )
            : HMR( new Parameters( aParameterList ) )
        {
            mDeleteParametersOnDestruction = true;
        }

// -----------------------------------------------------------------------------

        HMR::~HMR()
        {
            // delete T-Matrix objects
            this->delete_t_matrices();

            // delete B-Spline and Lagrange meshes
            this->delete_meshes();

            // delete Background Mesh
            delete mBackgroundMesh;

            /* for( auto tField: mFields )
            {
                delete tField;
            } */

            // delete parameters
            if ( mDeleteParametersOnDestruction )
            {
                delete mParameters;
            }
        }

// -----------------------------------------------------------------------------

        void
        HMR::delete_meshes()
        {

            // delete all pointers
            for( auto tMesh : mBSplineMeshes )
            {
                // delete this mesh
                delete tMesh;
            }

            // delete all pointers
            for( auto tMesh : mLagrangeMeshes )
            {
                // delete this mesh
                delete tMesh;
            }

        }
// -----------------------------------------------------------------------------

        void
        HMR::create_meshes()
        {
            // delete existing meshes
            this->delete_meshes();

            // create factory object
            Factory tFactory;

            // create BSpline meshes
            uint tNumberOfBSplineMeshes
                 = mParameters->get_number_of_bspline_meshes();

            // assign memory for B-Spline meshes
            mBSplineMeshes.resize ( tNumberOfBSplineMeshes, nullptr );

            for( uint k=0; k<tNumberOfBSplineMeshes; ++k )
            {
                mBSplineMeshes( k ) = tFactory.create_bspline_mesh(
                        mParameters,
                        mBackgroundMesh,
                        mParameters->get_bspline_pattern( k ),
                        mParameters->get_bspline_order( k ) );

                mBSplineMeshes( k )->set_index( k );
            }

            // create Lagrange meshes
            uint tNumberOfLagrangeMeshes
                = mParameters->get_number_of_lagrange_meshes();

            // assign memory for Lagrange meshes
            mLagrangeMeshes.resize ( tNumberOfLagrangeMeshes, nullptr );

            for( uint k=0; k<tNumberOfLagrangeMeshes; ++k )
            {
                        mLagrangeMeshes( k ) = tFactory.create_lagrange_mesh(
                        mParameters,
                        mBackgroundMesh,
                        mBSplineMeshes( mParameters->get_lagrange_to_bspline( k ) ),
                        mParameters->get_lagrange_pattern( k ),
                        mParameters->get_lagrange_order( k ) );

                        mLagrangeMeshes( k )->set_index( k );
            }
        }

// -----------------------------------------------------------------------------

        //@ fixme only activate meshes selectively to its patterns
        void
        HMR::update_meshes( const uint & aPattern )
        {

            // remember active pattern
            auto tActivePattern = mBackgroundMesh->get_activation_pattern();

            // activate output pattern
            mBackgroundMesh->set_activation_pattern( aPattern );

            // update all B-Spline meshes
            for( auto tMesh : mBSplineMeshes )
            {
                if( tMesh->get_activation_pattern() == aPattern )
                {
                    // synchronize mesh with background mesh
                    tMesh->update_mesh();
                }
            }

            // update all Lagrange meshes and link elements to their
            // B-Spline twins
            for( auto tMesh : mLagrangeMeshes )
            {
                if( tMesh->get_activation_pattern() == aPattern )
                {
                    // synchronize mesh with background mesh
                    tMesh->update_mesh();
                }
            }

            // reset pattern
            mBackgroundMesh->set_activation_pattern( tActivePattern );
        }

// -----------------------------------------------------------------------------

        Interface *
        HMR::create_mtk_interface()
        {
            //return this->create_mtk_interface(
            //        mParameters->get_output_pattern() );

            return new Interface( *this, mParameters->get_output_pattern() );
        }

// -----------------------------------------------------------------------------

        Interface *
        HMR::create_mtk_interface( const uint & aActivationPattern )
        {
            return new Interface( *this, aActivationPattern );
        }

// -----------------------------------------------------------------------------

        void
        HMR::init_t_matrices()
        {
            // get number of meshes
            uint tNumberOfMeshes
                = mParameters->get_number_of_lagrange_meshes();

            // allocate T-Matrix cell
            mTMatrix.resize( tNumberOfMeshes, nullptr );

            for( uint k=0; k<tNumberOfMeshes; ++k )
            {
                uint tBSplineMeshIndex
                    = mParameters->get_lagrange_to_bspline( k );

                // test if both meshes exist
                if ( mBSplineMeshes( tBSplineMeshIndex ) != NULL
                    && mLagrangeMeshes( k ) != NULL )
                {
                    // initialize T-Matrix object
                    mTMatrix( k ) = new T_Matrix( mParameters,
                            mBSplineMeshes( tBSplineMeshIndex ),
                            mLagrangeMeshes( k ) );
                }
            }
        }

// -----------------------------------------------------------------------------

        /**
         * Returns the pointer to a T-Matrix object.
         * Needed by Field constructor.
         */
        T_Matrix *
        HMR::get_t_matrix( const uint & aLagrangeMeshIndex )
        {
            return mTMatrix( aLagrangeMeshIndex );
        }

// -----------------------------------------------------------------------------

        void
        HMR::delete_t_matrices()
        {
            // delete pointers of calculation objects
            for( auto tTMatrix :  mTMatrix )
            {
                if ( tTMatrix != NULL )
                {
                    delete tTMatrix;
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        HMR::synchronize_t_matrix_flags()
        {
            // synchronize flags for T-Matrices on aura
            mBackgroundMesh->synchronize_t_matrix_flags();

            // get proc neighbors
            auto tMyProcNeighbors = mBackgroundMesh->get_proc_neigbors();

            uint tNumberOfNeighbors = tMyProcNeighbors.length();

            uint tMyRank = par_rank();

            for( auto tMesh: mLagrangeMeshes )
            {
                // ask mesh about number of basis per element
                auto tNumberOfBasisPerElement
                    = tMesh->get_number_of_basis_per_element();

                // ask mesh about pattern
                auto tLagrangePattern = tMesh->get_activation_pattern();
                auto tBSplinePattern  = tMesh->get_bspline_pattern();

                // loop over all procs
                for ( uint p=0; p<tNumberOfNeighbors; ++p )
                {
                    // only do this if there is a neighbor
                    if(        tMyProcNeighbors( p ) != gNoProcNeighbor
                            && tMyProcNeighbors( p ) != tMyRank )
                    {
                        Cell< Background_Element_Base* > tElements;

                        // get active elements from aura
                        mBackgroundMesh->collect_active_elements_from_aura(
                                p, 0, tElements );

                        // loop over all elements from aura
                        for( auto tElement : tElements )
                        {
                            // test if element is flagged
                            if( tElement->get_t_matrix_flag( tLagrangePattern ) )
                            {
                                // get pointer to B-Spline Element
                                auto tBElement = tMesh->get_element_by_memory_index(
                                        tElement->get_memory_index() );

                                // loop over all basis of B-Spline element
                                for( uint k=0; k<tNumberOfBasisPerElement; ++k )
                                {
                                    // get pointer to basis
                                    auto tBasis = tBElement->get_basis( k );

                                    // test if basis is owned by current proc
                                    if( tBasis->get_owner() == tMyRank )
                                    {
                                        // find first element that is owned
                                        auto tNumberOfElements = tBasis->get_element_counter();

                                        for( uint i=0; i<tNumberOfElements; ++i )
                                        {
                                            auto tOtherElement = tBasis->get_element( i );

                                            // test if this element is owned
                                            if ( tOtherElement->get_owner() == tMyRank )
                                            {
                                                // flag T-Matrix of this element
                                                tOtherElement->get_background_element()->set_t_matrix_flag( tBSplinePattern );

                                                // exit loop
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        HMR::activate_all_t_matrices()
        {
            for( auto tMesh : mLagrangeMeshes )
            {
                // activate pattern on background mesh
                tMesh->select_activation_pattern();

                auto tNumberOfElements = tMesh->get_number_of_elements();

                // loop over all elements
                for( luint e=0; e<tNumberOfElements; ++e )
                {
                    // get pointer to element
                    auto tLagrangeElement = tMesh->get_element( e );

                    // flag this element
                    tLagrangeElement->set_t_matrix_flag();
                }
            }

        }

// -----------------------------------------------------------------------------

        /**
         * creates a union of two patterns
         */
        void
        HMR::unite_patterns(
                const uint & aSourceA,
                const uint & aSourceB,
                const uint & aTarget )
        {
            tic tTimer;

            mBackgroundMesh->unite_patterns(
                    aSourceA,
                    aSourceB,
                    aTarget );

            for( uint k=0; k<gNumberOfPatterns; ++k )
            {
                this->update_meshes( k );
            }
            // this->update_meshes( aTarget );

            // create output messahe
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"%s United patterns %lu and %lu to %lu.\n               Calculation took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( long unsigned int ) aSourceA,
                        ( long unsigned int ) aSourceB,
                        ( long unsigned int ) aTarget,
                        ( double ) tElapsedTime / 1000 );

            }


        }
// -----------------------------------------------------------------------------

        /**
         * copies a source pattern to a target pattern
         */
        void
        HMR::copy_pattern(
                const uint & aSource,
                const uint & aTarget )
        {
            tic tTimer;

            mBackgroundMesh->copy_pattern(
                    aSource,
                    aTarget );

            this->update_meshes( aTarget );

            // create output messahe
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output

                std::fprintf( stdout,"%s Copied pattern %lu to %lu.\n               Calculation took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( long unsigned int ) aSource,
                        ( long unsigned int ) aTarget,
                        ( double ) tElapsedTime / 1000 );

            }


        }

// -----------------------------------------------------------------------------
        void
        HMR::finalize()
        {

            // synchronize flags for T-Matrices with other procs
            this->synchronize_t_matrix_flags();

            // remember active pattern
            auto tActivePattern = mBackgroundMesh->get_activation_pattern();


            // test if max polynomial is 3
            if ( mParameters->get_max_polynomial() > 2 )
            {
                // activate extra pattern for exodus
                this->add_extra_refinement_step_for_exodus();
            }

            // get number of Lagrange meshes
            uint tNumberOfLagrangeMeshes = mLagrangeMeshes.size();

            // loop over all meshes
            for( uint l=0; l<tNumberOfLagrangeMeshes; ++l )
            {
                mTMatrix( l )->evaluate();
            }

            // create communication table
            this->create_communication_table();

            for( auto tMesh : mBSplineMeshes )
            {
                tMesh->calculate_basis_indices( mCommunicationTable );
            }

            for( auto tMesh: mLagrangeMeshes )
            {
                tMesh->calculate_node_indices();
            }



            // reset active pattern
            if ( mBackgroundMesh->get_activation_pattern() != tActivePattern )
            {
                mBackgroundMesh->set_activation_pattern( tActivePattern );
            }



        }

// -----------------------------------------------------------------------------

        void
        HMR::create_communication_table()
        {
            uint tParSize = par_size();
            uint tMyRank  = par_rank();

            if( tParSize > 1 )
            {
                // in a first step, we identify all processers this proc wants
                // to talk to

                // this is a Bool-like matrix
                Mat< uint > tColumn( tParSize, 1, 0 );

                // test owners of B-Splines
                for( auto tMesh: mBSplineMeshes )
                {
                    // get number of active B-Splines
                    auto tNumberOfBSplines = tMesh->get_number_of_active_basis_on_proc();

                    // loop over all active basis on this mesh
                    for( uint k=0; k<tNumberOfBSplines; ++k )
                    {
                        // get pointer to basis
                        auto tBasis = tMesh->get_active_basis( k );

                        // test if flag of basis is set
                        if ( tBasis->is_flagged() )
                        {
                            // set flag for this proc
                            tColumn( tBasis->get_owner() ) = 1;
                        }
                    }
                }

                // remove self from row
                tColumn( tMyRank ) = 0;

                // communication table
                Mat< uint > tCommTable;

                // matrices to send
                Cell< Mat< uint > > tSend;

                // matrices to receive
                Cell< Mat< uint > > tRecv;

                if( tMyRank != 0 )
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
                    for( uint k=1; k<tParSize; ++k )
                    {
                        tCommTable( k ) = k;
                    }

                    // nothing to send
                    Mat< uint > tEmpty;
                    tSend.resize( tParSize, tEmpty );
                }

                // exchange matrices
                communicate_mats( tCommTable, tSend, tRecv );

                // process information on master proc
                if ( tMyRank == 0 )
                {
                    // create communication matrix
                    Mat< uint > tCommMatrix( tParSize, tParSize, 0 );

                    // process first row
                    tRecv( 0 ) = tColumn;

                    // loop over all procs and create comm matrix
                    for( uint j=0; j<tParSize; ++j )
                    {
                        for( uint i=0; i<tParSize; ++i )
                        {
                            if ( tRecv( j )( i, 0 ) != 0 )
                            {
                                tCommMatrix( i, j ) = 1;
                                tCommMatrix( j, i ) = 1;
                            }
                        }
                    }

                    // remove diagonal
                    for( uint i=0; i<tParSize; ++i )
                    {
                        tCommMatrix( i, i ) = 0;
                    }

                    // create sending list
                    Mat< uint > tEmpty;
                    tSend.resize( tParSize, tEmpty );

                    for( uint j=0; j<tParSize; ++j )
                    {
                        // count nonzero entries
                        uint tCount = 0;
                        for( uint i=0; i<tParSize; ++i )
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
                        for( uint i=0; i<tParSize; ++i )
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
            else // if run in serial
            {
                // communication table is empty
                mCommunicationTable.set_size( 0, 1 );
            }

        }
// -----------------------------------------------------------------------------

        /* Field *
        HMR::create_field(
                     const std::string & aLabel,
                     const uint        & aLagrangeIndex )
        {
            // create new field
            Field * aField = new Field(
                    this,
                    aLabel,
                    aLagrangeIndex );

            return aField;
        } */

// -----------------------------------------------------------------------------

        /* void
        HMR::add_field( Field * aField )
        {
            // add to database
            mFields.push_back( aField );

        } */

// ---------------------------------------------------------------------------
        void
        HMR::save_to_exodus( const std::string & aPath )
        {

            this->HMR::save_to_exodus(
                    mParameters->get_output_pattern(),
                    aPath );

        }

// -----------------------------------------------------------------------------

       /**
         * this function is for testing purpose only. Data is always copied.
         * This is not an efficient way to do things!
         */
        void
        HMR::save_to_exodus( const uint & aBlock, const std::string & aPath )
        {

            /* mLagrangeMeshes( aBlock )->reset_fields();

           for( auto tField : mFields )
            {
                if ( tField->get_lagrange_index() == aBlock )
                {

                    mLagrangeMeshes( aBlock )->add_field(
                            tField->get_label(), tField->get_data() );
                }
            } */

            // create STK object
            STK * tSTK = mLagrangeMeshes( aBlock )->create_stk_object();

            // save MTK to exodus
            tSTK->save_to_file( aPath );

            // delete file
            delete tSTK;
        }

// -----------------------------------------------------------------------------

        void
        HMR::save_to_hdf5( const std::string & aPath )
        {
            // create file object
            File tHDF5;

            // create file on disk
            tHDF5.create( aPath );

            // store settings object
            tHDF5.save_settings( mParameters );

            // remember active pattern
            auto tActivePattern = mBackgroundMesh->get_activation_pattern();

            // save output pattern into file
            tHDF5.save_refinement_pattern(
                    mBackgroundMesh,
                    mParameters->get_output_pattern() );

            if( tActivePattern != mBackgroundMesh->get_activation_pattern() )
            {
                mBackgroundMesh->set_activation_pattern( tActivePattern );
            }

            // close hdf5 file
            tHDF5.close();
        }
// -----------------------------------------------------------------------------

        HMR::HMR( const std::string & aPath ) :
            mParameters( create_hmr_parameters_from_hdf5_file( aPath ) )
        {
            // create file object
            File tHDF5;

            // open file on disk
            tHDF5.open( aPath );

            // create factory
            Factory tFactory;

            // create background mesh object
            mBackgroundMesh = tFactory.create_background_mesh( mParameters );

            // reset all patterns
            for( uint k=0; k<gNumberOfPatterns; ++k )
            {
                mBackgroundMesh->reset_pattern( k );
            }

            // remember active pattern
            auto tActivePattern = mBackgroundMesh->get_activation_pattern();

            // load input pattern into file
            tHDF5.load_refinement_pattern( mBackgroundMesh, mParameters->get_input_pattern()  );

            // copy input pattern to output pattern
            this->copy_pattern( mParameters->get_input_pattern(), mParameters->get_output_pattern() );

            if( tActivePattern != mBackgroundMesh->get_activation_pattern() )
            {
                mBackgroundMesh->set_activation_pattern( tActivePattern );
            }


            // close hdf5 file
            tHDF5.close();

            // initialize mesh objects
            this->create_meshes();

            // initialize T-Matrix objects
            this->init_t_matrices();
        }
// -----------------------------------------------------------------------------

        /**
         * Project a source field to its target
         */
        /*   Field *
        HMR::map_field_to_output_mesh( Field * aField )
        {

            // start timer
            //tic tTimer;

            MORIS_ERROR( aField->get_activation_pattern() == mParameters->get_input_pattern(),
                    "Source Field must be on input pattern" );

            //this->set_activation_pattern( mParameters->get_input_pattern() );
            //this->activate_all_t_matrices(); // < -- this is a problem and needs fixing


            // create unity of both patterns
            this->unite_patterns(
                    mParameters->get_input_pattern(),
                    mParameters->get_output_pattern(),
                    mParameters->get_union_pattern() );

            // create a field on the union mesh
            Field * tUnion = this->create_field(
                    aField->get_label(),
                    mParameters->get_union_mesh( aField->get_order() ) );

            // project input field to output field
            this->interpolate_field( aField, tUnion );

            // create output mesh
            Field * aOutput = this->create_field(
                    aField->get_label(),
                    mParameters->get_output_mesh( aField->get_order() ) );

            // create mesh interface
            auto tMesh
                = this->create_mtk_interface( mParameters->get_union_mesh( aField->get_order() ) );

            // calculate coefficients for output mesh

            // create IWG object
            moris::fem::IWG_L2 tIWG;

            // create model
            mdl::Model tModel( tMesh, tIWG, tUnion->get_data(), aOutput->get_coefficients() );

            // evaluate result on output mesh
            aOutput->evaluate_node_values();

            return aOutput;
        } */
// -----------------------------------------------------------------------------

        /**
         * Project a source field to its target. Alternative function with error testing
         */
        /*       Field *
        HMR::map_field_to_output_mesh(
                Field * aField,
                real & aIntegrationError,
                real (*aFunction)( const Mat< real > & aPoint ) )
        {

            // start timer
            //tic tTimer;

            MORIS_ERROR( aField->get_activation_pattern() == mParameters->get_input_pattern(),
                    "Source Field must be on input pattern" );

            //this->set_activation_pattern( mParameters->get_input_pattern() );
            //this->activate_all_t_matrices(); // < -- this is a problem and needs fixing


            // create unity of both patterns
            this->unite_patterns(
                    mParameters->get_input_pattern(),
                    mParameters->get_output_pattern(),
                    mParameters->get_union_pattern() );

            // create a field on the union mesh
            Field * tUnion = this->create_field(
                    aField->get_label(),
                    mParameters->get_union_mesh( aField->get_order() ) );

            // project input field to output field
            this->interpolate_field( aField, tUnion );
            //tUnion->evaluate_function( aFunction );

            // create output mesh
            Field * aOutput = this->create_field(
                    aField->get_label(),
                    mParameters->get_output_mesh( aField->get_order() ) );

            // create mesh interface
            auto tMesh
                = this->create_mtk_interface( mParameters->get_union_mesh( aField->get_order() ) );

            // calculate coefficients for output mesh

            // create IWG object
            moris::fem::IWG_L2 tIWG;

            // create model
            mdl::Model tModel( tMesh, tIWG, tUnion->get_data(), aOutput->get_coefficients() );



            // evaluate result on output mesh
            aOutput->evaluate_node_values();

            // calculate error
            aIntegrationError = tModel.compute_integration_error( aFunction );

            return aOutput;
        } */

// -----------------------------------------------------------------------------

        /**
         * aTarget must be a refined variant of aSource
         */
        void
        HMR::interpolate_field( Field * aSource, Field * aTarget )
        {

            // make sure that mesh orders match
            MORIS_ERROR(
                    aSource->get_interpolation_order()
                 == aTarget->get_interpolation_order(),
                    "Source and Target Field must have same interpolation order" );

            // make sure that both fields are scalar or of equal dimension
            MORIS_ERROR( aSource->get_number_of_dimensions() == aTarget->get_number_of_dimensions(),
                                "Source and Target Field must have same dimension" );

            // allocate memory for target values
            aTarget->get_node_values()->set_size(
                    aTarget->get_lagrange_mesh()->get_number_of_all_basis_on_proc(), 1 );

            // source mesh
            auto tSourceMesh = aSource->get_lagrange_mesh();

            // target mesh
            auto tTargetMesh = aTarget->get_lagrange_mesh();

            // target mesh index
            auto tTargetMeshIndex = tTargetMesh->get_index();

            // get source pattern
            auto tSourcePattern = tSourceMesh->get_activation_pattern();


            // unflag nodes on target
            tTargetMesh->unflag_all_basis();

            // number of elements on target mesh
            auto tNumberOfElements = tTargetMesh->get_number_of_elements();

            // number of nodes per element
            auto tNumberOfNodesPerElement = tTargetMesh->get_number_of_basis_per_element();

            // create unity matrix
            Mat< real > tEye( tNumberOfNodesPerElement, tNumberOfNodesPerElement, 0.0 );
            for( uint k=0; k<tNumberOfNodesPerElement; ++k )
            {
                tEye( k, k ) = 1.0;
            }

            // get values of source field
            Mat< real > & tSourceData = * aSource->get_node_values();
            Mat< real > & tTargetData = * aTarget->get_node_values();

            Mat< real > tElementSourceData( tNumberOfNodesPerElement, aSource->get_number_of_dimensions() );
            Mat< real > tElementTargetData( tNumberOfNodesPerElement, aTarget->get_number_of_dimensions() );

            // loop over all elements
            for( luint e=0; e<tNumberOfElements; ++e )
            {
                // get pointer to target element
                auto tTargetElement = tTargetMesh->get_element( e );

                // get backgrund element
                auto tBackgroundElement = tTargetElement->get_background_element();

                // initialize refinement Matrix
                Mat< real > tR( tEye );


                while( ! tBackgroundElement->is_active( tSourcePattern ) )
                {
                    // right multiply refinement matrix
                    tR = tR * mTMatrix( tTargetMeshIndex )->get_refinement_matrix(
                            tBackgroundElement->get_child_index() );

                    // jump to parent
                    tBackgroundElement = tBackgroundElement->get_parent();
                }

                // get pointer to source element
                auto tSourceElement = tSourceMesh->get_element_by_memory_index(
                        tBackgroundElement->get_memory_index() );

                // fill source data vector
                for( uint k=0; k<tNumberOfNodesPerElement; ++k )
                {
                    // get pointer to source node
                    auto tNode = tSourceElement->get_basis( k );
                    auto tIndex = tNode->get_index();

                    // copy data from source mesh
                    tElementSourceData.rows( k, k ) = tSourceData.rows( tIndex, tIndex );
                }

                // project data on target element
                //tElementTargetData = tR * tElementSourceData;

                // copy target data to target mesh
                for( uint k=0; k<tNumberOfNodesPerElement; ++k )
                {
                    // get pointer to target node
                    auto tNode = tTargetElement->get_basis( k );

                    // test if data has already been written to target
                    if ( ! tNode->is_flagged() )
                    {
                        // get node indes
                        auto tIndex = tNode->get_index();

                        // copy data to target
                        // tTargetData.rows( tIndex, tIndex ) = tElementTargetData.rows( k, k );

                        tTargetData.rows( tIndex, tIndex ) = tR.rows( k, k ) * tElementSourceData;

                        // flag this node
                        tNode->flag();
                    }
                }
            }
        }

// -----------------------------------------------------------------------------

        /**
         * Extract values from source and copy them to target.
         * Needed for testing
         * aSource must be a refined variant of aTarget
         */
        /* void
        HMR::extract_field( Field * aSource, Field* aTarget )
        {

            // link to input matrix
            Mat< real > & tSourceData = aSource->get_data();

            // link to output matrix
            Mat< real > & tTargetData = aTarget->get_data();

            // get pointer to input mesh
            auto tSource = aSource->get_mesh();

            // get pointer to output mesh
            auto tTarget = aTarget->get_mesh();

            tTarget->select_activation_pattern();

            // get number of elements on target
            auto tNumberOfElements = tTarget->get_number_of_elements();

            // assign memory of output matrix
            tTargetData.set_size( tTarget->get_number_of_nodes_on_proc(), 1 );

            // unflag all nodes on target
            tTarget->unflag_all_basis();

            // get number of nodes per element
            auto tNumberOfNodes = tTarget->get_number_of_basis_per_element();

            // loop over all elements on target
            for( luint e=0; e<tNumberOfElements; ++e )
            {
                // get pointer to element on target
                auto tTargetElement = tTarget->get_element( e );

                // get pointer to twin on source
                auto tSourceElement = tSource->get_element_by_memory_index(
                        tTargetElement->get_memory_index() );

                for( uint k=0; k<tNumberOfNodes; ++k )
                {
                    // get pointer to target basis
                    auto tTargetBasis = tTargetElement->get_basis( k );

                    // check if basis has been processed
                    if ( ! tTargetBasis->is_flagged() )
                    {
                        // get pointer to source basis
                        auto tSourceBasis = tSourceElement->get_basis( k );

                        // copy value from source to target
                        tTargetData( tTargetBasis->get_index() )
                        = tSourceData ( tSourceBasis->get_id() );

                        // flag target
                        tTargetBasis->flag();
                    }
                }
            }

        } */

// -----------------------------------------------------------------------------

        void
        HMR::add_extra_refinement_step_for_exodus()
        {
            // get refined pattern
            auto tPattern = mParameters->get_refined_output_pattern();

            // create refined pattern
            mBackgroundMesh->copy_pattern(
                    mParameters->get_output_pattern(),
                    tPattern );

            // activate output pattern
            mBackgroundMesh->set_activation_pattern( tPattern );

            // collect all active elements on background mesh
            mBackgroundMesh->collect_active_elements();

            // get number of elements
            luint tNumberOfElements = mBackgroundMesh->get_number_of_active_elements_on_proc();

            // flag all active elements
            for( luint e=0; e<tNumberOfElements; ++e )
            {
                mBackgroundMesh->get_element( e )->put_on_refinement_queue();
            }

            // perform refinement
            this->perform_refinement();
        }

// -----------------------------------------------------------------------------

        void
        HMR::flag_against_nodal_field(
                const Mat< real > & aNodalValues )
        {

            // number of nodes on input field
            uint tNumberOfNodes = aNodalValues.length();

            // number of meshes
            uint tNumberOfMeshes = mLagrangeMeshes.size();

            // index of target mesh
            uint tMeshIndex = tNumberOfMeshes;

            // find correct mesh
            for( uint k=0; k<tNumberOfMeshes; ++k )
            {
                if ( mLagrangeMeshes( k )->get_activation_pattern()
                      == mParameters->get_output_pattern()
                      && mLagrangeMeshes( k )->get_number_of_nodes_on_proc()
                         == tNumberOfNodes )
                {
                    tMeshIndex = k;
                    break;
                }
            }

            MORIS_ERROR( tMeshIndex < tNumberOfMeshes,
                    "Could not find any output mesh with matching number of nodes");

            // create refinement manager
            Refinement_Manager tRefMan( mLagrangeMeshes( tMeshIndex ) );

            // flag elements
            tRefMan.flag_against_nodal_field( aNodalValues );
        }

// -----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */
