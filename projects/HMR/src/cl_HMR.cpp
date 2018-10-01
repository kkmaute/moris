/*
 * cl_HMR.cpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_times.hpp" //LINALG/src
#include "fn_trans.hpp" //LINALG/src
#include "fn_eye.hpp" //LINALG/src

#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

#include "cl_HMR.hpp" //HMR/src
#include "cl_HMR_Mesh.hpp" //HMR/src
#include "cl_HMR_STK.hpp" //HMR/src
#include "cl_HMR_File.hpp" //HMR/src
#include "cl_MDL_Model.hpp"
#include "cl_FEM_IWG_L2.hpp"
#include "cl_HMR_Field.hpp"          //HMR/src
#include "cl_HMR_Background_Element_Base.hpp"

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

            // remember buffer
            // uint tBuffer = mParameters->get_buffer_size();

            // reset buffer
            //mBackgroundMesh->set_buffer_size( 0 );

            // load input pattern into file
            tHDF5.load_refinement_pattern( mBackgroundMesh, mParameters->get_input_pattern()  );

            // reset buffer
            // mBackgroundMesh->set_buffer_size( tBuffer );

            // copy input pattern to output pattern
            //this->copy_pattern( mParameters->get_input_pattern(), mParameters->get_output_pattern() );

            // close hdf5 file
            tHDF5.close();

            // initialize mesh objects
            this->create_meshes();

            // initialize T-Matrix objects
            this->init_t_matrices();

            // activate input pattern
            this->set_activation_pattern( mParameters->get_input_pattern() );
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

        void
        HMR::update_meshes()
        {

            // remember active pattern
            auto tActivePattern = mBackgroundMesh->get_activation_pattern();

            // update all B-Spline meshes
            for( auto tMesh : mBSplineMeshes )
            {
                // synchronize mesh with background mesh
                tMesh->update_mesh();
            }

            // update all Lagrange meshes and link elements to their
            // B-Spline twins
            for( auto tMesh : mLagrangeMeshes )
            {
                // synchronize mesh with background mesh
                tMesh->update_mesh();
            }

            // reset pattern
            mBackgroundMesh->set_activation_pattern( tActivePattern );


            // update T-Matrices
            this->finalize();
        }

// -----------------------------------------------------------------------------

        // fixme: this funciton should be obsolete
        Mesh *
        HMR::create_mtk_interface()
        {
            return new Mesh(
                    *this,
                    mParameters->get_max_polynomial(),
                    mParameters->get_output_pattern() );
        }

// -----------------------------------------------------------------------------

        Mesh *
        HMR::create_mtk_interface( const uint & aActivationPattern )
        {
            return new Mesh(
                    *this,
                    mParameters->get_max_polynomial(),
                    aActivationPattern );
        }

// -----------------------------------------------------------------------------

        Mesh *
        HMR::create_input_mesh()
        {

            return new Mesh( *this,
                    mParameters->get_max_polynomial(),
                    mParameters->get_input_pattern() );
        }


// -----------------------------------------------------------------------------

        Mesh *
        HMR::create_output_mesh()
        {
            return new Mesh(
                    *this,
                    mParameters->get_max_polynomial(),
                    mParameters->get_output_pattern() );
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

            // this->update_meshes();

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

            // remember active pattern
            auto tActivePattern = mBackgroundMesh->get_activation_pattern();


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


            if( mParameters->get_number_of_dimensions() == 3 )
            {
                mBackgroundMesh->create_faces_and_edges();
            }
            else
            {
                mBackgroundMesh->create_facets();
            }

            for( auto tMesh: mLagrangeMeshes )
            {
                tMesh->calculate_node_indices();


                // only needed for output mesh
                if( mParameters->get_output_pattern() == tMesh->get_activation_pattern() )
                {
                    // create facets
                    tMesh->create_facets();

                    // create edges
                    if( mParameters->get_number_of_dimensions() == 3 )
                    {
                        tMesh->create_edges();
                    }
                }

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
            moris_id tParSize = par_size();
            moris_id tMyRank  = par_rank();

            if( tParSize > 1 )
            {
                // in a first step, we identify all processers this proc wants
                // to talk to

                // this is a Bool-like matrix
                Matrix< IdMat > tColumn( tParSize, 1, 0 );

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
                Matrix< IdMat > tCommTable;

                // matrices to send
                Cell< Matrix< IdMat > > tSend;

                // matrices to receive
                Cell< Matrix< IdMat > > tRecv;

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
                    for( moris_id k=1; k<tParSize; ++k )
                    {
                        tCommTable( k ) = k;
                    }

                    // nothing to send
                    Matrix< IdMat > tEmpty;
                    tSend.resize( tParSize, tEmpty );
                }

                // exchange matrices
                communicate_mats( tCommTable, tSend, tRecv );

                // process information on master proc
                if ( tMyRank == 0 )
                {
                    // create communication matrix
                    Matrix< IdMat > tCommMatrix( tParSize, tParSize, 0 );

                    // process first row
                    tRecv( 0 ) = tColumn;

                    // loop over all procs and create comm matrix
                    for( moris_id j=0; j<tParSize; ++j )
                    {
                        for( moris_id i=0; i<tParSize; ++i )
                        {
                            if ( tRecv( j )( i, 0 ) != 0 )
                            {
                                tCommMatrix( i, j ) = 1;
                                tCommMatrix( j, i ) = 1;
                            }
                        }
                    }

                    // remove diagonal
                    for( moris_id i=0; i<tParSize; ++i )
                    {
                        tCommMatrix( i, i ) = 0;
                    }

                    // create sending list
                    Matrix< IdMat > tEmpty;
                    tSend.resize( tParSize, tEmpty );

                    for( moris_id j=0; j<tParSize; ++j )
                    {
                        // count nonzero entries
                        uint tCount = 0;
                        for( moris_id i=0; i<tParSize; ++i )
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
                        for( moris_id i=0; i<tParSize; ++i )
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
                real (*aFunction)( const Matrix< DDRMat > & aPoint ) )
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
        HMR::interpolate_field(
                const uint       & aSourcePattern,
                const mtk::Field * aSource,
                const uint       & aTargetPattern,
                      mtk::Field * aTarget )
        {

            // make sure that mesh orders match
            MORIS_ERROR(
                    aSource->get_interpolation_order()
                 == aTarget->get_interpolation_order(),
                    "Source and Target Field must have same interpolation order" );

            // make sure that both fields are scalar or of equal dimension
            MORIS_ERROR( aSource->get_number_of_dimensions() == aTarget->get_number_of_dimensions(),
                                "Source and Target Field must have same dimension" );


            // get interpolation order
            uint tOrder = aSource->get_interpolation_order();

            // pointer to mesh that is linked to input field
            Lagrange_Mesh_Base * tSourceMesh = nullptr;

            // find pointer to input mesh
            for( Lagrange_Mesh_Base *  tMesh : mLagrangeMeshes )
            {
                if (   tMesh->get_order() == tOrder
                    && tMesh->get_activation_pattern() == aSourcePattern )
                {
                    tSourceMesh = tMesh;
                    break;
                }
            }

            // pointer to mesh that is linked to output field
            Lagrange_Mesh_Base * tTargetMesh = nullptr;

            // find pointer to output mesh
            for( Lagrange_Mesh_Base *  tMesh : mLagrangeMeshes )
            {
                if (   tMesh->get_order() == tOrder
                        && tMesh->get_activation_pattern() == aTargetPattern )
                {
                    tTargetMesh = tMesh;
                    break;
                }
            }


            tTargetMesh->select_activation_pattern();

            // unflag nodes on target
            tTargetMesh->unflag_all_basis();

            // number of elements on target mesh
            auto tNumberOfElements = tTargetMesh->get_number_of_elements();

            // number of nodes per element
            auto tNumberOfNodesPerElement = tTargetMesh->get_number_of_basis_per_element();

            // create unity matrix
            Matrix< DDRMat > tEye( tNumberOfNodesPerElement, tNumberOfNodesPerElement, 0.0 );
            for( uint k=0; k<tNumberOfNodesPerElement; ++k )
            {
                tEye( k, k ) = 1.0;
            }

            // get values of source field
            const Matrix< DDRMat > & tSourceData = aSource->get_node_values();

            // get target data
            Matrix< DDRMat > & tTargetData = aTarget->get_node_values();

            // allocate value matrix
            tTargetData.set_size( tTargetMesh->get_number_of_all_basis_on_proc(), aTarget->get_number_of_dimensions() );

            // containers for source and target data
            Matrix< DDRMat > tElementSourceData( tNumberOfNodesPerElement, aSource->get_number_of_dimensions() );

            // target mesh index
            auto tTargetMeshIndex = tTargetMesh->get_index();




            // loop over all elements
            for( luint e=0; e<tNumberOfElements; ++e )
            {
                // get pointer to target element
                auto tTargetElement = tTargetMesh->get_element( e );

                // get backgrund element
                auto tBackgroundElement = tTargetElement->get_background_element();

                // initialize refinement Matrix
                Matrix< DDRMat > tR( tEye );

                while( ! tBackgroundElement->is_active( aSourcePattern ) )
                {
                    // right multiply refinement matrix
                    tR = tR.matrix_data() * mTMatrix( tTargetMeshIndex )->get_refinement_matrix(
                            tBackgroundElement->get_child_index() ).matrix_data();

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
                    tElementSourceData.set_row( k, tSourceData.get_row( tIndex ) );
                }

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

                        tTargetData.set_row( tIndex, tR.get_row( k ) * tElementSourceData );

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
            Matrix< DDRMat > & tSourceData = aSource->get_data();

            // link to output matrix
            Matrix< DDRMat > & tTargetData = aTarget->get_data();

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
            mBackgroundMesh->perform_refinement();
        }

// -----------------------------------------------------------------------------

        void
        HMR::flag_elements(
                      Cell< mtk::Cell* > & aElements,
                const uint                 aPattern )
        {

            // get  working pattern
            uint tWorkingPattern = mParameters->get_working_pattern();

            // check aPattern was set
            if( aPattern == MORIS_UINT_MAX )
            {
                // if default value is set, use output pattern
                mBackgroundMesh->set_activation_pattern( mParameters->get_input_pattern() );
            }
            else
            {
                // activate specified pattern on Background mesh
                mBackgroundMesh->set_activation_pattern( aPattern );
            }

            // loop over all active elements
            for( mtk::Cell* tCell :  aElements )
            {
                // get pointer to Background Element
                Background_Element_Base * tElement
                    = mBackgroundMesh->get_element( tCell->get_index() );

                // put this element on the list
                tElement->set_refined_flag( tWorkingPattern );

                // also flag all parents
                while( tElement->get_level() > 0 )
                {
                    // get parent of this element
                    tElement = tElement->get_parent();

                    // set flag for parent of element
                    tElement->set_refined_flag( tWorkingPattern );
                }
            }
        }

// -----------------------------------------------------------------------------

        mtk::Field *
        HMR::map_field_on_mesh( mtk::Field * aField, Mesh* aMesh )
        {
            // create pointer to output field
            mtk::Field * aOutField = aMesh->create_field( aField->get_label() );

            // create a temporary union mesh
            Mesh * tUnionMesh = new Mesh( *this,
                    aField->get_interpolation_order(),
                    mParameters->get_union_pattern() );

            // create a temporary untion field
            mtk::Field * tUnionField = tUnionMesh->create_field( aField->get_label() );

            // interpolate input field to union
            this->interpolate_field(
                    mParameters->get_input_pattern(),
                    aField,
                    mParameters->get_union_pattern(),
                    tUnionField );

            // perform L2 projection :

            // create IWG object
            moris::fem::IWG_L2 tIWG;

            // create model
            mdl::Model tModel(
                    tUnionMesh,
                    tIWG,
                    tUnionField->get_node_values(),
                    aOutField->get_coefficients() );

            // when the L2 projeciton is done, the union field is not needed anymore
            delete tUnionField;

            // neither is the union mesh
            delete tUnionMesh;

            // finally, the node values on the output mesh are evaluated
            aOutField->evaluate_node_values();

            // return output mesh
            return aOutField;
        }

// -----------------------------------------------------------------------------

        void
        HMR::perform_refinement( const bool aResetPattern )
        {

            // get pointer to working pattern
            uint tWorkingPattern = mParameters->get_working_pattern();

            // this function resets the output pattern
            if ( aResetPattern )
            {
                mBackgroundMesh->reset_pattern( mParameters->get_output_pattern() );
            }

            // activate input pattern
            mBackgroundMesh->set_activation_pattern( mParameters->get_output_pattern() );

            // get max level on this mesh
            uint tMaxLevel = mBackgroundMesh->get_max_level();

            // loop over all levels
            for( uint l=0; l<= tMaxLevel; ++l )
            {
                // container for elements on this level
                Cell< Background_Element_Base* > tElementList;

                // collect all elements on this level ( without aura )
                mBackgroundMesh
                    ->collect_elements_on_level( l, tElementList );

                // loop over all elements and flag them for refinement
                for( Background_Element_Base* tElement : tElementList )
                {
                    // test if element is marked on working pattern
                    if ( tElement->is_refined( tWorkingPattern ) )
                    {
                        // flag this element
                        tElement->put_on_refinement_queue();
                    }
                }
                mBackgroundMesh->perform_refinement();
            }

            // update max level
            tMaxLevel = mBackgroundMesh->get_max_level();

            // tidy up element flags
            for( uint l=0; l<= tMaxLevel; ++l )
            {
                // container for elements on this level
                Cell< Background_Element_Base* > tElementList;

                // collect all elements on this level ( with aura )
                mBackgroundMesh
                    ->collect_elements_on_level_including_aura( l, tElementList );

            }

            // tidy up working pattern
            mBackgroundMesh->reset_pattern( tWorkingPattern );

            // test if max polynomial is 3
            if ( mParameters->get_max_polynomial() > 2 )
            {
                // activate extra pattern for exodus
                this->add_extra_refinement_step_for_exodus();
            }

            // create union of input and output
            this->create_union_pattern();

            // update meshes according to new refinement patterns
            this->update_meshes();
        }

// -----------------------------------------------------------------------------

        uint
        HMR::flag_volume_and_surface_elements( const mtk::Field * aScalarField )
        {
            // the funciton returns the number of flagged elements
            uint aElementCounter = 0;

            // create geometry engine
            gen::Geometry_Engine tRefMan;

            // candidates for refinement
            Cell< mtk::Cell* > tCandidates;

            // elements to be flagged for refinement
            Cell< mtk::Cell* > tRefinementList;

            // get candidates from volume
            this->get_candidates_for_volume_refinement( tCandidates );

            // call refinement manager and get volume cells
            tRefMan.find_cells_within_levelset(
                    tRefinementList,
                    tCandidates,
                    aScalarField );

            // add length of list to counter
            aElementCounter += tRefinementList.size();

            // flag elements in hmr
            this->flag_elements( tRefinementList );

            // get candidates for surface
            this->get_candidates_for_surface_refinement( tCandidates );


            // call refinement manager and get intersected cells
            tRefMan.find_cells_intersected_by_levelset(
                    tRefinementList,
                    tCandidates,
                    aScalarField );

            // add length of list to counter
            aElementCounter += tRefinementList.size();

            // flag elements in hmr
            this->flag_elements( tRefinementList );


            // return number of flagged elements
            return aElementCounter;
        }

// -----------------------------------------------------------------------------

        void
        HMR::get_candidates_for_refinement(
                Cell< mtk::Cell* >   & aCandidates,
                const             uint aMaxLevel )
        {
            // reset candidate list
            aCandidates.clear();

            // make sure that input pattern is active
            this->set_activation_pattern( mParameters->get_input_pattern() );

            // get working pattern
            uint tWorkingPattern = mParameters->get_working_pattern();

            // number of active elements
            uint tNumberOfElements
                = mBackgroundMesh->get_number_of_active_elements_on_proc();

            // allocate output list
            aCandidates.resize(
                    tNumberOfElements,
                nullptr );

            // initialize counter
            uint tCount = 0;

            // pick first lagrange mesh on input pattern
            // fixme: add option to pick another one
            Lagrange_Mesh_Base * tMesh = mLagrangeMeshes( 0 );

            // make sure that this mesh uses correct pattern
            MORIS_ASSERT( tMesh->get_activation_pattern() ==  mParameters->get_input_pattern(),
                    "wrong pattern picked for get_candidates_for_refinement()");

            // loop over all elements
            for( uint e=0; e<tNumberOfElements; ++e )
            {
                // get pointer to background element
                Background_Element_Base * tElement
                    = mBackgroundMesh->get_element( e );

                // test if element is not flagged and below level
                if ( tElement->get_level() < aMaxLevel
                        && ! tElement->is_refined( tWorkingPattern ) )
                {
                    // add element to queue
                    aCandidates( tCount++ )
                            = tMesh->get_element_by_memory_index( tElement->get_memory_index() );
                }
            }

            // shrink output cell to fit
            aCandidates.resize( tCount );
        }

// -----------------------------------------------------------------------------

        void
        HMR::get_candidates_for_volume_refinement( Cell< mtk::Cell* > & aCandidates )
        {
            this->get_candidates_for_refinement(
                    aCandidates, mParameters->get_max_volume_level() );
        }

// -----------------------------------------------------------------------------

        void
        HMR::get_candidates_for_surface_refinement( Cell< mtk::Cell* > & aCandidates )
        {
            this->get_candidates_for_refinement(
                    aCandidates, mParameters->get_max_surface_level() );
        }

// -----------------------------------------------------------------------------

        void
        HMR::save_coeffs_to_binary_files(
                const std::string & aFilePath )
        {

            // loop over all meshes
            for( Lagrange_Mesh_Base * tMesh : mLagrangeMeshes )
            {
                // test if mesh links to output pattern
                if( tMesh->get_activation_pattern() == mParameters->get_output_pattern() )
                {
                    // calculate file path
                    std::string tFilePath =
                            aFilePath.substr(0,aFilePath.find_last_of(".")) // base path

                            // get order of lagrange mesh
                            + "_" + std::to_string( tMesh->get_order() )

                            // get order of bspline mesh
                            + "_" + std::to_string( tMesh->get_bspline_order() )

                            // finish path
                            +  aFilePath.substr( aFilePath.find_last_of("."), aFilePath.length() );

                    tMesh->save_coeffs_to_binary_file( tFilePath );
                }
            }
        }

// -----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */
