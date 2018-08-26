/*
 * cl_HMR_Data.cpp
 *
 *  Created on: Jun 25, 2018
 *      Author: messe
 */

#include "op_times.hpp"     //LNA/src
#include "fn_dot.hpp"     //LNA/src
#include "cl_HMR_Field.hpp" //HMR/src
#include "cl_HMR.cpp"        //HMR/src

#include "cl_MDL_Model.hpp"  // FEM/MDL/src
#include "cl_FEM_IWG_L2.hpp" // FEM/INT/src

namespace moris
{
    namespace hmr
        {
//------------------------------------------------------------------------------

        Field::Field(  HMR         * aHMR,
                const std::string  & aLabel,
                const uint         & aLagrangeMeshIndex ) :
                            mParameters( aHMR->get_parameters() ),
                            mHMR( aHMR ),
                            mMeshIndex( aLagrangeMeshIndex ),
                            mFieldIndex( aHMR->get_number_of_fields() ),
                            mMesh( aHMR->get_lagrange_mesh_by_index( aLagrangeMeshIndex ) ),
                            mTMatrix( aHMR->get_t_matrix( aLagrangeMeshIndex ) ),
                            mNodeValues( mMesh->create_field_data( aLabel ) ),
                            mLabel ( mMesh->get_field_label( mFieldIndex ) )
        {
            // assign datafield
            //
            //mCoefficients.set_size( mMesh->get_number_of_bsplines_on_proc(), 1, 0.0 );
        }

//-------------------------------------------------------------------------------

        /**
         * assigns memory for node values
         */
        void
        Field::allocate_node_values()
        {
            mNodeValues.set_size( mMesh->get_number_of_nodes_on_proc(), 1 );
        }

//------------------------------------------------------------------------------

        void
        Field::evaluate_function( real (*aFunction)( const Mat< real > & aPoint ) )
        {
            // make sure that correct pattern is selected
            mMesh->select_activation_pattern();

            // ask mesh for number of nodes
            luint tNumberOfNodes = mMesh->get_number_of_nodes_on_proc();

            mNodeValues.set_size( tNumberOfNodes, 1 );

            // loop over all nodes
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                // get node pointer
                auto tNode = mMesh->get_node_by_index( k );

                // evaluate function
                mNodeValues( k ) = aFunction( tNode->get_coords() );

            }
        }

//------------------------------------------------------------------------------

        void
        Field::evaluate_node_values( )
        {

            // start timer
            tic tTimer;

            // make sure that correct pattern is selected
            mMesh->select_activation_pattern();

            // ask mesh for number of nodes
            luint tNumberOfNodes = mMesh->get_number_of_nodes_on_proc();

            // ask mesh for number of elements
            luint tNumberOfElements = mMesh->get_number_of_elements();

            // flag all elements of interest
            for( luint e=0; e<tNumberOfElements; ++e )
            {
                mMesh->get_element( e )->set_t_matrix_flag();
            }

            mHMR->synchronize_t_matrix_flags();

            // FIXME : only scalar
            // set size of values matrix
            MORIS_ERROR( mDimension == 1, "n-dimension not supported yet" );

            mNodeValues.set_size( tNumberOfNodes, 1 );

            mTMatrix->evaluate();

            for( luint k=0; k<tNumberOfNodes; ++k )
            {
                // get pointer to node
                auto tNode = mMesh->get_node_by_index( k );

                // get PDOFs from node
                auto tBSplines = tNode->get_adof_pointers();

                // get T-Matrix
                Mat< real > tTMatrix( *tNode->get_t_matrix() );

                // get number of nodes
                uint tNumberOfCoeffs = tTMatrix.length();

                // fill coeffs vector
                Mat< real > tCoeffs( tNumberOfCoeffs, 1 );
                for( uint i=0; i<tNumberOfCoeffs; ++i )
                {
                    mtk::Vertex* tBSpline = tBSplines( i );

                    // get index of basis
                    auto tIndex = tBSpline->get_index();

                    tCoeffs( i ) = mCoefficients( tIndex );
                }

                // write value into solution
                mNodeValues( k ) = dot( tTMatrix, tCoeffs );
            }

            // create output messahe
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output

                std::fprintf( stdout,"%s Calculated values for field %s.\n               Calcuation took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        this->get_label().c_str(),
                        ( double ) tElapsedTime / 1000 );

            }

        }

//------------------------------------------------------------------------------

        /**
         * performs an L2 projection in order to calculate coefficients
         */
        void
        Field::l2_project_coefficients()
        {
            // start timer
            tic tTimer;

            // activate my pattern
            mHMR->set_active_pattern( mMesh->get_active_pattern() );



            // create mesh interface
            auto tMesh = mHMR->create_interface( mMesh->get_active_pattern() );

            // tell hmr to use all T-matrices
            // fixme: find out why this needs to be called
            //mHMR->activate_all_t_matrices();

            // create IWG object
            moris::fem::IWG_L2 tIWG;

            // create model
            mdl::Model tModel(
                    tMesh,
                    tIWG,
                    this->get_data(),
                    this->get_coefficients() );

            // create output messahe
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output

                std::fprintf( stdout,"%s Calculated coefficients for field %s.\n               L2 projection took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        this->get_label().c_str(),
                        ( double ) tElapsedTime / 1000 );

            }
        }

//-------------------------------------------------------------------------------

        //Field *
        //Field::l2_map_to_pattern( const uint & aPattern )
       // {

        //}
//-------------------------------------------------------------------------------

        } /* namespace hmr */
} /* namespace moris */
