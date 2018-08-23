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
            mNodeValues.set_size( mMesh->get_number_of_nodes_on_proc(), 1, 0.0 );
            mCoefficients.set_size( mMesh->get_number_of_bsplines_on_proc(), 1, 0.0 );
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
        Field::evaluate_node_values()
        {


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

                std::cout << "Node " << tNode->get_index() << std::endl;

                // fill coeffs vector
                Mat< real > tCoeffs( tNumberOfCoeffs, 1 );
                for( uint i=0; i<tNumberOfCoeffs; ++i )
                {
                    mtk::Vertex* tBSpline = tBSplines( i );

                    // get index of basis
                    auto tIndex = tBSpline->get_index();

                    tCoeffs( i ) = mCoefficients( tIndex );

                    std::cout << "    Basis : id " << tBSpline->get_id() << " index " << tBSpline->get_index() << " T " << tTMatrix( i ) << std::endl;
                }
                std::cout  << std::endl;
                // write value into solution
                mNodeValues( k ) = dot( tTMatrix, tCoeffs );
            }

        }

//------------------------------------------------------------------------------

        } /* namespace hmr */
} /* namespace moris */
