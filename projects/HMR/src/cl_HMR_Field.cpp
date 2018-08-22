/*
 * cl_HMR_Data.cpp
 *
 *  Created on: Jun 25, 2018
 *      Author: messe
 */

#include "op_times.hpp"     //LNA/src
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
                            mMeshIndex( aLagrangeMeshIndex ),
                            mFieldIndex( aHMR->get_number_of_fields() ),
                            mMesh( aHMR->get_lagrange_mesh_by_index( aLagrangeMeshIndex ) ),
                            mNodeValues( mMesh->create_field_data( aLabel ) ),
                            mLabel ( mMesh->get_field_label( mFieldIndex ) )
        {
            // assign datafield
            mNodeValues.set_size( mMesh->get_number_of_nodes_on_proc(), 1, 0.0 );
        }

//------------------------------------------------------------------------------

        void
        Field::evaluate_function( real (*aFunction)( const Mat< real > & aPoint ) )
        {
            // ask mesh for number of nodes
            luint tNumberOfNodes = mMesh->get_number_of_nodes_on_proc();

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

        } /* namespace hmr */
} /* namespace moris */
