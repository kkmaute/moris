/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Edge3.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_EDGE3_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_EDGE3_HPP_

#include "cl_HMR_Lagrange_Edge.hpp"
#include "typedefs.hpp"

namespace moris::hmr
{
// ----------------------------------------------------------------------------

    template <>
    inline
    mtk::Interpolation_Order Lagrange_Edge< 3 >::get_interpolation_order() const
    {
        return mtk::Interpolation_Order::QUADRATIC;
    }

// ----------------------------------------------------------------------------

    template <>
    inline
    void Lagrange_Edge< 3 >::copy_vertex_pointers()
    {
        // get pointer to leader
        Element * tLeader = mElements( mIndexOfLeader );

        // pick edges from corresponding edge
        switch ( this->get_index_on_leader() )
        {
            case( 0 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 1 );
                mVertices[ 1 ] = tLeader->get_basis_function( 0 );
                mVertices[ 2 ] = tLeader->get_basis_function( 8 );
                break;
            }
            case( 1 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 2 );
                mVertices[ 1 ] = tLeader->get_basis_function( 1 );
                mVertices[ 2 ] = tLeader->get_basis_function( 9 );
                break;
            }
            case( 2 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 3 );
                mVertices[ 1 ] = tLeader->get_basis_function( 2 );
                mVertices[ 2 ] = tLeader->get_basis_function( 10 );
                break;
            }
            case( 3 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 3 );
                mVertices[ 1 ] = tLeader->get_basis_function( 0 );
                mVertices[ 2 ] = tLeader->get_basis_function( 11 );
                break;
            }
            case( 4 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 0 );
                mVertices[ 1 ] = tLeader->get_basis_function( 4 );
                mVertices[ 2 ] = tLeader->get_basis_function( 12 );
                break;
            }
            case( 5 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 1 );
                mVertices[ 1 ] = tLeader->get_basis_function( 5 );
                mVertices[ 2 ] = tLeader->get_basis_function( 13 );
                break;
            }
            case( 6 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 2 );
                mVertices[ 1 ] = tLeader->get_basis_function( 6 );
                mVertices[ 2 ] = tLeader->get_basis_function( 14 );
                break;
            }
            case( 7 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 3 );
                mVertices[ 1 ] = tLeader->get_basis_function( 7 );
                mVertices[ 2 ] = tLeader->get_basis_function( 15 );
                break;
            }
            case( 8 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 5 );
                mVertices[ 1 ] = tLeader->get_basis_function( 4 );
                mVertices[ 2 ] = tLeader->get_basis_function( 16 );
                break;
            }
            case( 9 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 5 );
                mVertices[ 1 ] = tLeader->get_basis_function( 6 );
                mVertices[ 2 ] = tLeader->get_basis_function( 17 );
                break;
            }
            case( 10 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 6 );
                mVertices[ 1 ] = tLeader->get_basis_function( 7 );
                mVertices[ 2 ] = tLeader->get_basis_function( 18 );
                break;
            }
            case( 11 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis_function( 7 );
                mVertices[ 1 ] = tLeader->get_basis_function( 4 );
                mVertices[ 2 ] = tLeader->get_basis_function( 19 );
                break;
            }
            default :
            {
                 MORIS_ERROR( false, "Unknown edge index");
            }
        }
    }

// ----------------------------------------------------------------------------
} /* namespace moris */
#endif /* SRC_HMR_CL_HMR_LAGRANGE_EDGE3_HPP_ */

