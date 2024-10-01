/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Edge2.hpp
 *
 */

#pragma once

#include "cl_HMR_Lagrange_Edge.hpp"
#include "moris_typedefs.hpp"

namespace moris::hmr
{
// ----------------------------------------------------------------------------

    template <>
    inline
    mtk::Interpolation_Order Lagrange_Edge< 2 >::get_interpolation_order() const
    {
        return mtk::Interpolation_Order::LINEAR;
    }

// ----------------------------------------------------------------------------

    template <>
    inline
    void Lagrange_Edge< 2 >::copy_vertex_pointers()
    {
        // get pointer to leader
        Element * tLeader = mElements( mIndexOfLeader );

        // pick edges from corresponding edge
        switch ( this->get_index_on_leader() )
        {
            case( 0 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 1 );
                mVertices[ 1 ] = tLeader->get_basis( 0 );
                break;
            }
            case( 1 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 2 );
                mVertices[ 1 ] = tLeader->get_basis( 1 );
                break;
            }
            case( 2 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 3 );
                mVertices[ 1 ] = tLeader->get_basis( 2 );
                break;
            }
            case( 3 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 3 );
                mVertices[ 1 ] = tLeader->get_basis( 0 );
                break;
            }
            case( 4 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 0 );
                mVertices[ 1 ] = tLeader->get_basis( 4 );
                break;
            }
            case( 5 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 1 );
                mVertices[ 1 ] = tLeader->get_basis( 5 );
                break;
            }
            case( 6 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 2 );
                mVertices[ 1 ] = tLeader->get_basis( 6 );
                break;
            }
            case( 7 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 3 );
                mVertices[ 1 ] = tLeader->get_basis( 7 );
                break;
            }
            case( 8 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 5 );
                mVertices[ 1 ] = tLeader->get_basis( 4 );
                break;
            }
            case( 9 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 5 );
                mVertices[ 1 ] = tLeader->get_basis( 6 );
                break;
            }
            case( 10 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 6 );
                mVertices[ 1 ] = tLeader->get_basis( 7 );
                break;
            }
            case( 11 ) :
            {
                mVertices[ 0 ] = tLeader->get_basis( 7 );
                mVertices[ 1 ] = tLeader->get_basis( 4 );
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
