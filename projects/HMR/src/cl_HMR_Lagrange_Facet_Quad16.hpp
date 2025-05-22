/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Facet_Quad16.hpp
 *
 */

#pragma once

#include "cl_HMR_Lagrange_Facet.hpp"

namespace moris::hmr
{

// ----------------------------------------------------------------------------

    template<>
    inline
    mtk::Geometry_Type
    Lagrange_Facet< 3, 16 >::get_geometry_type() const
    {
        return mtk::Geometry_Type::QUAD;
    }

// ----------------------------------------------------------------------------

    template<>
    inline
    mtk::Interpolation_Order
    Lagrange_Facet< 3, 16 >::get_interpolation_order() const
    {
        return mtk::Interpolation_Order::CUBIC;
    }

// ----------------------------------------------------------------------------

    template<>
    inline
    void
    Lagrange_Facet< 3, 16 >::copy_vertex_pointers( uint aIndex )
    {
        // pick side of parent element
        switch( aIndex )
        {
            case( 0 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 0 );
                mVertices[ 1 ] = mLeader->get_basis( 1 );
                mVertices[ 2 ] = mLeader->get_basis( 5 );
                mVertices[ 3 ] = mLeader->get_basis( 4 );
                mVertices[ 4 ] = mLeader->get_basis( 8 );
                mVertices[ 5 ] = mLeader->get_basis( 9 );
                mVertices[ 6 ] = mLeader->get_basis( 16 );
                mVertices[ 7 ] = mLeader->get_basis( 17 );
                mVertices[ 8 ] = mLeader->get_basis( 25 );
                mVertices[ 9 ] = mLeader->get_basis( 24 );
                mVertices[ 10 ] = mLeader->get_basis( 13 );
                mVertices[ 11 ] = mLeader->get_basis( 12 );
                mVertices[ 12 ] = mLeader->get_basis( 36 );
                mVertices[ 13 ] = mLeader->get_basis( 37 );
                mVertices[ 14 ] = mLeader->get_basis( 38 );
                mVertices[ 15 ] = mLeader->get_basis( 39 );
                break;
            }
            case( 1 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 1 );
                mVertices[ 1 ] = mLeader->get_basis( 2 );
                mVertices[ 2 ] = mLeader->get_basis( 6 );
                mVertices[ 3 ] = mLeader->get_basis( 5 );
                mVertices[ 4 ] = mLeader->get_basis( 14 );
                mVertices[ 5 ] = mLeader->get_basis( 15 );
                mVertices[ 6 ] = mLeader->get_basis( 20 );
                mVertices[ 7 ] = mLeader->get_basis( 21 );
                mVertices[ 8 ] = mLeader->get_basis( 29 );
                mVertices[ 9 ] = mLeader->get_basis( 28 );
                mVertices[ 10 ] = mLeader->get_basis( 17 );
                mVertices[ 11 ] = mLeader->get_basis( 16 );
                mVertices[ 12 ] = mLeader->get_basis( 44 );
                mVertices[ 13 ] = mLeader->get_basis( 45 );
                mVertices[ 14 ] = mLeader->get_basis( 46 );
                mVertices[ 15 ] = mLeader->get_basis( 47 );
                break;
            }
            case( 2 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 2 );
                mVertices[ 1 ] = mLeader->get_basis( 3 );
                mVertices[ 2 ] = mLeader->get_basis( 7 );
                mVertices[ 3 ] = mLeader->get_basis( 6 );
                mVertices[ 4 ] = mLeader->get_basis( 18 );
                mVertices[ 5 ] = mLeader->get_basis( 19 );
                mVertices[ 6 ] = mLeader->get_basis( 22 );
                mVertices[ 7 ] = mLeader->get_basis( 23 );
                mVertices[ 8 ] = mLeader->get_basis( 31 );
                mVertices[ 9 ] = mLeader->get_basis( 30 );
                mVertices[ 10 ] = mLeader->get_basis( 21 );
                mVertices[ 11 ] = mLeader->get_basis( 20 );
                mVertices[ 12 ] = mLeader->get_basis( 48 );
                mVertices[ 13 ] = mLeader->get_basis( 49 );
                mVertices[ 14 ] = mLeader->get_basis( 50 );
                mVertices[ 15 ] = mLeader->get_basis( 51 );
                break;
            }
            case( 3 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 3 );
                mVertices[ 1 ] = mLeader->get_basis( 0 );
                mVertices[ 2 ] = mLeader->get_basis( 4 );
                mVertices[ 3 ] = mLeader->get_basis( 7 );
                mVertices[ 4 ] = mLeader->get_basis( 11 );
                mVertices[ 5 ] = mLeader->get_basis( 10 );
                mVertices[ 6 ] = mLeader->get_basis( 12 );
                mVertices[ 7 ] = mLeader->get_basis( 13 );
                mVertices[ 8 ] = mLeader->get_basis( 26 );
                mVertices[ 9 ] = mLeader->get_basis( 27 );
                mVertices[ 10 ] = mLeader->get_basis( 23 );
                mVertices[ 11 ] = mLeader->get_basis( 22 );
                mVertices[ 12 ] = mLeader->get_basis( 43 );
                mVertices[ 13 ] = mLeader->get_basis( 40 );
                mVertices[ 14 ] = mLeader->get_basis( 41 );
                mVertices[ 15 ] = mLeader->get_basis( 42 );
                break;
            }
            case( 4 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 3 );
                mVertices[ 1 ] = mLeader->get_basis( 2 );
                mVertices[ 2 ] = mLeader->get_basis( 1 );
                mVertices[ 3 ] = mLeader->get_basis( 0 );
                mVertices[ 4 ] = mLeader->get_basis( 19 );
                mVertices[ 5 ] = mLeader->get_basis( 18 );
                mVertices[ 6 ] = mLeader->get_basis( 15 );
                mVertices[ 7 ] = mLeader->get_basis( 14 );
                mVertices[ 8 ] = mLeader->get_basis( 9 );
                mVertices[ 9 ] = mLeader->get_basis( 8 );
                mVertices[ 10 ] = mLeader->get_basis( 10 );
                mVertices[ 11 ] = mLeader->get_basis( 11 );
                mVertices[ 12 ] = mLeader->get_basis( 33 );
                mVertices[ 13 ] = mLeader->get_basis( 34 );
                mVertices[ 14 ] = mLeader->get_basis( 35 );
                mVertices[ 15 ] = mLeader->get_basis( 32 );
                break;
            }
            case( 5 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 4 );
                mVertices[ 1 ] = mLeader->get_basis( 5 );
                mVertices[ 2 ] = mLeader->get_basis( 6 );
                mVertices[ 3 ] = mLeader->get_basis( 7 );
                mVertices[ 4 ] = mLeader->get_basis( 24 );
                mVertices[ 5 ] = mLeader->get_basis( 25 );
                mVertices[ 6 ] = mLeader->get_basis( 28 );
                mVertices[ 7 ] = mLeader->get_basis( 29 );
                mVertices[ 8 ] = mLeader->get_basis( 30 );
                mVertices[ 9 ] = mLeader->get_basis( 31 );
                mVertices[ 10 ] = mLeader->get_basis( 27 );
                mVertices[ 11 ] = mLeader->get_basis( 26 );
                mVertices[ 12 ] = mLeader->get_basis( 52 );
                mVertices[ 13 ] = mLeader->get_basis( 53 );
                mVertices[ 14 ] = mLeader->get_basis( 54 );
                mVertices[ 15 ] = mLeader->get_basis( 55 );
                break;
            }
        }
    }

// ----------------------------------------------------------------------------
} /* namespace moris */
