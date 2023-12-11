/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Facet_Quad16.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_FACET_QUAD16_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_FACET_QUAD16_HPP_

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
                mVertices[ 0 ] = mLeader->get_basis_function( 0 );
                mVertices[ 1 ] = mLeader->get_basis_function( 1 );
                mVertices[ 2 ] = mLeader->get_basis_function( 5 );
                mVertices[ 3 ] = mLeader->get_basis_function( 4 );
                mVertices[ 4 ] = mLeader->get_basis_function( 8 );
                mVertices[ 5 ] = mLeader->get_basis_function( 9 );
                mVertices[ 6 ] = mLeader->get_basis_function( 16 );
                mVertices[ 7 ] = mLeader->get_basis_function( 17 );
                mVertices[ 8 ] = mLeader->get_basis_function( 25 );
                mVertices[ 9 ] = mLeader->get_basis_function( 24 );
                mVertices[ 10 ] = mLeader->get_basis_function( 13 );
                mVertices[ 11 ] = mLeader->get_basis_function( 12 );
                mVertices[ 12 ] = mLeader->get_basis_function( 36 );
                mVertices[ 13 ] = mLeader->get_basis_function( 37 );
                mVertices[ 14 ] = mLeader->get_basis_function( 38 );
                mVertices[ 15 ] = mLeader->get_basis_function( 39 );
                break;
            }
            case( 1 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 1 );
                mVertices[ 1 ] = mLeader->get_basis_function( 2 );
                mVertices[ 2 ] = mLeader->get_basis_function( 6 );
                mVertices[ 3 ] = mLeader->get_basis_function( 5 );
                mVertices[ 4 ] = mLeader->get_basis_function( 14 );
                mVertices[ 5 ] = mLeader->get_basis_function( 15 );
                mVertices[ 6 ] = mLeader->get_basis_function( 20 );
                mVertices[ 7 ] = mLeader->get_basis_function( 21 );
                mVertices[ 8 ] = mLeader->get_basis_function( 29 );
                mVertices[ 9 ] = mLeader->get_basis_function( 28 );
                mVertices[ 10 ] = mLeader->get_basis_function( 17 );
                mVertices[ 11 ] = mLeader->get_basis_function( 16 );
                mVertices[ 12 ] = mLeader->get_basis_function( 44 );
                mVertices[ 13 ] = mLeader->get_basis_function( 45 );
                mVertices[ 14 ] = mLeader->get_basis_function( 46 );
                mVertices[ 15 ] = mLeader->get_basis_function( 47 );
                break;
            }
            case( 2 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 2 );
                mVertices[ 1 ] = mLeader->get_basis_function( 3 );
                mVertices[ 2 ] = mLeader->get_basis_function( 7 );
                mVertices[ 3 ] = mLeader->get_basis_function( 6 );
                mVertices[ 4 ] = mLeader->get_basis_function( 18 );
                mVertices[ 5 ] = mLeader->get_basis_function( 19 );
                mVertices[ 6 ] = mLeader->get_basis_function( 22 );
                mVertices[ 7 ] = mLeader->get_basis_function( 23 );
                mVertices[ 8 ] = mLeader->get_basis_function( 31 );
                mVertices[ 9 ] = mLeader->get_basis_function( 30 );
                mVertices[ 10 ] = mLeader->get_basis_function( 21 );
                mVertices[ 11 ] = mLeader->get_basis_function( 20 );
                mVertices[ 12 ] = mLeader->get_basis_function( 48 );
                mVertices[ 13 ] = mLeader->get_basis_function( 49 );
                mVertices[ 14 ] = mLeader->get_basis_function( 50 );
                mVertices[ 15 ] = mLeader->get_basis_function( 51 );
                break;
            }
            case( 3 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 3 );
                mVertices[ 1 ] = mLeader->get_basis_function( 0 );
                mVertices[ 2 ] = mLeader->get_basis_function( 4 );
                mVertices[ 3 ] = mLeader->get_basis_function( 7 );
                mVertices[ 4 ] = mLeader->get_basis_function( 11 );
                mVertices[ 5 ] = mLeader->get_basis_function( 10 );
                mVertices[ 6 ] = mLeader->get_basis_function( 12 );
                mVertices[ 7 ] = mLeader->get_basis_function( 13 );
                mVertices[ 8 ] = mLeader->get_basis_function( 26 );
                mVertices[ 9 ] = mLeader->get_basis_function( 27 );
                mVertices[ 10 ] = mLeader->get_basis_function( 23 );
                mVertices[ 11 ] = mLeader->get_basis_function( 22 );
                mVertices[ 12 ] = mLeader->get_basis_function( 43 );
                mVertices[ 13 ] = mLeader->get_basis_function( 40 );
                mVertices[ 14 ] = mLeader->get_basis_function( 41 );
                mVertices[ 15 ] = mLeader->get_basis_function( 42 );
                break;
            }
            case( 4 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 3 );
                mVertices[ 1 ] = mLeader->get_basis_function( 2 );
                mVertices[ 2 ] = mLeader->get_basis_function( 1 );
                mVertices[ 3 ] = mLeader->get_basis_function( 0 );
                mVertices[ 4 ] = mLeader->get_basis_function( 19 );
                mVertices[ 5 ] = mLeader->get_basis_function( 18 );
                mVertices[ 6 ] = mLeader->get_basis_function( 15 );
                mVertices[ 7 ] = mLeader->get_basis_function( 14 );
                mVertices[ 8 ] = mLeader->get_basis_function( 9 );
                mVertices[ 9 ] = mLeader->get_basis_function( 8 );
                mVertices[ 10 ] = mLeader->get_basis_function( 10 );
                mVertices[ 11 ] = mLeader->get_basis_function( 11 );
                mVertices[ 12 ] = mLeader->get_basis_function( 33 );
                mVertices[ 13 ] = mLeader->get_basis_function( 34 );
                mVertices[ 14 ] = mLeader->get_basis_function( 35 );
                mVertices[ 15 ] = mLeader->get_basis_function( 32 );
                break;
            }
            case( 5 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 4 );
                mVertices[ 1 ] = mLeader->get_basis_function( 5 );
                mVertices[ 2 ] = mLeader->get_basis_function( 6 );
                mVertices[ 3 ] = mLeader->get_basis_function( 7 );
                mVertices[ 4 ] = mLeader->get_basis_function( 24 );
                mVertices[ 5 ] = mLeader->get_basis_function( 25 );
                mVertices[ 6 ] = mLeader->get_basis_function( 28 );
                mVertices[ 7 ] = mLeader->get_basis_function( 29 );
                mVertices[ 8 ] = mLeader->get_basis_function( 30 );
                mVertices[ 9 ] = mLeader->get_basis_function( 31 );
                mVertices[ 10 ] = mLeader->get_basis_function( 27 );
                mVertices[ 11 ] = mLeader->get_basis_function( 26 );
                mVertices[ 12 ] = mLeader->get_basis_function( 52 );
                mVertices[ 13 ] = mLeader->get_basis_function( 53 );
                mVertices[ 14 ] = mLeader->get_basis_function( 54 );
                mVertices[ 15 ] = mLeader->get_basis_function( 55 );
                break;
            }
        }
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_FACET_QUAD16_HPP_ */

