/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Facet_Quad9.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_FACET_QUAD9_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_FACET_QUAD9_HPP_

#include "cl_HMR_Lagrange_Facet.hpp"

namespace moris::hmr
{

// ----------------------------------------------------------------------------

    template<>
    inline
    mtk::Geometry_Type
    Lagrange_Facet< 3, 9 >::get_geometry_type() const
    {
        return mtk::Geometry_Type::QUAD;
    }

// ----------------------------------------------------------------------------

    template<>
    inline
    mtk::Interpolation_Order
    Lagrange_Facet< 3, 9 >::get_interpolation_order() const
    {
        return mtk::Interpolation_Order::QUADRATIC;
    }

// ----------------------------------------------------------------------------

    template<>
    inline
    void
    Lagrange_Facet< 3, 9 >::copy_vertex_pointers( uint aIndex )
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
                mVertices[ 5 ] = mLeader->get_basis_function( 13 );
                mVertices[ 6 ] = mLeader->get_basis_function( 16 );
                mVertices[ 7 ] = mLeader->get_basis_function( 12 );
                mVertices[ 8 ] = mLeader->get_basis_function( 25 );
                break;
            }
            case( 1 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 1 );
                mVertices[ 1 ] = mLeader->get_basis_function( 2 );
                mVertices[ 2 ] = mLeader->get_basis_function( 6 );
                mVertices[ 3 ] = mLeader->get_basis_function( 5 );
                mVertices[ 4 ] = mLeader->get_basis_function( 9 );
                mVertices[ 5 ] = mLeader->get_basis_function( 14 );
                mVertices[ 6 ] = mLeader->get_basis_function( 17 );
                mVertices[ 7 ] = mLeader->get_basis_function( 13 );
                mVertices[ 8 ] = mLeader->get_basis_function( 24 );
                break;
            }
            case( 2 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 2 );
                mVertices[ 1 ] = mLeader->get_basis_function( 3 );
                mVertices[ 2 ] = mLeader->get_basis_function( 7 );
                mVertices[ 3 ] = mLeader->get_basis_function( 6 );
                mVertices[ 4 ] = mLeader->get_basis_function( 10 );
                mVertices[ 5 ] = mLeader->get_basis_function( 15 );
                mVertices[ 6 ] = mLeader->get_basis_function( 18 );
                mVertices[ 7 ] = mLeader->get_basis_function( 14 );
                mVertices[ 8 ] = mLeader->get_basis_function( 26 );
                break;
            }
            case( 3 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 3 );
                mVertices[ 1 ] = mLeader->get_basis_function( 0 );
                mVertices[ 2 ] = mLeader->get_basis_function( 4 );
                mVertices[ 3 ] = mLeader->get_basis_function( 7 );
                mVertices[ 4 ] = mLeader->get_basis_function( 11 );
                mVertices[ 5 ] = mLeader->get_basis_function( 12 );
                mVertices[ 6 ] = mLeader->get_basis_function( 19 );
                mVertices[ 7 ] = mLeader->get_basis_function( 15 );
                mVertices[ 8 ] = mLeader->get_basis_function( 23 );
                break;
            }
            case( 4 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 3 );
                mVertices[ 1 ] = mLeader->get_basis_function( 2 );
                mVertices[ 2 ] = mLeader->get_basis_function( 1 );
                mVertices[ 3 ] = mLeader->get_basis_function( 0 );
                mVertices[ 4 ] = mLeader->get_basis_function( 10 );
                mVertices[ 5 ] = mLeader->get_basis_function( 9 );
                mVertices[ 6 ] = mLeader->get_basis_function( 8 );
                mVertices[ 7 ] = mLeader->get_basis_function( 11 );
                mVertices[ 8 ] = mLeader->get_basis_function( 21 );
                break;
            }
            case( 5 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 4 );
                mVertices[ 1 ] = mLeader->get_basis_function( 5 );
                mVertices[ 2 ] = mLeader->get_basis_function( 6 );
                mVertices[ 3 ] = mLeader->get_basis_function( 7 );
                mVertices[ 4 ] = mLeader->get_basis_function( 16 );
                mVertices[ 5 ] = mLeader->get_basis_function( 17 );
                mVertices[ 6 ] = mLeader->get_basis_function( 18 );
                mVertices[ 7 ] = mLeader->get_basis_function( 19 );
                mVertices[ 8 ] = mLeader->get_basis_function( 22 );
                break;
            }
        }
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_FACET_QUAD9_HPP_ */

