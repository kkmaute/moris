/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Facet_Line2.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_FACET_LINE2_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_FACET_LINE2_HPP_

#include "cl_HMR_Lagrange_Facet.hpp"

namespace moris::hmr
{

// ----------------------------------------------------------------------------

    template<>
    inline
    mtk::Geometry_Type
    Lagrange_Facet< 2, 2 >::get_geometry_type() const
    {
        return mtk::Geometry_Type::LINE;
    }

// ----------------------------------------------------------------------------

    template<>
    inline
    mtk::Interpolation_Order
    Lagrange_Facet< 2, 2 >::get_interpolation_order() const
    {
        return mtk::Interpolation_Order::LINEAR;
    }

// ----------------------------------------------------------------------------

    template<>
    inline
    void
    Lagrange_Facet< 2, 2 >::copy_vertex_pointers( uint aIndex )
    {
        // pick side of parent element
        switch( aIndex )
        {
            case( 0 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 0 );
                mVertices[ 1 ] = mLeader->get_basis_function( 1 );
                break;
            }
            case( 1 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 1 );
                mVertices[ 1 ] = mLeader->get_basis_function( 2 );
                break;
            }
            case( 2 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 2 );
                mVertices[ 1 ] = mLeader->get_basis_function( 3 );
                break;
            }
            case( 3 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis_function( 3 );
                mVertices[ 1 ] = mLeader->get_basis_function( 0 );
                break;
            }
        }
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_FACET_LINE2_HPP_ */

