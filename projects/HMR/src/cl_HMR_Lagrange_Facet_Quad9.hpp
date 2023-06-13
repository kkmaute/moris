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
                mVertices[ 0 ] = mLeader->get_basis( 0 );
                mVertices[ 1 ] = mLeader->get_basis( 1 );
                mVertices[ 2 ] = mLeader->get_basis( 5 );
                mVertices[ 3 ] = mLeader->get_basis( 4 );
                mVertices[ 4 ] = mLeader->get_basis( 8 );
                mVertices[ 5 ] = mLeader->get_basis( 13 );
                mVertices[ 6 ] = mLeader->get_basis( 16 );
                mVertices[ 7 ] = mLeader->get_basis( 12 );
                mVertices[ 8 ] = mLeader->get_basis( 25 );
                break;
            }
            case( 1 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 1 );
                mVertices[ 1 ] = mLeader->get_basis( 2 );
                mVertices[ 2 ] = mLeader->get_basis( 6 );
                mVertices[ 3 ] = mLeader->get_basis( 5 );
                mVertices[ 4 ] = mLeader->get_basis( 9 );
                mVertices[ 5 ] = mLeader->get_basis( 14 );
                mVertices[ 6 ] = mLeader->get_basis( 17 );
                mVertices[ 7 ] = mLeader->get_basis( 13 );
                mVertices[ 8 ] = mLeader->get_basis( 24 );
                break;
            }
            case( 2 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 2 );
                mVertices[ 1 ] = mLeader->get_basis( 3 );
                mVertices[ 2 ] = mLeader->get_basis( 7 );
                mVertices[ 3 ] = mLeader->get_basis( 6 );
                mVertices[ 4 ] = mLeader->get_basis( 10 );
                mVertices[ 5 ] = mLeader->get_basis( 15 );
                mVertices[ 6 ] = mLeader->get_basis( 18 );
                mVertices[ 7 ] = mLeader->get_basis( 14 );
                mVertices[ 8 ] = mLeader->get_basis( 26 );
                break;
            }
            case( 3 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 3 );
                mVertices[ 1 ] = mLeader->get_basis( 0 );
                mVertices[ 2 ] = mLeader->get_basis( 4 );
                mVertices[ 3 ] = mLeader->get_basis( 7 );
                mVertices[ 4 ] = mLeader->get_basis( 11 );
                mVertices[ 5 ] = mLeader->get_basis( 12 );
                mVertices[ 6 ] = mLeader->get_basis( 19 );
                mVertices[ 7 ] = mLeader->get_basis( 15 );
                mVertices[ 8 ] = mLeader->get_basis( 23 );
                break;
            }
            case( 4 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 3 );
                mVertices[ 1 ] = mLeader->get_basis( 2 );
                mVertices[ 2 ] = mLeader->get_basis( 1 );
                mVertices[ 3 ] = mLeader->get_basis( 0 );
                mVertices[ 4 ] = mLeader->get_basis( 10 );
                mVertices[ 5 ] = mLeader->get_basis( 9 );
                mVertices[ 6 ] = mLeader->get_basis( 8 );
                mVertices[ 7 ] = mLeader->get_basis( 11 );
                mVertices[ 8 ] = mLeader->get_basis( 21 );
                break;
            }
            case( 5 ) :
            {
                mVertices[ 0 ] = mLeader->get_basis( 4 );
                mVertices[ 1 ] = mLeader->get_basis( 5 );
                mVertices[ 2 ] = mLeader->get_basis( 6 );
                mVertices[ 3 ] = mLeader->get_basis( 7 );
                mVertices[ 4 ] = mLeader->get_basis( 16 );
                mVertices[ 5 ] = mLeader->get_basis( 17 );
                mVertices[ 6 ] = mLeader->get_basis( 18 );
                mVertices[ 7 ] = mLeader->get_basis( 19 );
                mVertices[ 8 ] = mLeader->get_basis( 22 );
                break;
            }
        }
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_FACET_QUAD9_HPP_ */

