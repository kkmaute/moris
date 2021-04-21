/*
 * cl_HMR_Lagrange_Facet_Quad16.hpp
 *
 *  Created on: September 25, 2018
 *  using MATLAB
 */
 
#ifndef SRC_HMR_CL_HMR_LAGRANGE_FACET_QUAD16_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_FACET_QUAD16_HPP_

#include "cl_HMR_Lagrange_Facet.hpp"

namespace moris
{
    namespace hmr
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
        Lagrange_Facet< 3, 16 >::copy_vertex_pointers( const uint & aIndex )
        {
            // pick side of parent element
            switch( aIndex )
            {
                case( 0 ) :
                {
                    mVertices[ 0 ] = mMaster->get_basis( 0 );
                    mVertices[ 1 ] = mMaster->get_basis( 1 );
                    mVertices[ 2 ] = mMaster->get_basis( 5 );
                    mVertices[ 3 ] = mMaster->get_basis( 4 );
                    mVertices[ 4 ] = mMaster->get_basis( 8 );
                    mVertices[ 5 ] = mMaster->get_basis( 9 );
                    mVertices[ 6 ] = mMaster->get_basis( 16 );
                    mVertices[ 7 ] = mMaster->get_basis( 17 );
                    mVertices[ 8 ] = mMaster->get_basis( 25 );
                    mVertices[ 9 ] = mMaster->get_basis( 24 );
                    mVertices[ 10 ] = mMaster->get_basis( 13 );
                    mVertices[ 11 ] = mMaster->get_basis( 12 );
                    mVertices[ 12 ] = mMaster->get_basis( 36 );
                    mVertices[ 13 ] = mMaster->get_basis( 37 );
                    mVertices[ 14 ] = mMaster->get_basis( 38 );
                    mVertices[ 15 ] = mMaster->get_basis( 39 );
                    break;
                }
                case( 1 ) :
                {
                    mVertices[ 0 ] = mMaster->get_basis( 1 );
                    mVertices[ 1 ] = mMaster->get_basis( 2 );
                    mVertices[ 2 ] = mMaster->get_basis( 6 );
                    mVertices[ 3 ] = mMaster->get_basis( 5 );
                    mVertices[ 4 ] = mMaster->get_basis( 14 );
                    mVertices[ 5 ] = mMaster->get_basis( 15 );
                    mVertices[ 6 ] = mMaster->get_basis( 20 );
                    mVertices[ 7 ] = mMaster->get_basis( 21 );
                    mVertices[ 8 ] = mMaster->get_basis( 29 );
                    mVertices[ 9 ] = mMaster->get_basis( 28 );
                    mVertices[ 10 ] = mMaster->get_basis( 17 );
                    mVertices[ 11 ] = mMaster->get_basis( 16 );
                    mVertices[ 12 ] = mMaster->get_basis( 44 );
                    mVertices[ 13 ] = mMaster->get_basis( 45 );
                    mVertices[ 14 ] = mMaster->get_basis( 46 );
                    mVertices[ 15 ] = mMaster->get_basis( 47 );
                    break;
                }
                case( 2 ) :
                {
                    mVertices[ 0 ] = mMaster->get_basis( 2 );
                    mVertices[ 1 ] = mMaster->get_basis( 3 );
                    mVertices[ 2 ] = mMaster->get_basis( 7 );
                    mVertices[ 3 ] = mMaster->get_basis( 6 );
                    mVertices[ 4 ] = mMaster->get_basis( 18 );
                    mVertices[ 5 ] = mMaster->get_basis( 19 );
                    mVertices[ 6 ] = mMaster->get_basis( 22 );
                    mVertices[ 7 ] = mMaster->get_basis( 23 );
                    mVertices[ 8 ] = mMaster->get_basis( 31 );
                    mVertices[ 9 ] = mMaster->get_basis( 30 );
                    mVertices[ 10 ] = mMaster->get_basis( 21 );
                    mVertices[ 11 ] = mMaster->get_basis( 20 );
                    mVertices[ 12 ] = mMaster->get_basis( 48 );
                    mVertices[ 13 ] = mMaster->get_basis( 49 );
                    mVertices[ 14 ] = mMaster->get_basis( 50 );
                    mVertices[ 15 ] = mMaster->get_basis( 51 );
                    break;
                }
                case( 3 ) :
                {
                    mVertices[ 0 ] = mMaster->get_basis( 3 );
                    mVertices[ 1 ] = mMaster->get_basis( 0 );
                    mVertices[ 2 ] = mMaster->get_basis( 4 );
                    mVertices[ 3 ] = mMaster->get_basis( 7 );
                    mVertices[ 4 ] = mMaster->get_basis( 11 );
                    mVertices[ 5 ] = mMaster->get_basis( 10 );
                    mVertices[ 6 ] = mMaster->get_basis( 12 );
                    mVertices[ 7 ] = mMaster->get_basis( 13 );
                    mVertices[ 8 ] = mMaster->get_basis( 26 );
                    mVertices[ 9 ] = mMaster->get_basis( 27 );
                    mVertices[ 10 ] = mMaster->get_basis( 23 );
                    mVertices[ 11 ] = mMaster->get_basis( 22 );
                    mVertices[ 12 ] = mMaster->get_basis( 43 );
                    mVertices[ 13 ] = mMaster->get_basis( 40 );
                    mVertices[ 14 ] = mMaster->get_basis( 41 );
                    mVertices[ 15 ] = mMaster->get_basis( 42 );
                    break;
                }
                case( 4 ) :
                {
                    mVertices[ 0 ] = mMaster->get_basis( 3 );
                    mVertices[ 1 ] = mMaster->get_basis( 2 );
                    mVertices[ 2 ] = mMaster->get_basis( 1 );
                    mVertices[ 3 ] = mMaster->get_basis( 0 );
                    mVertices[ 4 ] = mMaster->get_basis( 19 );
                    mVertices[ 5 ] = mMaster->get_basis( 18 );
                    mVertices[ 6 ] = mMaster->get_basis( 15 );
                    mVertices[ 7 ] = mMaster->get_basis( 14 );
                    mVertices[ 8 ] = mMaster->get_basis( 9 );
                    mVertices[ 9 ] = mMaster->get_basis( 8 );
                    mVertices[ 10 ] = mMaster->get_basis( 10 );
                    mVertices[ 11 ] = mMaster->get_basis( 11 );
                    mVertices[ 12 ] = mMaster->get_basis( 33 );
                    mVertices[ 13 ] = mMaster->get_basis( 34 );
                    mVertices[ 14 ] = mMaster->get_basis( 35 );
                    mVertices[ 15 ] = mMaster->get_basis( 32 );
                    break;
                }
                case( 5 ) :
                {
                    mVertices[ 0 ] = mMaster->get_basis( 4 );
                    mVertices[ 1 ] = mMaster->get_basis( 5 );
                    mVertices[ 2 ] = mMaster->get_basis( 6 );
                    mVertices[ 3 ] = mMaster->get_basis( 7 );
                    mVertices[ 4 ] = mMaster->get_basis( 24 );
                    mVertices[ 5 ] = mMaster->get_basis( 25 );
                    mVertices[ 6 ] = mMaster->get_basis( 28 );
                    mVertices[ 7 ] = mMaster->get_basis( 29 );
                    mVertices[ 8 ] = mMaster->get_basis( 30 );
                    mVertices[ 9 ] = mMaster->get_basis( 31 );
                    mVertices[ 10 ] = mMaster->get_basis( 27 );
                    mVertices[ 11 ] = mMaster->get_basis( 26 );
                    mVertices[ 12 ] = mMaster->get_basis( 52 );
                    mVertices[ 13 ] = mMaster->get_basis( 53 );
                    mVertices[ 14 ] = mMaster->get_basis( 54 );
                    mVertices[ 15 ] = mMaster->get_basis( 55 );
                    break;
                }
            }
        }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_FACET_QUAD16_HPP_ */
