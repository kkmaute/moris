/*
 * cl_HMR_Lagrange_Facet_Line2.hpp
 *
 *  Created on: September 25, 2018
 *  using MATLAB
 */
 
#ifndef SRC_HMR_CL_HMR_LAGRANGE_FACET_LINE2_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_FACET_LINE2_HPP_

#include "cl_HMR_Lagrange_Facet.hpp"

namespace moris
{
    namespace hmr
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
        Lagrange_Facet< 2, 2 >::copy_vertex_pointers( const uint & aIndex )
        {
            // pick side of parent element
            switch( aIndex )
            {
                case( 0 ) :
                {
                    mVertices[ 0 ] = mMaster->get_basis( 0 );
                    mVertices[ 1 ] = mMaster->get_basis( 1 );
                    break;
                }
                case( 1 ) :
                {
                    mVertices[ 0 ] = mMaster->get_basis( 1 );
                    mVertices[ 1 ] = mMaster->get_basis( 2 );
                    break;
                }
                case( 2 ) :
                {
                    mVertices[ 0 ] = mMaster->get_basis( 2 );
                    mVertices[ 1 ] = mMaster->get_basis( 3 );
                    break;
                }
                case( 3 ) :
                {
                    mVertices[ 0 ] = mMaster->get_basis( 3 );
                    mVertices[ 1 ] = mMaster->get_basis( 0 );
                    break;
                }
            }
        }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_FACET_LINE2_HPP_ */
