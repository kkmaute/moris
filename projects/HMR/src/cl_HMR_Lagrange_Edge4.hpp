/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Edge4.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_EDGE4_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_EDGE4_HPP_

#include "cl_HMR_Lagrange_Edge.hpp"
#include "typedefs.hpp"

namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------

        template <>
        inline
        mtk::Interpolation_Order Lagrange_Edge< 4 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::CUBIC;
        }

// ----------------------------------------------------------------------------

        template <>
        inline
        void Lagrange_Edge< 4 >::copy_vertex_pointers()
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
                    mVertices[ 2 ] = tLeader->get_basis( 9 );
                    mVertices[ 3 ] = tLeader->get_basis( 8 );
                    break;
                }
                case( 1 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 2 );
                    mVertices[ 1 ] = tLeader->get_basis( 1 );
                    mVertices[ 2 ] = tLeader->get_basis( 15 );
                    mVertices[ 3 ] = tLeader->get_basis( 14 );
                    break;
                }
                case( 2 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 3 );
                    mVertices[ 1 ] = tLeader->get_basis( 2 );
                    mVertices[ 2 ] = tLeader->get_basis( 19 );
                    mVertices[ 3 ] = tLeader->get_basis( 18 );
                    break;
                }
                case( 3 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 3 );
                    mVertices[ 1 ] = tLeader->get_basis( 0 );
                    mVertices[ 2 ] = tLeader->get_basis( 11 );
                    mVertices[ 3 ] = tLeader->get_basis( 10 );
                    break;
                }
                case( 4 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 0 );
                    mVertices[ 1 ] = tLeader->get_basis( 4 );
                    mVertices[ 2 ] = tLeader->get_basis( 12 );
                    mVertices[ 3 ] = tLeader->get_basis( 13 );
                    break;
                }
                case( 5 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 1 );
                    mVertices[ 1 ] = tLeader->get_basis( 5 );
                    mVertices[ 2 ] = tLeader->get_basis( 16 );
                    mVertices[ 3 ] = tLeader->get_basis( 17 );
                    break;
                }
                case( 6 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 2 );
                    mVertices[ 1 ] = tLeader->get_basis( 6 );
                    mVertices[ 2 ] = tLeader->get_basis( 20 );
                    mVertices[ 3 ] = tLeader->get_basis( 21 );
                    break;
                }
                case( 7 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 3 );
                    mVertices[ 1 ] = tLeader->get_basis( 7 );
                    mVertices[ 2 ] = tLeader->get_basis( 22 );
                    mVertices[ 3 ] = tLeader->get_basis( 23 );
                    break;
                }
                case( 8 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 5 );
                    mVertices[ 1 ] = tLeader->get_basis( 4 );
                    mVertices[ 2 ] = tLeader->get_basis( 25 );
                    mVertices[ 3 ] = tLeader->get_basis( 24 );
                    break;
                }
                case( 9 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 5 );
                    mVertices[ 1 ] = tLeader->get_basis( 6 );
                    mVertices[ 2 ] = tLeader->get_basis( 28 );
                    mVertices[ 3 ] = tLeader->get_basis( 29 );
                    break;
                }
                case( 10 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 6 );
                    mVertices[ 1 ] = tLeader->get_basis( 7 );
                    mVertices[ 2 ] = tLeader->get_basis( 30 );
                    mVertices[ 3 ] = tLeader->get_basis( 31 );
                    break;
                }
                case( 11 ) :
                {
                    mVertices[ 0 ] = tLeader->get_basis( 7 );
                    mVertices[ 1 ] = tLeader->get_basis( 4 );
                    mVertices[ 2 ] = tLeader->get_basis( 27 );
                    mVertices[ 3 ] = tLeader->get_basis( 26 );
                    break;
                }
                default :
                {
                     MORIS_ERROR( false, "Unknown edge index");
                }
            }
        }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
#endif /* SRC_HMR_CL_HMR_LAGRANGE_EDGE4_HPP_ */

