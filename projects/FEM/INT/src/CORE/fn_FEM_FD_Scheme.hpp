/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_FD_Scheme.hpp
 *
 */

#ifndef SRC_FEM_FN_FEM_FD_SCHEME_HPP_
#define SRC_FEM_FN_FEM_FD_SCHEME_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_FEM_Enums.hpp"

namespace moris
{
//------------------------------------------------------------------------------
    namespace fem
    {
//------------------------------------------------------------------------------
        inline
        void fd_scheme(
                enum fem::FDScheme_Type              aFDSchemeType,
                Vector< Vector< real > > & aFDScheme )
        {
            // set size for aFDScheme
            aFDScheme.resize( 3 );

            // switch on geometry type
            switch( aFDSchemeType )
            {
                case fem::FDScheme_Type::POINT_1_FORWARD :
                {
                    // fill points
                    aFDScheme( 0 ) = { { 0.0 }, { 1.0 } };

                    // fill coeffs
                    aFDScheme( 1 ) = { { -1.0 }, { 1.0 } };

                    // fill deno
                    aFDScheme( 2 ) = { 1.0 };

                    break;
                }
                case fem::FDScheme_Type::POINT_1_BACKWARD :
                {
                    // fill points
                    aFDScheme( 0 ) = { { 0.0 }, { -1.0 } };

                    // fill coeffs
                    aFDScheme( 1 ) = { { 1.0 }, { -1.0 } };

                    // fill deno
                    aFDScheme( 2 ) = { 1.0 };

                    break;
                }
                case fem::FDScheme_Type::POINT_3_CENTRAL :
                {
                    // fill points
                    aFDScheme( 0 ) = { { -1.0 }, { 1.0 } };

                    // fill coeffs
                    aFDScheme( 1 ) = { { -1.0 }, { 1.0 } };

                    // fill deno
                    aFDScheme( 2 ) = { 2.0 };

                    break;
                }
                case fem::FDScheme_Type::POINT_5 :
                {
                    // fill points
                    aFDScheme( 0 ) = { { -2.0 }, { -1.0 }, { 1.0 }, { 2.0 } };

                    // fill coeffs
                    aFDScheme( 1 ) = { { 1.0 }, { -8.0 }, { 8.0 }, { -1.0 } };

                    // fill deno
                    aFDScheme( 2 ) = { 12.0 };

                    break;
                }
                default:
                {
                    MORIS_ERROR( false, " fd_scheme - unknown fd scheme type ");
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_FN_FEM_FD_SCHEME_HPP_ */

