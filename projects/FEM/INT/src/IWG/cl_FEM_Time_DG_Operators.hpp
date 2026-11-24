/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#ifndef MORIS_CL_FEM_TIME_DG_OPERATORS_HPP
#define MORIS_CL_FEM_TIME_DG_OPERATORS_HPP

#include <cmath>
#include <initializer_list>
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_inv.hpp"

namespace moris::fem
{
    /** Energy-corrected DG(1), report Sections 7 and 17. Temporal coordinates
     * below are xi in [0,1]; MORIS field interpolators use tau = 2*xi-1.
     * Rows 0 and 1 are respectively the constant weak test and collocation.
     */
    struct Time_DG_Operators
    {
        Matrix< DDRMat > H;
        Matrix< DDRMat > G;
        Matrix< DDRMat > e;
        Matrix< DDRMat > A;
        Matrix< DDRMat > b;
        Matrix< DDRMat > GA;
        Matrix< DDRMat > Gb;

        Time_DG_Operators()
                : H( 2, 2, 0.0 )
                , G( 2, 2, 0.0 )
                , e( 2, 1, 0.0 )
        {
            // Generate the reduced rows by Gauss quadrature, Eq. (26)-(28).
            for ( real tXi : { 0.5 - 0.5 / std::sqrt( 3.0 ), 0.5 + 0.5 / std::sqrt( 3.0 ) } )
            {
                H( 0, 0 ) += 0.5 * ( 1.0 - tXi );
                H( 0, 1 ) += 0.5 * tXi;
                G( 0, 0 ) -= 0.5;
                G( 0, 1 ) += 0.5;
            }
            G( 0, 0 ) += 1.0;
            H( 1, 0 ) = 1.0;
            G( 1, 0 ) = -1.0 + 4.0;
            G( 1, 1 ) = 1.0;
            e( 0 )    = 1.0;
            e( 1 )    = 4.0;

            // Only a fixed 2x2 temporal inverse is used, never a spatial mass inverse.
            Matrix< DDRMat > tInverse = inv( H );
            A                         = tInverse * G;
            b                         = tInverse * e;
            GA                        = G * A;
            Gb                        = G * b;
        }

        static const Time_DG_Operators& get()
        {
            static const Time_DG_Operators tOperators;
            return tOperators;
        }
    };
}    // namespace moris::fem
#endif
