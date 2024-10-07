/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Check.hpp
 *
 */

#ifndef SRC_FEM_FN_FEM_CHECK_HPP_
#define SRC_FEM_FN_FEM_CHECK_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"

//------------------------------------------------------------------------------
namespace moris::fem
{
    //------------------------------------------------------------------------------
    inline bool
    check(
            Matrix< DDRMat >& aMatrixCheck,
            Matrix< DDRMat >& aMatrixRef,
            real              aEpsilon,
            bool              aShowDifferences    = true,
            bool              aShowMaxDifferences = false,
            real              aFDtolerance        = -1.0 )
    {
        // set tolerance on reference value above which relative difference is computed
        const real tRelDiffThreshold = 1e-12;

        // check if FD tolerance is defined
        // aFDtolerance is the absolute error do to finite differencing,
        // i.e. the absolute error that will accepted by the check, even if the relative error check fails
        if ( aFDtolerance < 0.0 )
        {
            // if FD tolerance is not set by user, use epsilon as default
            aFDtolerance = aEpsilon;
        }

        // check that matrices to compare have same size
        MORIS_ERROR(
                ( aMatrixCheck.n_rows() == aMatrixRef.n_rows() ) && ( aMatrixCheck.n_cols() == aMatrixRef.n_cols() ),
                "fem::check() - matrices to check do not share same dimensions - check: %zu x %zu ref: %zu x %zu.",
                aMatrixCheck.n_rows(),
                aMatrixCheck.n_cols(),
                aMatrixRef.n_rows(),
                aMatrixRef.n_cols() );

        // define a boolean for check
        bool tCheckBool = true;

        // define a real for absolute difference
        real tAbsolute    = 0.0;
        real tAbsoluteMax = 0.0;

        // define a real for relative difference
        real tRelative    = 0.0;
        real tRelativeMax = 0.0;

        for ( uint ii = 0; ii < aMatrixCheck.n_rows(); ii++ )
        {
            for ( uint jj = 0; jj < aMatrixCheck.n_cols(); jj++ )
            {
                // get absolute difference
                tAbsolute    = std::abs( aMatrixCheck( ii, jj ) - aMatrixRef( ii, jj ) );
                tAbsoluteMax = std::max( tAbsoluteMax, tAbsolute );

                // get relative difference
                real tAbsReference = std::abs( aMatrixRef( ii, jj ) );

                tRelative = tAbsReference > tRelDiffThreshold ? tAbsolute / tAbsReference : 0.0;

                // store maximum of absolute and relative difference if larger than FD tolerance
                if ( tAbsolute > aFDtolerance )
                {
                    tRelativeMax = std::max( tRelativeMax, tRelative );
                }

                // update check value
                tCheckBool = tCheckBool and ( ( tAbsolute < aFDtolerance ) or ( tRelative < aEpsilon ) );

                // for debug
                if ( ( tAbsolute > aFDtolerance ) and ( tRelative > aEpsilon ) and aShowDifferences )
                {
                    std::cout << "ii " << ii << " - jj " << jj << "\n";
                    std::cout << "aMatrixCheck( ii, jj ) " << aMatrixCheck( ii, jj ) << "\n";
                    std::cout << "aMatrixRef( ii, jj )   " << aMatrixRef( ii, jj ) << "\n";
                    std::cout << "Absolute difference " << tAbsolute << "\n";
                    std::cout << "Relative difference " << tRelative << "\n"
                              << std::flush;
                }
            }
        }

        if ( aShowMaxDifferences )
        {
            std::cout << "Max absolute difference: " << tAbsoluteMax << " \n"
                      << std::flush;
            std::cout << "Max relative difference: " << tRelativeMax << " \n"
                      << std::flush;
        }

        return tCheckBool;
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_FN_FEM_CHECK_HPP_ */
