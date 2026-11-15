/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "fn_matrix_difference_bound.hpp"
#include "fn_trans.hpp"

namespace moris
{
    // Relative symmetry comparison for square real matrices, using the same
    // conservative Frobenius bounds as matrix_difference_bound().
    // Bound relative asymmetry ||A - A^T||_2 / ||A||_2. An exactly zero matrix
    // is symmetric; nonfinite values produce a failing bound rather than a pass.
    inline real matrix_symmetry_bound( const Matrix< DDRMat >& aMatrix )
    {
        MORIS_ASSERT( aMatrix.n_rows() == aMatrix.n_cols(),
                "matrix_symmetry_bound() requires a square matrix." );

        const Matrix< DDRMat > tTranspose = trans( aMatrix );
        return matrix_difference_bound( aMatrix, tTranspose, 1.0, 0.0 );
    }
}    // namespace moris
