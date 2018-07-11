#ifndef MORIS_LINALG_FN_PINV_HPP_
#define MORIS_LINALG_FN_PINV_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief
     * Moore-Penrose pseudo-inverse of matrix A
     * and overrides the default tolerance, max(m,n)*max_sv*%datum::eps.
     *
     * @note Deleted function because Eigen does not support the pseudo-inverse computation.
     *
     * @param[in] A Matrix.
     * @param[in] tol Tolerance.
     *
     * @return pinv(A).
     */
    template< typename T1, typename T2 >
    void
    pinv(
            const moris::Mat< T1 > & A,
            const T2 &               tol = 0.01 ) = delete;
}

#endif  /* MORIS_LINALG_FN_PINV_HPP_ */
