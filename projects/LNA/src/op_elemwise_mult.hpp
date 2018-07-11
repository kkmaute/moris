#ifndef MORIS_LINALG_OP_ELEMWISE_MULT_HPP_
#define MORIS_LINALG_OP_ELEMWISE_MULT_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

namespace arma_Math
{
    template< typename T >
    auto
    operator%(
            const moris::Mat< T > & A,
            const moris::Mat< T > & B )
    -> decltype( A.data() % B.data() )
    {
        return A.data() % B.data();
    }
}

// ----------------------------------------------------------------------------

namespace eigen_Math
{
    template< typename T >
    auto
    operator%(
            const moris::Mat< T > & A,
            const moris::Mat< T > & B )
    -> decltype( A.data().cwiseProduct( B.data() ) )
    {
        return A.data().cwiseProduct( B.data() ); // cwiseProduct is an Eigen functionality for element-wise multiplication of two matrices
    }
}

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Element wise multiplication operator.
     *
     * @param[in] A Matrix.
     * @param[in] B Matrix.
     *
     * @return Creates a matrix corresponding to element-wise
     * multiplication of the two input matrices.
     *
     * Example:
     * @include LNA/src/op_elemwise_mult.inc
     */
    template< typename T >
    auto
    operator%(
            const moris::Mat< T > & A,
            const moris::Mat< T > & B )
    -> decltype( moris::Math::operator%( A, B ) )
    {
        return moris::Math::operator%( A, B );
    }
}

#endif  /* MORIS_LINALG_OP_ELEMWISE_MULT_HPP_ */
