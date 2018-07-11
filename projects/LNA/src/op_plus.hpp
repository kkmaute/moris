#ifndef MORIS_LINALG_OP_PLUS_HPP_
#define MORIS_LINALG_OP_PLUS_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"
#include "cl_Sp_Mat.hpp"

namespace moris
{
    /**
     * @brief Addition of two arrays
     *
     * @param[in] aA Input matrix
     * @param[in] aB Input matrix
     *
     * This function returns a matrix such that @f$ \mathbf{C}_{ij} = \mathbf{A}_{ij} + \mathbf{B}_{ij} @f$
     * where C is sum of A and B. If A is an m-by-p and B is a m-by-p matrix, then C is an m-by-p matrix.
     *
     * Example:
     * @include LNA/src/op_plus.inc
     *
     */
    template< typename T1, typename T2 >
    auto
    operator+(
            moris::Base_Mat< T1 > const & aA,
            moris::Base_Mat< T2 > const & aB )
    ->decltype( aA.data() + aB.data() )
    {
        return aA.data() + aB.data();
    }

    template< typename T1, typename T2,
    bool M = ( moris::is_Mat< T2 >::value || moris::is_Sp_Mat< T2 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator+(
            moris::Base_Mat< T1 > const & aA,
            T2                    const & aB )
    ->decltype( aA.data() + aB )
    {
        return aA.data() + aB;
    }

    template< typename T1, typename T2,
    bool M = ( moris::is_Mat< T1 >::value || moris::is_Sp_Mat< T1 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator+(
            T1                    const & aA,
            moris::Base_Mat< T2 > const & aB )
    ->decltype( aA + aB.data() )
    {
        return aA + aB.data();
    }
}

#endif  /* MORIS_LINALG_OP_PLUS_HPP_ */
