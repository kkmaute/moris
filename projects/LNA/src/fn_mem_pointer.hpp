#ifndef MORIS_LINALG_FN_MEM_POINTER_HPP_
#define MORIS_LINALG_FN_MEM_POINTER_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    T*
    mem_pointer(
            moris::Mat< T > & aA )
    {
        return ( aA.data() ).memptr();
    }

    template< typename T >
    const T*
    mem_pointer(
            const moris::Mat< T > & aA )
    {
        return ( aA.data() ).memptr();
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    T*
    mem_pointer(
            moris::Mat< T > & aA )
    {
        return ( aA.data() ).data();
    }

    template< typename T >
    const T*
    mem_pointer(
            const moris::Mat< T > & aA )
    {
        return ( aA.data() ).data();
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Converts a moris::Mat to a pointer
     *
     * @param[in] aA Matrix
     *
     * @return Pointer to values in matrix aA.
     *
     * This function can only be used with dense matrices. The pointers
     * are ordered column major.
     *
     */
    template< typename T >
    T*
    mem_pointer(
            moris::Mat< T > & aA )
    {
        return moris::Math::mem_pointer( aA );
    }

    template< typename T >
    const T*
    mem_pointer(
            const moris::Mat< T > & aA )
    {
        return moris::Math::mem_pointer( aA );
    }
}

#endif  /* MORIS_LINALG_FN_MEM_POINTER_HPP_ */
