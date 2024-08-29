/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_equal_to.tpp
 *
 */

namespace moris {
// ----------------------------------------------------------------------------

template< class T >
inline
bool
equal_to(
                T           const & a,
                T           const & b,
                moris::real const & error_factor )
{
    return (
            a == b ||
            std::abs( a - b ) < std::abs( std::min( a, b ) ) *
            std::numeric_limits< T >::epsilon() * error_factor  ||
            std::abs( a - b ) < std::numeric_limits< T >::epsilon() * error_factor );
}

// ----------------------------------------------------------------------------

template<>
inline
bool
equal_to(
                int         const & a,
                int         const & b,
                moris::real const &  )
{
    return ( a == b  );
}

// ----------------------------------------------------------------------------

template<>
inline
bool
equal_to(
                uint        const & a,
                uint        const & b,
                moris::real const &  )
{
    return ( a == b  );
}

// ----------------------------------------------------------------------------

template<>
inline
bool
equal_to(
                luint       const & a,
                luint       const & b,
                moris::real const &  )
{
    return ( a == b  );
}

// ----------------------------------------------------------------------------

template<>
inline
bool
equal_to(
                long long unsigned int const & a,
                long long unsigned int const & b,
                moris::real            const &  )
{
    return ( a == b  );
}

// ----------------------------------------------------------------------------

template<>
inline
bool
equal_to(
                long long int const & a,
                long long int const & b,
                moris::real   const &  )
{
    return ( a == b  );
}

// ----------------------------------------------------------------------------

template< class T, class A >
inline
bool
equal_to(
                T           const & a,
                A           const & b,
                moris::real const & error_factor )
{
    return (
            a == (T) b ||
            std::abs( a - (T) b ) < std::abs( std::min( a, (T) b ) ) *
            std::numeric_limits< T >::epsilon() * error_factor ||
            std::abs( a - (T) b ) < std::numeric_limits< T >::epsilon() * error_factor );
}

// ----------------------------------------------------------------------------

template< >
inline
bool
equal_to(
                uint        const & a,
                int         const & b,
                moris::real const &  )
{
     return ( (int) a == b );
}

// ----------------------------------------------------------------------------

template< >
inline
bool
equal_to(
                uint        const & a,
                luint       const & b,
                moris::real const &  )
{
     return ( a == b );
}

// ----------------------------------------------------------------------------

template< >
inline
bool
equal_to(
                luint       const & a,
                uint        const & b,
                moris::real const &  )
{
     return ( a == b );
}

// ----------------------------------------------------------------------------

template< >
inline
bool
equal_to(
                uint                   const & a,
                long long unsigned int const & b,
                moris::real            const &  )
{
     return ( a == b );
}

// ----------------------------------------------------------------------------

template< >
inline
bool
equal_to(
                long long unsigned int const & a,
                uint                   const & b,
                moris::real            const &  )
{
     return ( a == b );
}

// ----------------------------------------------------------------------------

template< >
inline
bool
equal_to(
                long long unsigned int const & a,
                int                    const & b,
                moris::real            const &  )
{
     return ( (lint) a == (lint) b );
}

// ----------------------------------------------------------------------------

template< >
inline
bool
equal_to(
                luint       const & a,
                int         const & b,
                moris::real const &  )
{
     return ( (lint) a == (lint) b );
}

// ----------------------------------------------------------------------------

template< >
inline
bool
equal_to(
                uint        const & a,
                real        const & b,
                moris::real const &  )
{
     return ( (real) a == b );
}

// ----------------------------------------------------------------------------

template< >
inline
bool
equal_to(
                luint       const & a,
                real        const & b,
                moris::real const &  )
{
     return ( (real) a == b );
}
}

