
// ----------------------------------------------------------------------------

template< class T >
bool
moris::equal_to(
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

template< class T, class A >
bool
moris::equal_to(
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
