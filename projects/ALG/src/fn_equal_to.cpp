#include "fn_equal_to.hpp"

// ----------------------------------------------------------------------------

template<>
bool
moris::equal_to<moris::cplx>(
        moris::cplx const & a,
        moris::cplx const & b,
        moris::real const & error_factor )
{
    moris::real a_t = a.real();
    moris::real b_t = b.real();

    return (
            a_t == b_t ||
            std::abs( a_t - b_t ) < std::abs( std::min( a_t, b_t ) ) *
            std::numeric_limits< moris::real >::epsilon() * error_factor );
}
