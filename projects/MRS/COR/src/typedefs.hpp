#ifndef MORIS_TYPEDEFS_HPP_
#define MORIS_TYPEDEFS_HPP_

// C++ header files.
#include <cstddef>
#include <vector>
#include <complex>
#include <limits>
#include <cfloat>
#include <climits>
/**
 * @morisdescription
 */
namespace moris
{
    /**
     * Definitions for basic variable declarations.
     *
     * These typedefs allow the user to switch number precision.
     * For example, the user could use:
     * @code{.cpp}
     *     typedef float moris::real;
     * @endcode
     * On x86 systems, float is typically 4 bytes long
     *
     * @code{.cpp}
     *     typedef double moris::real;
     * @endcode

     * On x86 systems, moris::real is 8 bytes long.
     *
     * For some discussion on floats and doubles, see this StackOverflow
     *     <a href="http://stackoverflow.com/questions/14511910/switching-between-float-and-double-precision-at-compile-time">
     * Question</a>.
     * For a discussion on the difference between long doubles and doubles, see this
     *     <a href="http://stackoverflow.com/questions/14221612/difference-between-long-double-and-double-in-c-and-c">
     * Question</a>.
     * Finally, for a discussion on the different variants of int, see this
     *     <a href="http://stackoverflow.com/questions/9696660/what-is-the-difference-between-int-int16-int32-and-int64">
     * Question</a>.
     *
     *
     * @note
     * The typedef definition 'typedef double Real' has been commented out,
     * because switching between single and double precision would not be straightforward.
     * Trilinos and several other Third-party header files use double precision,
     * so they would not work with single precision.
     *
     * Trilinos uses long long instead of int for global IDs in 64-bit systems.
     * MORIS will use moris::uint for indexing arrays and for local IDs,
     * while int (for 32-bit) and long long (for 64-bit) for global IDs.
     *
     * MORIS will start typdef-ing the STL containers.
     * This is done in the hopes of generating a pseudo-language, a la armadillo.
     * Also, it is done in hopes of maintaining backwards compatibility,
     * so if we ever need to change the vector implementation,
     * the user is already using namespace moris.
     * However, this is an experimental feature and may not work.
     * If the vector implementation needs to be changed,
     * multiple functions need to be replicated,
     * which might prove to be challenging.
     */
    typedef std::size_t               size_t;
    typedef bool                      bool_t;
    typedef unsigned short int        short_uint;  // Short unsigned integer (generally should not be used)
    typedef long unsigned int         luint;

    typedef int                       moris_id;
    typedef int                       moris_index;

    const moris_index                 gNoIndex = -1;
    const moris_id                    gNoID    = -1;

#ifdef MORIS_USE_32BIT
    typedef int                       lint;  // long int in 64bit, int in 32bit.
    typedef int                       sint;  // signed int in 32bit and 64bit.
    typedef unsigned int              uint;  // unsigned int in 32bit and 64bit.
    typedef double                    real;  // real
    typedef std::complex<double>      cplx;  // complex
#elif MORIS_USE_64BIT
    typedef long int                  lint;
    typedef int                       sint;
    typedef long unsigned int         uint;
    typedef long double               real;
    typedef std::complex<long double> cplx;
#endif

#define MORIS_LUINT_MAX   ULONG_MAX
#define MORIS_ID_MAX      std::numeric_limits<moris::moris_id>::max()
#define MORIS_INDEX_MAX   std::numeric_limits<moris::moris_index>::max()

    // define maximum numbers
#if MORIS_USE_64BIT
#define MORIS_LINT_MAX   LONG_MAX
#define MORIS_UINT_MAX   ULONG_MAX
#define MORIS_REAL_MAX   LDBL_MAX
#define MORIS_REAL_MIN   LDBL_MIN
#define MORIS_SINT_MAX   INT_MAX
#define MORIS_REAL_EPS   LDBL_EPSILON
#else
#define MORIS_LINT_MAX   INT_MAX
#define MORIS_UINT_MAX   UINT_MAX
#define MORIS_REAL_MAX   DBL_MAX
#define MORIS_REAL_MIN   DBL_MIN
#define MORIS_SINT_MAX   INT_MAX
#define MORIS_REAL_EPS   DBL_EPSILON
#endif

    // define maximum capacity (in bytes) of moris containers
#define MORIS_MAX_CELL_CAPACITY                 1e8   // maximum allowable size of matrix in bytes
#define MORIS_CELL_RESIZE_CHECK_LIMIT           0.01  // fraction of max allowable cell size that triggers check for reasonable resize
#define MORIS_CELL_RESIZE_FRACTION_LIMIT        0.1   // lower limit on fraction of new over original cell size in check for reasonable resize
#define MORIS_CELL_UTILIZATION_FRACTION_LIMIT   0.75  // lower limit on fraction of cell size vs cell capacity at cell destruction
#define MORIS_CELL_IMPLICIT_RESIZE_CALL_LIMIT   5     // maximum allowable number of implicit resize calls as part of push_back
#define MORIS_CELL_RESIZE_CALL_LIMIT            500   // maximum allowable number of resize calls

#define MORIS_MAX_MATRIX_SIZE               1e8   // maximum allowable size of dense matrix in bytes
#define MORIS_MATRIX_RESIZE_CHECK_LIMIT     0.01  // fraction of max allowable matrix size that triggers check for reasonable resize
#define MORIS_MATRIX_RESIZE_FRACTION_LIMIT  0.1   // lower limit on fraction of new over original matrix size in check for reasonable resize
#define MORIS_MATRIX_RESIZE_CALL_LIMIT      100   // maximum allowable number of resize calls
}

/**
 * Wrapper for Armadillo functions.
 *
 */
namespace arma_Math
{
}

/**
 * Wrapper for Eigen functions.
 */
namespace eigen_Math
{
}

namespace moris
{
    /**
     * Alias for third-party linear algebra packages.
     */
#ifdef MORIS_USE_ARMA
    namespace Math = arma_Math;
#elif  MORIS_USE_EIGEN
    namespace Math = eigen_Math;
#endif
}

#endif /* MORIS_TYPEDEFS_HPP_ */
