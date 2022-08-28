/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * xtk_typedefs.hpp
 *
 */

#ifndef SRC_CORE_XTK_TYPEDEFS_HPP_
#define SRC_CORE_XTK_TYPEDEFS_HPP_

// C++ header files.
#include <cstddef>
#include <vector>
#include <complex>

/**
 * @xtkdescription
 */
namespace xtk
{
/**
 * Definitions for basic variable declarations.
 *
 * These typedefs allow the user to switch number precision.
 * For example, the user could use:
 * @code{.cpp}
 *     typedef float xtk::real;
 * @endcode
 * On x86 systems, float is typically 4 bytes long
 *
 * @code{.cpp}
 *     typedef double xtk::real;
 * @endcode

 * On x86 systems, xtk::real is 8 bytes long.
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
 * MORIS will use xtk::uint for indexing arrays and for local IDs,
 * while int (for 32-bit) and long long (for 64-bit) for global IDs.
 *
 * MORIS will start typdef-ing the STL containers.
 * This is done in the hopes of generating a pseudo-language, a la armadillo.
 * Also, it is done in hopes of maintaining backwards compatibility,
 * so if we ever need to change the vector implementation,
 * the user is already using namespace xtk.
 * However, this is an experimental feature and may not work.
 * If the vector implementation needs to be changed,
 * multiple functions need to be replicated,
 * which might prove to be challenging.
 */

    typedef std::size_t               size_t;
    typedef bool                      bool_t;

//#ifdef XTK_USE_32BIT
    typedef int                       lint;  // long int in 64bit, int in 32bit.
    typedef int                       sint;  // signed int in 32bit and 64bit.
    typedef unsigned int              uint;  // unsigned int in 32bit and 64bit.
    typedef double                    real;  // real
    typedef std::complex<double>      cplx;  // complex
//#elif XTK_USE_64BIT
//    typedef long int                  lint;
//    typedef int                       sint;
//    typedef long unsigned int         uint;
//    typedef long double               real;
//    typedef std::complex<long double> cplx;
//#endif
}

#ifdef XTK_USE_EIGEN
#include <Eigen/Dense>
    typedef Eigen::Matrix<xtk::real, Eigen::Dynamic, Eigen::Dynamic> moris::DDRMat;
    typedef Eigen::Matrix<xtk::size_t, Eigen::Dynamic, Eigen::Dynamic> moris::DDSTMat;
    typedef Eigen::Matrix<xtk::sint, Eigen::Dynamic, Eigen::Dynamic> moris::DDSTMat_Sint;
#endif

#endif /* SRC_CORE_XTK_TYPEDEFS_HPP_ */

