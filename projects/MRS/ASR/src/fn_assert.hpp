/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_assert.hpp
 *
 */

#ifndef MORIS_ASSERT_FN_ASSERT_HPP_
#define MORIS_ASSERT_FN_ASSERT_HPP_

// MORIS header files.
#include "fn_log_assert.hpp"    // ASR/src
#include "exceptions.hpp"

/**
 * @brief Assertion, only active if DEBUG is defined.
 */
#if defined( MORIS_HAVE_DEBUG )
#define MORIS_ASSERT( check, ... )                                                                \
    do                                                                                            \
    {                                                                                             \
        if ( !( check ) )                                                                         \
        {                                                                                         \
            moris::assert::moris_assert( __FILE__, __LINE__, __FUNCTION__, #check, __VA_ARGS__ ); \
        }                                                                                         \
    } while ( false )
#else
#define MORIS_ASSERT( check, ... )
#endif

/**
 * @brief Assertion for memory usage, only active if CHECK_MEMORY is defined.
 */
#if defined( CHECK_MEMORY )
#define MORIS_CHECK_MEMORY( check, ... )                                                           \
    do                                                                                             \
    {                                                                                              \
        if ( !( check ) )                                                                          \
        {                                                                                          \
            moris::assert::moris_warning( __FILE__, __LINE__, __FUNCTION__, #check, __VA_ARGS__ ); \
        }                                                                                          \
    } while ( false )
#else
#define MORIS_CHECK_MEMORY( check, ... )
#endif

/**
 * @brief Exit, active in debug and opt mode.
 */
#define MORIS_ERROR( check, ... )                                                                 \
    do                                                                                            \
    {                                                                                             \
        if ( !( check ) )                                                                         \
        {                                                                                         \
            moris::assert::moris_assert( __FILE__, __LINE__, __FUNCTION__, #check, __VA_ARGS__ ); \
        }                                                                                         \
    } while ( false )

#endif /* MORIS_ASSERT_FN_ASSERT_HPP_ */
