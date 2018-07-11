#ifndef MORIS_ASSERT_FN_ASSERT_HPP_
#define MORIS_ASSERT_FN_ASSERT_HPP_

// MORIS header files.
#include "fn_log_assert.hpp" // ASR/src
#include "exceptions.hpp"

/**
 * @brief Assertion, only active if DEBUG is defined.
 */
#if !defined(NDEBUG) || defined(DEBUG)
#define MORIS_ASSERT(check, msg) \
    do \
    { \
        if (! (check)) \
        { \
            moris::assert::moris_assert(__FILE__, __LINE__, __FUNCTION__, #check, msg); \
        } \
    } while (false)
#else
#define MORIS_ASSERT(check, msg)
#endif

/**
 * @brief Exit, active in debug and opt mode.
 */
#if !defined(NDEBUG) || defined(DEBUG)
#define MORIS_ERROR(check, msg) \
    do \
    { \
        if (! (check)) \
        { \
            moris::assert::moris_assert(__FILE__, __LINE__, __FUNCTION__, #check, msg); \
        } \
    } while (false)
#else
#define MORIS_ERROR(check, msg) \
    do \
    { \
        if (! (check)) \
        { \
            moris::assert::moris_assert(__FILE__, __LINE__, __FUNCTION__, #check, msg); \
        } \
    } while (false)
#endif

#endif /* MORIS_ASSERT_FN_ASSERT_HPP_ */
