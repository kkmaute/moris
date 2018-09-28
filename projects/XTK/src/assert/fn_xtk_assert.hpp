/*
 * fn_assert.hpp
 *
 *  Created on: Jun 19, 2017
 *      Author: ktdoble
 */

#ifndef SRC_ASSERT_FN_ASSERT_HPP_
#define SRC_ASSERT_FN_ASSERT_HPP_

// XTK Header files
#include "fn_xtk_log_assert.hpp"

/**
 * @brief Assertion, only active if DEBUG is defined.
 */
#if !defined(NDEBUG) || defined(DEBUG)
#define XTK_ASSERT(check, msg) \
    do \
    { \
        if (! (check)) \
        { \
            xtk::assert::xtk_assert(__FILE__, __LINE__, __FUNCTION__, #check, msg); \
            MPI_Abort(MPI_COMM_WORLD,0);\
        } \
    } while (false)
#else
#define XTK_ASSERT(check, msg) \
    do \
    { \
        if (! (check)) \
        { \
            xtk::assert::xtk_assert(__FILE__, __LINE__, __FUNCTION__, #check, msg); \
            MPI_Abort(MPI_COMM_WORLD,0);\
        } \
    } while (false)
#endif


#endif /* SRC_ASSERT_FN_ASSERT_HPP_ */
