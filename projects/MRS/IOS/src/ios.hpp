/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ios.hpp
 *
 */

#ifndef MORIS_IOS_HPP_
#define MORIS_IOS_HPP_

#include "core.hpp"

/**
 * IOS macro constants.
 * These macro constants expands to integral expressions
 * corresponding to the size needed for an array of char elements
 * to hold the longest file name or string name allowed by the library.
 *
 * Modeled after FILENAME_MAX in C++.
 * @see [FILENAME_MAX](http://www.cplusplus.com/reference/cstdio/FILENAME_MAX)
 */
const moris::size_t MORIS_FILENAME_MAX = 255;
const moris::size_t MORIS_STRING_MAX   = 31;
const moris::size_t MORIS_FIELD_MAX    = 63;
const moris::size_t MORIS_LINE_MAX     = 255;

#include "cl_Logger.hpp"
#include "fn_Logger.hpp"
#include "fn_to_stdio.hpp"

#endif /* MORIS_IOS_HPP_ */

