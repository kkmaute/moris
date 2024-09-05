/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Waitbar.hpp
 *
 */

#ifndef MORIS_CHRONOS_CL_WAITBAR_HPP_
#define MORIS_CHRONOS_CL_WAITBAR_HPP_

// Third-party header files.
// Boost Timer library requires this fix
// to use auto_cpu_timer and progress_display together.
// @see [How to use boost auto_cpu_timer and progress_display at the same time?]
// (http://stackoverflow.com/questions/12760828/how-to-use-boost-auto-cpu-timer-and-progress-display-at-the-same-time)
#define timer boost_timer
#include <boost/timer/progress_display.hpp>

namespace moris
{
    namespace chronos
    {
        /**
         * Waitbar API wrapper for Boost progress_display.
         * A wait bar is a figure that displays what percentage of a calculation
         * is complete as the calculation proceeds
         * by progressively filling a bar from left to right.
         */
        using Waitbar = boost::timer::progress_display;

    }    // namespace chronos
}    // namespace moris

#endif /* MORIS_CHRONOS_CL_WAITBAR_HPP_ */
