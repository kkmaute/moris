/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Stopwatch.hpp
 *
 */

#ifndef MORIS_CHRONOS_CL_STOPWATCH_HPP_
#define MORIS_CHRONOS_CL_STOPWATCH_HPP_

// C++ header files.
#include <chrono>

// Third-party header files.
#include <boost/timer/timer.hpp>

// MORIS header files.
#include "typedefs.hpp" // COR/src

namespace moris
{
    namespace chronos
    {
        // Boost cpu_times API.
        using cpu_times = boost::timer::cpu_times;

        // STL time units API.
        using nanoseconds  = std::chrono::nanoseconds;
        using microseconds = std::chrono::microseconds;
        using milliseconds = std::chrono::milliseconds;
        using seconds      = std::chrono::seconds;
        using minutes      = std::chrono::minutes;
        using hours        = std::chrono::hours;

        /**
         * Retrieving elapsed and wall-clock time is problematic
         * due to platform-dependent and language-dependent issues.
         * To avoid those problems, MORIS furnishes the moris::chronos::Stopwatch class.
         *
         * It provides the typical start, stop, resume and reset features of a stopwatch.
         *
         * For further documentation, visit <a href ='http://www.boost.org/doc/libs/1_48_0/libs/timer/doc/cpu_timers.html' >this</a>
         */
        class Stopwatch
        {
        public:

            /**
             * moris::chronos::Stopwatch default destructor.
             */
            Stopwatch() noexcept
            : timer_m()
              {
              }

            /**
             * moris::chronos::Stopwatch destructor.
             */
            ~Stopwatch() = default;

        private:

            /**
             * Boost CPU timer.
             *
             * @see [Boost Timer Library]
             * (http://www.boost.org/doc/libs/1_58_0/libs/timer/doc/index.html)
             */
            boost::timer::cpu_timer timer_m;

        public:

            /**
             * @brief Reset start time.
             *
             * Resets the start time for the moris::chronos::Stopwatch
             * to the current time.
             * A code section can be timed by putting it
             * between a call to reset() and toc().
             */
            void
            reset() noexcept
            {
                timer_m.start();
            }

            /**
             * If !is_stopped(), stops accumulating elapsed time.
             */
            void
            stop() noexcept
            {
                timer_m.stop();
            }

            /**
             * Is the moris::chronos::Stopwatch stopped?
             *
             * @return true if stop() was the most recent action function called,
             *         otherwise false.
             */
            bool
            is_stopped() const noexcept
            {
                return timer_m.is_stopped();
            }

            /**
             * @brief Restarts the timer, accumulating additional elapsed time.
             *
             * If is_stopped(), resumes accumulating additional elapsed time,
             * as of the current time values. Otherwise, no effect.
             */
            void
            resume() noexcept
            {
                timer_m.resume();
            }

            /**
             * @brief Elapsed time function.
             *
             * Returns the elapsed time in seconds
             * since the moris::chronos::Stopwatch object was constructed,
             * or since the toc() function was called.
             *
             * A code section can be timed by putting it
             * between the moris::chronos::Stopwatch constructor
             * and a call to toc(),
             * or between a call to reset() and toc().
             *
             * @return Elapsed time.
             */
            template<typename T = moris::chronos::nanoseconds>
            moris::chronos::cpu_times
            toc() const noexcept
            {
                auto elapsed_time = timer_m.elapsed();

                moris::chronos::nanoseconds wall_time(elapsed_time.wall);
                moris::chronos::nanoseconds user_time(elapsed_time.user);
                moris::chronos::nanoseconds system_time(elapsed_time.system);

                return moris::chronos::cpu_times{
                        std::chrono::duration_cast<T>(wall_time).count(),
                        std::chrono::duration_cast<T>(user_time).count(),
                        std::chrono::duration_cast<T>(system_time).count()};
            }

            /**
             * @brief Same as toc(), but returns formatted timing information.
             *
             * @todo[chvillanuevap@gmail.com]
             * The Boost cpu_timer only prints timing information in seconds.
             * Add capability to display in other units of time.
             *
             * @param[in] format_a Timing display format.
             * @return Formatted timing info.
             */
            template<typename T = moris::chronos::nanoseconds>
            std::string
            elapsed(
                    std::string const & format_a =
                            "%ws wall, %us user + %ss system = %ts CPU (%p%)") const
            {
                auto elapsed_time = this->toc<T>();

                return ::boost::timer::format(elapsed_time, 9, format_a);
            }
        };

    }// namespace chronos

    // Alias for a Matlab-like tic-toc functionality.
    using tic = moris::chronos::Stopwatch;

}// namespace moris

#endif/* MORIS_CHRONOS_CL_STOPWATCH_HPP_ */

