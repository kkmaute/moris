#ifndef MORIS_IOS_CL_LOGGER_SINK_HPP_
#define MORIS_IOS_CL_LOGGER_SINK_HPP_

// C++ header files.
#include <string>

// Third-party header files.
#include <iosfwd>
#include <boost/iostreams/categories.hpp>

// MORIS header files.
#include "core.hpp"
#include "cl_Logger.hpp"

namespace moris
{
namespace ios
{
    /**
     * Logger sink.
     */
    class Logger_Sink
    {
    public:

        /**
         * Type of characters handled by the device.
         * Required by Boost Iostreams library.
         */
        typedef char char_type;

        /**
         * Fundamental I/O operations supported by the device.
         * Required by Boost Iostreams library.
         */
        typedef boost::iostreams::sink_tag category;

        /**
         * moris::ios::Logger_Sink default constructor.
         */
        Logger_Sink() = default;

        /**
         * moris::ios::Logger_Sink default destructor.
         */
        virtual
        ~Logger_Sink() = default;

        /**
         * Write up to n characters from the buffer s into the moris logger,
         * returning the number of characters written.
         *
         * @param[in] s What to write.
         * @param[in] n How much to write.
         * @return Number of characters written.
         */
        virtual
        std::streamsize
        write(
                char_type const * s,
                std::streamsize   n) = 0;
    };

    /**
     * Logger sink for a trace severity message.
     *
     * @include "cl_Logger_Sink/clog.inc"
     */
    class clog_Sink: public moris::ios::Logger_Sink
    {
    public:

        /**
         * Write a trace message to the logger.
         */
        std::streamsize
        write(
                char_type const * s,
                std::streamsize   n)
        {
            MORIS_LOG_TRACE << std::string(s, n);
            return n;
        }
    };

    /**
     * Logger sink for an informational severity message.
     *
     * @include "cl_Logger_Sink/cout.inc"
     */
    class cout_Sink: public moris::ios::Logger_Sink
    {
    public:

        /**
         * Write an information message to the logger.
         */
        std::streamsize
        write(
                char_type const * s,
                std::streamsize   n)
        {
            MORIS_LOG_INFO << std::string(s, n);
            return n;
        }
    };

    /**
     * Logger sink for an error severity message.
     *
     * @include "cl_Logger_Sink/cerr.inc"
     */
    class cerr_Sink: public moris::ios::Logger_Sink
    {
    public:

        /**
         * Write an error message to the logger.
         *
         */
        std::streamsize
        write(
                char_type const * s,
                std::streamsize   n)
        {
            MORIS_LOG_ERROR << std::string(s, n);
            return n;
        }
    };

}    // namespace ios
}    // namespace moris

#endif    /* MORIS_IOS_CL_LOGGER_SINK_HPP_ */
