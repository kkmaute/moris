#ifndef MORIS_IOS_FN_CLOG_HPP_
#define MORIS_IOS_FN_CLOG_HPP_

// MORIS header files.
#include "cl_Logger_Sink.hpp"

namespace moris
{
    /**
     * Object of class std::ostream that represents
     * the trace logging stream oriented to narrow characters
     * (of type moris::char_t).
     *
     * @include fn_iostream/clog.inc
     */
    extern std::ostream clog;

    /**
     * Object of class std::ostream that represents
     * the informational output stream oriented to narrow characters
     * (of type moris::char_t).
     *
     * @include "fn_iostream/cout.inc"
     */
    extern std::ostream cout;

    /**
     * Object of class std::ostream that represents
     * the error stream oriented to narrow characters
     * (of type moris::char_t).
     *
     * @include "fn_iostream/cerr.inc"
     */
    extern std::ostream cerr;

}    // namespace moris

#endif    /* MORIS_IOS_FN_CLOG_HPP_ */
