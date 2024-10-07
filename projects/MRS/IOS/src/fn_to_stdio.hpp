/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_to_stdio.hpp
 *
 */

#ifndef MORIS_IOS_FN_TO_STDIO_HPP_
#define MORIS_IOS_FN_TO_STDIO_HPP_

// C++ header files.
#include <ios>

// MORIS header files.
#include "cl_Logger.hpp"

//extern moris::Logger gLogger;

namespace moris::ios
{
    /**
     * @brief Convert an std::ios file access mode to stdio.
     *
     * +----------------------------------------------------------------------+
     * | std::ios        stdio    Effect                                      |
     * |----------------------------------------------------------------------+
     * | in              "r"      Open text file for reading only.            |
     * |----------------------------------------------------------------------+
     * | out|trunc       "w"      Truncate to 0 length, if existent,          |
     * | out             "w"      or create text file for writing only.       |
     * |----------------------------------------------------------------------+
     * | out|app         "a"      Append; open or create text file            |
     * |                          only for writing at end of file.            |
     * |----------------------------------------------------------------------+
     * | in|out          "r+"     Open text file for update                   |
     * |                          (reading and writing).                      |
     * |----------------------------------------------------------------------+
     * | in|out|trunc    "w+"     Truncate to 0 length, if existent,          |
     * |                          or create text file for update.             |
     * |----------------------------------------------------------------------+
     * | in|out|app      "a+"     Append; open or create text file            |
     * |                          for update, writing at end of file.         |
     * |----------------------------------------------------------------------+
     *
     * @param[in] mode std::ios file access mode.
     * @return stdio file access mode.
     */
    inline std::string
    to_stdio(
            std::ios::openmode const &mode )
    {
        if ( mode == ( std::ios::in ) )
            return "r";
        else if ( mode == ( std::ios::out | std::ios::trunc ) )
            return "w";
        else if ( mode == ( std::ios::app | std::ios::out ) )
            return "a";
        else if ( mode == ( std::ios::in | std::ios::out ) )
            return "r+";
        else if ( mode == ( std::ios::in | std::ios::out | std::ios::trunc ) )
            return "w+";
        else if ( mode == ( std::ios::app | std::ios::in | std::ios::out ) )
            return "a+";

        // I would prefer to throw an exception here,
        // but the standard only wants a NULL result.
        else
        {
            MORIS_LOG_ERROR( "File access mode is not defined." );
            return std::string();
        }
    }

}    // namespace moris::ios

#endif    /* MORIS_IOS_FN_TO_STDIO_HPP_ */

