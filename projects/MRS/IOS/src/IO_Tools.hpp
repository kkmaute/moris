/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * IO_Tools.hpp
 *
 */

#ifndef PROJECTS_MRS_IOS_SRC_IO_TOOLS_HPP_
#define PROJECTS_MRS_IOS_SRC_IO_TOOLS_HPP_

#include <cstdio>
#include <string>
#include <fstream>
#include <memory>

template< typename... Args >
std::string
print_log( const char* aFormat, const Args... aArgs )
{
    // consider the case that print argument is not empty
    if constexpr ( sizeof...( Args ) != 0 )
    {
        // Determine size of string
        auto tSize = snprintf( nullptr, 0, aFormat, aArgs... );

        // create char pointer with size of string length + 1 for \0
        std::unique_ptr< char[] > tMsg( new char[ tSize + 1 ] );

        // write string into buffered char pointer
        snprintf( tMsg.get(), tSize + 1, aFormat, aArgs... );

        return ( std::string( tMsg.get(), tMsg.get() + tSize ).c_str() );
    }

    return aFormat;
}

#endif /* PROJECTS_MRS_IOS_SRC_IO_TOOLS_HPP_ */
