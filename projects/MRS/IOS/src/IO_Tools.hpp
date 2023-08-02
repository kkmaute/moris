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
#include <cstdarg>
#include <string>
#include <fstream>
#include <memory>

#ifdef __GNUC__
__attribute__ ((format (printf, 1, 2)))
#endif
inline
std::string
print_log( const char* aFormat, ... )
{
    va_list args, args2;
    va_start(args, aFormat);

    // Determine size of string
    va_copy(args2, args);
    auto tSize = vsnprintf( nullptr, 0, aFormat, args2 );
    va_end(args2);

    // create char pointer with size of string length + 1 for \0
    std::unique_ptr< char[] > tMsg( new char[ tSize + 1 ] );

    // write string into buffered char pointer
    vsnprintf( tMsg.get(), tSize + 1, aFormat, args );
    va_end(args);

    return ( std::string( tMsg.get(), tMsg.get() + tSize ).c_str() );
}

#endif /* PROJECTS_MRS_IOS_SRC_IO_TOOLS_HPP_ */
