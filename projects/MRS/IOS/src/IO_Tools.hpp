/*
 * IO_Tools.hpp
 *
 *  Created on: Jan 7, 2019
 *      Author: schmidt
 */
#ifndef PROJECTS_MRS_IOS_SRC_IO_TOOLS_HPP_
#define PROJECTS_MRS_IOS_SRC_IO_TOOLS_HPP_

#include <cstdio>
#include <string>
#include <fstream>
#include <memory>

template< typename ... Args >
std::string print_log( const Args ... aArgs )
{
    // Determine size of string
    auto tSize = snprintf( nullptr, 0, aArgs ... );

    // create char pointer with size of string length + 1 for \0
    std::unique_ptr< char[] > tMsg( new char[ tSize + 1 ]);

    // write string into buffered char pointer
    snprintf( tMsg.get(), tSize + 1, aArgs ... );

    return( std::string( tMsg.get(), tMsg.get() + tSize ).c_str() );
}

#endif /* PROJECTS_MRS_IOS_SRC_IO_TOOLS_HPP_ */
