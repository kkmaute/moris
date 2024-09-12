/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_Parsing_Tools.cpp
 *
 */

#include "fn_Parsing_Tools.hpp"

namespace moris
{
    // -----------------------------------------------------------------------------

    Vector< std::string > split_string(
            const std::string& aString,
            const std::string& aDelim )
    {
        // create empty cell of strings
        Vector< std::string > tCellOfStrings;

        size_t start;
        size_t end = 0;

        while ( ( start = aString.find_first_not_of( aDelim, end ) ) != std::string::npos )
        {
            end = aString.find( aDelim, start );
            tCellOfStrings.push_back( aString.substr( start, end - start ) );
        }

        return tCellOfStrings;
    }

    // -----------------------------------------------------------------------------

    void ltrim_string( std::string& aString )
    {
        auto it = std::find_if( aString.begin(), aString.end(), []( char c ) {
            return !std::isspace< char >( c, std::locale::classic() );
        } );
        aString.erase( aString.begin(), it );
    }

    // -----------------------------------------------------------------------------

    void rtrim_string( std::string& aString )
    {
        auto it = std::find_if( aString.rbegin(), aString.rend(), []( char c ) {
            return !std::isspace< char >( c, std::locale::classic() );
        } );
        aString.erase( it.base(), aString.end() );
    }

    // -----------------------------------------------------------------------------

    void trim_string( std::string& aString )
    {
        // if aString is empty return
        if ( aString.length() == 0 )
        {
            return;
        }

        // trim off trailing whitespaces
        rtrim_string( aString );

        // trim off leading whitespaces
        ltrim_string( aString );
    }

    // -----------------------------------------------------------------------------

    void split_trim_string(
            std::string&       aString,
            const std::string& aDelimiter )
    {
        char tDelim[] = " ";

        // if aString is empty return
        if ( aString.empty() )
        {
            return;
        }
        // if no delimiter is provided just trim the given string
        if ( aDelimiter.empty() )
        {
            // trim substring
            trim_string( aString );
            return;
        }

        // loop over all delimiter characters
        for ( uint id = 0; id < aDelimiter.length(); id++ )
        {
            // check if delimiter is first element of string
            bool isFirst = aString.front() == aDelimiter[ id ];

            // check if delimiter is last element of string
            bool isLast = aString.back() == aDelimiter[ id ];

            // get char* of current delimiter
            strncpy( tDelim, &aDelimiter[ id ], 1 );

            // allocate vector for storing substrings
            std::vector< std::string > tSubStringVec;

            // get first substring
            char* token = strtok( const_cast< char* >( aString.c_str() ), tDelim );

            // loop over all substrings
            while ( token != nullptr )
            {
                // create substring
                std::string tSubstring( token );

                // trim substring
                trim_string( tSubstring );

                // check if substring is empty or has only white spaces
                if ( !tSubstring.empty() )
                {
                    if ( !std::all_of( tSubstring.begin(), tSubstring.end(), isspace ) )
                    {
                        tSubStringVec.push_back( tSubstring );
                    }
                }

                // get next substring substring
                token = strtok( nullptr, tDelim );
            }

            // reassemble substrings
            aString.clear();

            if ( isFirst ) aString = tDelim;

            if ( tSubStringVec.size() > 0 )
            {
                aString += tSubStringVec[ 0 ];

                for ( uint is = 1; is < tSubStringVec.size(); ++is )
                {
                    aString += tDelim;
                    aString += tSubStringVec[ is ];
                }
            }
            if ( isLast ) aString += tDelim;
        }
    }

    // -----------------------------------------------------------------------------

    template<>
    void string_to_cell< std::string >(
            const std::string&     aString,
            Vector< std::string >& aCell )
    {
        // check that vector is empty
        MORIS_ASSERT( aCell.size() == 0, "string_to_cell - vector needs to be empty" );

        // convert string to string stream and a sub string
        std::stringstream tStringStream( aString );

        // separate the string by the delimiter
        std::string tSubString;
        while ( std::getline( tStringStream, tSubString, ',' ) )
        {
            // Add to the cell
            aCell.push_back( tSubString );
        }
    }

    // -----------------------------------------------------------------------------
}    // namespace moris
