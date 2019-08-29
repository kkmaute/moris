/*
 * fn_Parsing_Tools.hpp
 *
 *  Created on: Aug 20, 2019
 *      Author: schmidt
 */

#ifndef MORIS_IOS_FN_PARSING_TOOLS_HPP_
#define MORIS_IOS_FN_PARSING_TOOLS_HPP_
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>

#include "typedefs.hpp"
#include "IO_Tools.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
// -----------------------------------------------------------------------------

template < typename T >
void string_to_mat( const std::string & aString,
                          Matrix< T > & aMat )
{
    if( aString.size() > 0 )
    {
        uint tCount = std::count( aString.begin(), aString.end(), ',') + 1;

        std::string tString( aString );

        // allocate memory
        aMat.set_size( tCount, 1 );

        // reset counter
        tCount = 0;

        // reset position
        size_t tPos = 0;

        // reset string
        tString = aString;

        while( tPos < tString.size() )
        {
            // find string
            tPos = tString.find( "," );

            // copy value into output matrix
            if( tPos <  tString.size() )
            {
                aMat( tCount++ ) = stod(  tString.substr( 0, tPos ) );
                tString =  tString.substr( tPos+1, tString.size() );
            }
        }

        // copy value into output matrix
        aMat( tCount++ ) = stod( tString );
    }
    else
    {
        aMat.set_size( 0, 1 );
    }
}

// -----------------------------------------------------------------------------

template < typename T >
void mat_to_string( const Matrix< T > & aMat,
                          std::string & aString )
{
    aString = "";

    uint tLength = aMat.length();

    for( uint k=0; k<tLength; ++k )
    {
        if( k > 0 )
        {
            aString = aString + ", ";
        }
        aString = aString + std::to_string( aMat( k ) );
    }
}

// -----------------------------------------------------------------------------

template < typename T >
void string_to_cell_mat( const std::string         & aString,
                               Cell< Matrix< T > > & aCellMat )
{
    if( aString.size() > 0 )
    {
        uint tCellCount = std::count( aString.begin(), aString.end(), ';') + 1;

        aCellMat.resize( tCellCount );

        std::string tString( aString );

        uint tCount = 0;

        // reset position
        size_t tPos = 0;

        bool tBool = true;

        while( tBool )
        {
            // find string
            tPos = tString.find( ";" );

            if( tPos == std::string::npos )
            {
                tBool = false;
            }

            // reset position
            size_t tPosSubString = 0;

            uint tCount1 = 0;

            if( tBool )
            {
                std::string tStringMat =  tString.substr( 0, tPos );

                uint tCountNum = std::count( tStringMat.begin(), tStringMat.end(), ',') + 1;

                aCellMat( tCount ).set_size( tCountNum, 1 );

                while( tPosSubString < tStringMat.size() )
                {
                    // find string
                    tPosSubString = tStringMat.find( "," );

                    // copy value into output matrix
                    if( tPosSubString < tStringMat.size() )
                    {
                        aCellMat( tCount )( tCount1++ ) = stod( tStringMat.substr( 0, tPosSubString ) );
                        tStringMat =  tStringMat.substr( tPosSubString+1, tStringMat.size() );
                    }
                }

                // copy value into output matrix
                aCellMat( tCount )( tCount1++ ) = stod( tStringMat );

                tString =  tString.substr( tPos+1, tString.size() );
            }
            else
            {
                uint tCountNum = std::count( tString.begin(), tString.end(), ',') + 1;

                aCellMat( tCount ).set_size( tCountNum, 1 );

                while( tPosSubString < tString.size() )
                {
                    // find string
                    tPosSubString = tString.find( "," );

                    // copy value into output matrix
                    if( tPosSubString < tString.size() )
                    {
                        aCellMat( tCount )( tCount1++ ) = stod( tString.substr( 0, tPosSubString ) );
                        tString =  tString.substr( tPosSubString+1, tString.size() );
                    }
                }

                // copy value into output matrix
                aCellMat( tCount )( tCount1++ ) = stod( tString );
            }
            tCount++;
        }
    }
    else
    {
        aCellMat( 0 );
    }
}

// -----------------------------------------------------------------------------

//    inline
//    bool string_to_bool( const std::string & aString )
//    {
//        // locale
//        std::locale loc;
//
//        // lower string of aString
//        std::string tLowerString( aString );
//        for( uint i=0; i < aString.length(); ++i)
//        {
//            tLowerString[ i ] = std::tolower( aString[i] );
//        }
//
//        return ( tLowerString == "true" || tLowerString == "on"
//              || tLowerString == "yes"  || tLowerString == "1" ) ;
//    }

// -----------------------------------------------------------------------------

}

#endif	/* MORIS_IOS_FN_PARSING_TOOLS_HPP_ */
