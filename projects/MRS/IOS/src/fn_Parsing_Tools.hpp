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
//#include <map>

#include "typedefs.hpp"
#include "IO_Tools.hpp"
#include "cl_Matrix.hpp"
#include "cl_Map.hpp"

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

template < typename T >
void string_to_cell_mat_2( const std::string         & aString,
                                 Cell< Matrix< T > > & aCellMat )
{
    // if non-empty string
    if( aString.size() > 0 )
    {
        // get number of cells
        uint tCellCount = std::count( aString.begin(), aString.end(), '/') + 1;

        // set size for the list
        aCellMat.resize( tCellCount );

        // get the string
        std::string tString( aString );

        // init cell count
        uint tCount = 0;

        // set global reading position
        size_t tPos = 0;

        // not last cell
        bool tBool = true;

        // while not last cell
        while( tBool )
        {
            // find end of first matrix
            tPos = tString.find( "/" );

            // if end of string
            if( tPos == std::string::npos )
            {
                // last cell
                tBool = false;
            }

            // set local reading position
            size_t tPosSubString = 0;

            // set row counter
            uint tCount1 = 0;

            // if not last cell
            if( tBool )
            {
                // get string for matrix
                std::string tStringMat =  tString.substr( 0, tPos );

                // get number of row
                uint tCountRow = std::count( tStringMat.begin(), tStringMat.end(), ';') + 1;

                // get first row
                std::string tStringFirstRow = tString.substr( 0, tStringMat.find( ";" ) );

                // get number of cols
                uint tCountCol = std::count( tStringFirstRow.begin(), tStringFirstRow.end(), ',') + 1;

                // set matrix size
                aCellMat( tCount ).set_size( tCountRow, tCountCol );

                // not last row
                bool tBoolMat = true;

                while( tBoolMat )
                {
                    // find end of first row
                    tPosSubString = tStringMat.find( ";" );

                    // if end of mat string
                    if( tPosSubString == std::string::npos )
                    {
                        // last row
                        tBoolMat = false;
                    }

                    // set local reading position
                    size_t tPosSubSubString = 0;

                   // set col counter
                   uint tCount2 = 0;

                   if( tBoolMat )
                   {
                       // get string for row
                       std::string tStringRow =  tStringMat.substr( 0, tPosSubString );

                        while( tPosSubSubString < tStringRow.size() )
                        {
                            // find string
                            tPosSubSubString = tStringRow.find( "," );

                            // copy value into output matrix
                            if( tPosSubSubString < tStringRow.size() )
                            {
                                aCellMat( tCount )( tCount1, tCount2++ ) = stod( tStringRow.substr( 0, tPosSubSubString ) );
                                tStringRow =  tStringRow.substr( tPosSubSubString+1, tStringRow.size() );
                            }
                        }
                        // copy value into output matrix
                        aCellMat( tCount )( tCount1, tCount2++ ) = stod( tStringRow );
                        tStringMat =  tStringMat.substr( tPosSubString+1, tStringMat.size() );
                    }
                   else
                   {
                       while( tPosSubSubString < tStringMat.size() )
                       {
                           // find string
                           tPosSubSubString = tStringMat.find( "," );

                           // copy value into output matrix
                           if( tPosSubSubString < tStringMat.size() )
                           {
                               aCellMat( tCount )( tCount1, tCount2++ ) = stod( tStringMat.substr( 0, tPosSubSubString ) );
                               tStringMat =  tStringMat.substr( tPosSubSubString+1, tStringMat.size() );
                           }
                       }
                       // copy value into output matrix
                       aCellMat( tCount )( tCount1, tCount2++ ) = stod( tStringMat );
                       tStringMat = tStringMat.substr( tPosSubString+1, tStringMat.size() );
                    }
                   tCount1++;
                }
                tString =  tString.substr( tPos+1, tString.size() );
            }
            else
            {
                // get number of row in last matrix
                uint tCountRow = std::count( tString.begin(), tString.end(), ';') + 1;

                // get first row
                std::string tStringFirstRow = tString.substr( 0, tString.find( ";" ) );

                // get number of cols
                uint tCountCol = std::count( tStringFirstRow.begin(), tStringFirstRow.end(), ',') + 1;

                // set matrix size
                aCellMat( tCount ).set_size( tCountRow, tCountCol );

                // not last row
                bool tBoolMat = true;

                while( tBoolMat )
                {
                    // find end of first row
                    tPosSubString = tString.find( ";" );

                    // if end of mat string
                    if( tPosSubString == std::string::npos )
                    {
                        // last row
                        tBoolMat = false;
                    }

                    // set local reading position
                    size_t tPosSubSubString = 0;

                   // set col counter
                   uint tCount2 = 0;

                   if( tBoolMat )
                   {
                       // get string for row
                       std::string tStringRow =  tString.substr( 0, tPosSubString );

                        while( tPosSubSubString < tStringRow.size() )
                        {
                            // find string
                            tPosSubSubString = tStringRow.find( "," );

                            // copy value into output matrix
                            if( tPosSubSubString < tStringRow.size() )
                            {
                                aCellMat( tCount )( tCount1, tCount2++ ) = stod( tStringRow.substr( 0, tPosSubSubString ) );
                                tStringRow =  tStringRow.substr( tPosSubSubString+1, tStringRow.size() );
                            }
                        }
                        // copy value into output matrix
                        aCellMat( tCount )( tCount1, tCount2++ ) = stod( tStringRow );
                        tString =  tString.substr( tPosSubString+1, tString.size() );
                    }
                   else
                   {
                       while( tPosSubSubString < tString.size() )
                       {
                           // find string
                           tPosSubSubString = tString.find( "," );

                           // copy value into output matrix
                           if( tPosSubSubString < tString.size() )
                           {
                               aCellMat( tCount )( tCount1, tCount2++ ) = stod( tString.substr( 0, tPosSubSubString ) );
                               tString =  tString.substr( tPosSubSubString+1, tString.size() );
                           }
                       }
                       // copy value into output matrix
                       aCellMat( tCount )( tCount1, tCount2++ ) = stod( tString );
                    }
                   tCount1++;
                }
            }
            tCount++;
        }
    }
//    // if empty string
//    else
//    {
//        aCellMat( 0 );
//    }
}

// -----------------------------------------------------------------------------

template < typename T >
void string_to_cell_of_cell( const std::string                     & aString,
                                   moris::Cell< moris::Cell< T > > & aCellCell,
                                   moris::map< std::string, T >    & aMap )
{
    if( aString.size() > 0 )
    {
        uint tCellCount = std::count( aString.begin(), aString.end(), ';') + 1;

        aCellCell.resize( tCellCount );

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

                aCellCell( tCount ).resize( tCountNum );

                while( tPosSubString < tStringMat.size() )
                {
                    // find string
                    tPosSubString = tStringMat.find( "," );

                    // copy value into output matrix
                    if( tPosSubString < tStringMat.size() )
                    {
                        T tComponent = aMap.find( tStringMat.substr( 0, tPosSubString ) );
                        aCellCell( tCount )( tCount1++ ) = tComponent;
                        tStringMat =  tStringMat.substr( tPosSubString+1, tStringMat.size() );
                    }
                }

                // copy value into output matrix
                T tComponent = aMap.find( tStringMat );
                aCellCell( tCount )( tCount1++ ) = tComponent;
                tString =  tString.substr( tPos+1, tString.size() );
            }
            else
            {
                uint tCountNum = std::count( tString.begin(), tString.end(), ',') + 1;

                aCellCell( tCount ).resize( tCountNum );

                while( tPosSubString < tString.size() )
                {
                    // find string
                    tPosSubString = tString.find( "," );

                    // copy value into output matrix
                    if( tPosSubString < tString.size() )
                    {
                        T tComponent = aMap.find( tString.substr( 0, tPosSubString ) );
                        aCellCell( tCount )( tCount1++ ) = tComponent;
                        tString =  tString.substr( tPosSubString+1, tString.size() );
                    }
                }

                // copy value into output matrix
                T tComponent = aMap.find( tString );
                aCellCell( tCount )( tCount1++ ) = tComponent;
            }
            tCount++;
        }
    }
    else
    {
        aCellCell.resize( 0 );
    }
}

template < typename T >
void string_to_cell_of_cell( const std::string                     & aString,
                                   moris::Cell< moris::Cell< T > > & aCellCell )
{
    if( aString.size() > 0 )
    {
        uint tCellCount = std::count( aString.begin(), aString.end(), ';') + 1;

        aCellCell.resize( tCellCount );

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

                aCellCell( tCount ).resize( tCountNum );

                while( tPosSubString < tStringMat.size() )
                {
                    // find string
                    tPosSubString = tStringMat.find( "," );

                    // copy value into output matrix
                    if( tPosSubString < tStringMat.size() )
                    {
                        aCellCell( tCount )( tCount1++ ) = tStringMat.substr( 0, tPosSubString );
                        tStringMat =  tStringMat.substr( tPosSubString+1, tStringMat.size() );
                    }
                }

                // copy value into output matrix
                aCellCell( tCount )( tCount1++ ) = tStringMat;
                tString =  tString.substr( tPos+1, tString.size() );
            }
            else
            {
                uint tCountNum = std::count( tString.begin(), tString.end(), ',') + 1;

                aCellCell( tCount ).resize( tCountNum );

                while( tPosSubString < tString.size() )
                {
                    // find string
                    tPosSubString = tString.find( "," );

                    // copy value into output matrix
                    if( tPosSubString < tString.size() )
                    {
                        aCellCell( tCount )( tCount1++ ) = tString.substr( 0, tPosSubString );
                        tString =  tString.substr( tPosSubString+1, tString.size() );
                    }
                }

                // copy value into output matrix
                aCellCell( tCount )( tCount1++ ) = tString;
            }
            tCount++;
        }
    }
    else
    {
        aCellCell.resize( 0 );
    }
}

template < typename T >
void string_to_cell( const std::string & aString,
                     moris::Cell< T >  & aCell,
                     moris::map< std::string, T >    & aMap )
{
    if( aString.size() > 0 )
    {
        uint tCellCount = std::count( aString.begin(), aString.end(), ',') + 1;

        aCell.resize( tCellCount );

        std::string tString( aString );

        // reset position
        size_t tPos = 0;

        uint tCount = 0;

        bool tBool = true;

        while( tBool )
        {
            // find string
            tPos = tString.find( "," );

            if( tPos == std::string::npos )
            {
                tBool = false;
            }

            if( tBool )
            {
                std::string tStringMat = tString.substr( 0, tPos );

                // copy value into output matrix
                T tComponent = aMap.find( tStringMat );
                aCell( tCount++ ) = tComponent;
                tString =  tString.substr( tPos+1, tString.size() );
            }
            else
            {
                // copy value into output matrix
                T tComponent = aMap.find( tString );
                aCell( tCount++ ) = tComponent;
            }
        }
    }
    else
    {
        aCell.resize( 0 );
    }
}

template < typename T >
void string_to_cell( const std::string & aString,
                     moris::Cell< T >  & aCell )
{
    if( aString.size() > 0 )
    {
        uint tCellCount = std::count( aString.begin(), aString.end(), ',') + 1;

        aCell.resize( tCellCount );

        std::string tString( aString );

        // reset position
        size_t tPos = 0;

        uint tCount = 0;

        bool tBool = true;

        while( tBool )
        {
            // find string
            tPos = tString.find( "," );

            if( tPos == std::string::npos )
            {
                tBool = false;
            }

            if( tBool )
            {
                std::string tStringMat = tString.substr( 0, tPos );
                // copy value into output matrix
                aCell( tCount++ ) = tStringMat;
                tString =  tString.substr( tPos+1, tString.size() );
            }
            else
            {
                // copy value into output matrix
                aCell( tCount++ ) = tString;
            }
        }
    }
//    else
//    {
//        aCell.resize( 0 );
//    }
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
