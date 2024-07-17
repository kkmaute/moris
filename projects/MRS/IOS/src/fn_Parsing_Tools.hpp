/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_Parsing_Tools.hpp
 *
 */

#pragma once

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>

#include "moris_typedefs.hpp"
#include "IO_Tools.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "cl_Map.hpp"

namespace moris
{
    // -----------------------------------------------------------------------------
    /**
     * splits string into substrings based on delimiter
     *
     * @param[in] aString  input string
     * @param[in] aDelim   delimiter
     *
     * @return  cell of strings
     */
    Vector< std::string >
    split_string(
            const std::string &aString,
            const std::string &aDelim );

    // -----------------------------------------------------------------------------
    /**
     * removes leading whitespace of string
     *
     * @param s std::string
     */
    void ltrim_string( std::string &aString );

    // -----------------------------------------------------------------------------
    /**
     * removes trailing whitespace of string
     *
     * @param aString std::string
     */
    void rtrim_string( std::string &aString );

    // -----------------------------------------------------------------------------
    /**
     * removes leading and trailing whitespace of string
     *
     * @param aString std::string
     */
    void trim_string( std::string &aString );

    // -----------------------------------------------------------------------------
    /**
     * splits a string in sub-strings defined by a set of delimiters, removes leading and trailing whitespace
     *  of sub-string, and reassembles sub-strings into single string; if not delimiter is provides input string
     *  will be trimmed.
     *
     * @param aString std::string
     * @param aDelimiter const std::string
     */
    void split_trim_string(
            std::string       &aString,
            const std::string &aDelimiter );

    // -----------------------------------------------------------------------------
    /**
     * splits a string in sub-strings based on a delimiter and then converts substrings
     * to numerical values which are stored in a vector
     *
     * @param[in] aString    const std::string
     * @param[in] aDelimiter const std::string (size = 1)
     *
     * @param[out] aMat      matrix
     */
    template< typename T >
    void
    delimited_string_to_mat(
            const std::string &aString,
            const std::string &aDelim,
            Matrix< T >       &aMat )
    {
        // check for empty string
        if ( aString.empty() )
        {
            aMat.set_size( 0, 0 );
            return;
        }

        // split string into substrings
        Vector< std::string > tSubStringVec = split_string( aString, aDelim );

        // get number of substrings
        uint tNumberOfSubStrings = tSubStringVec.size();

        if ( tNumberOfSubStrings == 0 )
        {
            aMat.set_size( 0, 0 );
            return;
        }

        // set size of matrix
        aMat.set_size( tNumberOfSubStrings, 1 );

        // convert strings to numerical value
        for ( uint tStrgIndex = 0; tStrgIndex < tNumberOfSubStrings; ++tStrgIndex )
        {
            aMat( tStrgIndex ) = stod( tSubStringVec( tStrgIndex ) );
        }
    }

    // -----------------------------------------------------------------------------
    /**
     * Converts a given string into a matrix of templated type;
     * rows are separated by ";", columns are separated by ","
     *
     * @tparam T Type of moris::Matrix
     * @param aString String to be converted into a matrix
     * @param aMat Matrix converted from string, returned as reference
     */

    template< typename T >
    void
    string_to_mat(
            const std::string &aString,
            Matrix< T >       &aMat )
    {
        if ( !aString.empty() )
        {
            uint tCountRow = std::count( aString.begin(), aString.end(), ';' ) + 1;
            uint tCountCol = ( std::count( aString.begin(), aString.end(), ',' ) / tCountRow ) + 1;

            std::string tString( aString );

            // allocate memory
            aMat.set_size( tCountRow, tCountCol );

            // reset position
            size_t tPos;

            // reset string
            tString = aString;

            // Create output matrix
            uint tRowIndex;
            uint tColIndex;
            for ( tRowIndex = 0; tRowIndex < tCountRow; tRowIndex++ )
            {
                for ( tColIndex = 0; tColIndex < tCountCol - 1; tColIndex++ )
                {
                    // find string
                    tPos = tString.find( ',' );

                    // copy value into output matrix
                    if ( tPos < tString.size() )
                    {
                        aMat( tRowIndex, tColIndex ) = stod( tString.substr( 0, tPos ) );
                        tString                      = tString.substr( tPos + 1, tString.size() );
                    }
                }

                // Last value before next row
                tPos = tString.find( ';' );

                // copy value into output matrix
                if ( tPos < tString.size() )
                {
                    aMat( tRowIndex, tColIndex ) = stod( tString.substr( 0, tPos ) );
                    tString                      = tString.substr( tPos + 1, tString.size() );
                }
            }

            // copy value into output matrix
            aMat( tCountRow - 1, tCountCol - 1 ) = stod( tString );
        }
        else
        {
            aMat.set_size( 0, 1 );
        }
    }

    // -----------------------------------------------------------------------------
    /**
     * Converts a given string into a matrix of templated type, which is returned instead of passed by reference
     * rows are separated by ";", columns are separated by ","
     *
     * @tparam T Type of moris::Matrix
     * @param aString String to be converted into a matrix
     * @return Matrix converted from string
     */

    template< typename T >
    Matrix< T >
    string_to_mat( const std::string &aString )
    {
        Matrix< T > tMat;
        string_to_mat( aString, tMat );

        return tMat;
    }

    // -----------------------------------------------------------------------------
    /**
     * Converts a matrix into a string with the matrix components separated by ","
     *
     * @tparam T Type of moris::Matrix
     * @param aMat Matrix to be converted into a string
     * @return String converted from matrix
     */

    template< typename T >
    void
    mat_to_string(
            const Matrix< T > &aMat,
            std::string       &aString )
    {
        aString = "";

        uint tLength = aMat.length();

        for ( uint k = 0; k < tLength; ++k )
        {
            if ( k > 0 )
            {
                aString += ", ";
            }
            aString = aString + std::to_string( aMat( k ) );
        }
    }

    // -----------------------------------------------------------------------------
    /**
     * Converts a string into a cell of vectors (1-D matrices)
     * Cells are separated by ";", vector components are separated by ","
     *
     * @tparam T Type of moris::Matrix
     * @param aString String to be converted into cell of vectors
     * @return aCellMat cell of vectors converted from string
     */

    template< typename T >
    void
    string_to_cell_mat(
            const std::string     &aString,
            Vector< Matrix< T > > &aCellMat )
    {
        if ( !aString.empty() )
        {
            uint tCellCount = std::count( aString.begin(), aString.end(), ';' ) + 1;

            aCellMat.resize( tCellCount );

            std::string tString( aString );

            uint tCount = 0;

            // reset position
            size_t tPos;

            bool tBool = true;

            while ( tBool )
            {
                // find string
                tPos = tString.find( ';' );

                if ( tPos == std::string::npos )
                {
                    tBool = false;
                }

                // reset position
                size_t tPosSubString = 0;

                uint tCount1 = 0;

                if ( tBool )
                {
                    std::string tStringMat = tString.substr( 0, tPos );

                    uint tCountNum = std::count( tStringMat.begin(), tStringMat.end(), ',' ) + 1;

                    aCellMat( tCount ).set_size( tCountNum, 1 );

                    while ( tPosSubString < tStringMat.size() )
                    {
                        // find string
                        tPosSubString = tStringMat.find( ',' );

                        // copy value into output matrix
                        if ( tPosSubString < tStringMat.size() )
                        {
                            aCellMat( tCount )( tCount1++ ) = stod( tStringMat.substr( 0, tPosSubString ) );
                            tStringMat                      = tStringMat.substr( tPosSubString + 1, tStringMat.size() );
                        }
                    }

                    // copy value into output matrix
                    aCellMat( tCount )( tCount1++ ) = stod( tStringMat );

                    tString = tString.substr( tPos + 1, tString.size() );
                }
                else
                {
                    uint tCountNum = std::count( tString.begin(), tString.end(), ',' ) + 1;

                    aCellMat( tCount ).set_size( tCountNum, 1 );

                    while ( tPosSubString < tString.size() )
                    {
                        // find string
                        tPosSubString = tString.find( ',' );

                        // copy value into output matrix
                        if ( tPosSubString < tString.size() )
                        {
                            aCellMat( tCount )( tCount1++ ) = stod( tString.substr( 0, tPosSubString ) );
                            tString                         = tString.substr( tPosSubString + 1, tString.size() );
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

    template< typename T >
    inline void
    string_to_cell_mat_2(
            const std::string     &aString,
            Vector< Matrix< T > > &aCellMat )
    {
        // if non-empty string
        if ( !aString.empty() )
        {
            // get the string to stream over
            std::istringstream tStringStream( aString );

            // get number of cells
            uint tCellCount = std::count( aString.begin(), aString.end(), '/' ) + 1;

            // set size for the list and initialize a counter
            aCellMat.resize( tCellCount );
            uint iCellCounter = 0;

            // initialize a string for individual matrices
            std::string tMatrixString;

            // initialize the value that will be inserted as the entries of the matrix
            real tValue;

            // loop ove the individual cells
            while ( std::getline( tStringStream, tMatrixString, '/' ) )
            {
                // get the matrix string to stream over
                std::istringstream tMatrixStringStream( tMatrixString );

                // set the size of the matrix
                uint tNumRows = std::count( tMatrixString.begin(), tMatrixString.end(), ';' ) + 1;
                uint tNumCols = std::count( tMatrixString.begin(), tMatrixString.end(), ',' ) + 1;
                aCellMat( iCellCounter ).set_size( tNumRows, tNumCols );

                // row  and column counter
                uint iRow = 0;
                uint iCol = 0;

                // string to traverse through each row
                std::string tRowString;

                // loop over the rows of the matrix
                while ( std::getline( tMatrixStringStream, tRowString, ';' ) )
                {
                    // convert row string to streamable string
                    std::istringstream tRowStringStream( tRowString );

                    // initialize a rwo string to store the values
                    std::string tColString;

                    // reset the column counter
                    iCol = 0;

                    // loop over the columns of each row
                    while ( std::getline( tRowStringStream, tColString, ',' ) )
                    {
                        // stream the value into a real value
                        std::istringstream tColStringStream( tColString );
                        if ( tColStringStream >> tValue )
                        {
                            // store the value
                            aCellMat( iCellCounter )( iRow, iCol++ ) = tValue;
                        }
                    }

                    // increment the row
                    iRow++;
                }

                // resize the matrix to correct size
                aCellMat( iCellCounter ).resize( iRow, iCol );

                // increment the cell counter
                iCellCounter++;
            }
        }
    }

    template< typename T >
    [[nodiscard]] inline Vector< Matrix< T > >
    string_to_cell_mat_2( const std::string &aString )
    {
        Vector< Matrix< T > > tCellMat;
        string_to_cell_mat_2( aString, tCellMat );
        return tCellMat;
    }

    // -----------------------------------------------------------------------------

    template< typename T >
    void
    string_to_cell_of_cell(
            const std::string                  &aString,
            Vector< Vector< T > >              &aCellCell,
            moris::map< std::string, T > const &aMap )
    {
        if ( !aString.empty() )
        {
            uint tCellCount = std::count( aString.begin(), aString.end(), ';' ) + 1;

            aCellCell.resize( tCellCount );

            std::string tString( aString );

            uint tCount = 0;

            // reset position
            size_t tPos;

            bool tBool = true;

            while ( tBool )
            {
                // find string
                tPos = tString.find( ';' );

                if ( tPos == std::string::npos )
                {
                    tBool = false;
                }

                // reset position
                size_t tPosSubString = 0;

                uint tCount1 = 0;

                if ( tBool )
                {
                    std::string tStringMat = tString.substr( 0, tPos );

                    uint tCountNum = std::count( tStringMat.begin(), tStringMat.end(), ',' ) + 1;

                    aCellCell( tCount ).resize( tCountNum );

                    while ( tPosSubString < tStringMat.size() )
                    {
                        // find string
                        tPosSubString = tStringMat.find( ',' );

                        // copy value into output matrix
                        if ( tPosSubString < tStringMat.size() )
                        {
                            T tComponent                     = aMap.find( tStringMat.substr( 0, tPosSubString ) );
                            aCellCell( tCount )( tCount1++ ) = tComponent;
                            tStringMat                       = tStringMat.substr( tPosSubString + 1, tStringMat.size() );
                            tPosSubString                    = tStringMat.find( ',' );
                        }
                    }

                    // copy value into output matrix
                    T tComponent                     = aMap.find( tStringMat );
                    aCellCell( tCount )( tCount1++ ) = tComponent;
                    tString                          = tString.substr( tPos + 1, tString.size() );
                }
                else
                {
                    uint tCountNum = std::count( tString.begin(), tString.end(), ',' ) + 1;

                    aCellCell( tCount ).resize( tCountNum );

                    while ( tPosSubString < tString.size() )
                    {
                        // find string
                        tPosSubString = tString.find( ',' );

                        // copy value into output matrix
                        if ( tPosSubString < tString.size() )
                        {
                            T tComponent                     = aMap.find( tString.substr( 0, tPosSubString ) );
                            aCellCell( tCount )( tCount1++ ) = tComponent;
                            tString                          = tString.substr( tPosSubString + 1, tString.size() );
                            tPosSubString                    = tString.find( ',' );
                        }
                    }

                    // copy value into output matrix
                    T tComponent                     = aMap.find( tString );
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

    template< typename T >
    [[nodiscard]] Vector< Vector< T > >
    string_to_cell_of_cell(
            const std::string                  &aString,
            moris::map< std::string, T > const &aMap )
    {
        Vector< Vector< T > > tCellCell;
        string_to_cell_of_cell( aString, tCellCell, aMap );
        return tCellCell;
    }

    // -----------------------------------------------------------------------------

    template< typename T >
    void
    string_to_cell_of_cell(
            const std::string     &aString,
            Vector< Vector< T > > &aCellCell )
    {
        if ( !aString.empty() )
        {
            uint tCellCount = std::count( aString.begin(), aString.end(), ';' ) + 1;

            aCellCell.resize( tCellCount );

            std::string tString( aString );

            uint tCount = 0;

            // reset position
            size_t tPos;

            bool tBool = true;

            while ( tBool )
            {
                // find string
                tPos = tString.find( ';' );

                if ( tPos == std::string::npos )
                {
                    tBool = false;
                }

                // reset position
                size_t tPosSubString = 0;

                uint tCount1 = 0;

                if ( tBool )
                {
                    std::string tStringMat = tString.substr( 0, tPos );

                    uint tCountNum = std::count( tStringMat.begin(), tStringMat.end(), ',' ) + 1;

                    aCellCell( tCount ).resize( tCountNum );

                    while ( tPosSubString < tStringMat.size() )
                    {
                        // find string
                        tPosSubString = tStringMat.find( ',' );

                        // copy value into output matrix
                        if ( tPosSubString < tStringMat.size() )
                        {
                            aCellCell( tCount )( tCount1++ ) = tStringMat.substr( 0, tPosSubString );
                            tStringMat                       = tStringMat.substr( tPosSubString + 1, tStringMat.size() );
                            tPosSubString                    = tStringMat.find( ',' );
                        }
                    }

                    // copy value into output matrix
                    aCellCell( tCount )( tCount1++ ) = tStringMat;
                    tString                          = tString.substr( tPos + 1, tString.size() );
                }
                else
                {
                    uint tCountNum = std::count( tString.begin(), tString.end(), ',' ) + 1;

                    aCellCell( tCount ).resize( tCountNum );

                    while ( tPosSubString < tString.size() )
                    {
                        // find string
                        tPosSubString = tString.find( ',' );

                        // copy value into output matrix
                        if ( tPosSubString < tString.size() )
                        {
                            aCellCell( tCount )( tCount1++ ) = tString.substr( 0, tPosSubString );
                            tString                          = tString.substr( tPosSubString + 1, tString.size() );
                            tPosSubString                    = tString.find( ',' );
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

    template< typename T >
    [[nodiscard]] Vector< Vector< T > >
    string_to_cell_of_cell( const std::string &aString )
    {
        Vector< Vector< T > > tCellCell;
        string_to_cell_of_cell( aString, tCellCell );
        return tCellCell;
    }

    // -----------------------------------------------------------------------------

    template< typename T >
    void
    string_to_cell(
            const std::string                  &aString,
            Vector< T >                        &aCell,
            moris::map< std::string, T > const &aMap )
    {
        if ( !aString.empty() )
        {
            uint tCellCount = std::count( aString.begin(), aString.end(), ',' ) + 1;

            aCell.resize( tCellCount );

            std::string tString( aString );

            // reset position
            size_t tPos;

            uint tCount = 0;

            bool tBool = true;

            while ( tBool )
            {
                // find string
                tPos = tString.find( ',' );

                if ( tPos == std::string::npos )
                {
                    tBool = false;
                }

                if ( tBool )
                {
                    std::string tStringMat = tString.substr( 0, tPos );

                    // check that output type is member of map
                    MORIS_ERROR( aMap.key_exists( tStringMat ),
                            "fn_Parsing_Tools::string_to_cell - key does not exist: %s",
                            tString.c_str() );

                    // copy value into output matrix
                    T tComponent = aMap.find( tStringMat );

                    aCell( tCount++ ) = tComponent;

                    tString = tString.substr( tPos + 1, tString.size() );
                }
                else
                {
                    // check that output type is member of map
                    MORIS_ERROR( aMap.key_exists( tString ),
                            "fn_Parsing_Tools::string_to_cell - key does not exist: %s",
                            tString.c_str() );

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

    template< typename T >
    [[nodiscard]] Vector< T >
    string_to_cell(
            const std::string                  &aString,
            moris::map< std::string, T > const &aMap )
    {
        Vector< T > tCell;
        string_to_cell( aString, tCell, aMap );
        return tCell;
    }

    /**
     * Converts an input string into values to be pushed back into a cell with "," delimiter.
     *
     * @tparam T Cell data type
     * @param aString Input string
     * @param aCell Cell of converted data
     */
    template< typename T >
    void string_to_cell(
            const std::string &aString,
            Vector< T >       &aCell )
    {
        // convert string to string stream and a sub string
        std::stringstream tStringStream( aString );

        // seperate the string by the delimiter
        std::string tSubString;
        while ( std::getline( tStringStream, tSubString, ',' ) )
        {
            // convert the sub string to a stream
            std::stringstream tSubStringStream( tSubString );

            // Extract the substring into the value
            T tValue;
            tSubStringStream >> tValue;

            // Add to the cell
            aCell.push_back( tValue );
        }
    }

    /**
     * Converts an input string into values to be pushed back into a cell of strings with "," delimiter.
     *
     * @param aString Input string
     * @param aCell Cell of converted data
     */
    template<>
    void string_to_cell< std::string >(
            const std::string     &aString,
            Vector< std::string > &aCell );

    /**
     * Converts an input string into a new cell to be returned.
     *
     * @tparam T Cell data type
     * @param aString Input string
     * @return Cell of converted data
     */
    template< typename T >
    [[nodiscard]] Vector< T > string_to_cell( const std::string &aString )
    {
        Vector< T > tCell;
        string_to_cell( aString, tCell );
        return tCell;
    }

    /**
     * This function inverts little endian to big endian and vice versa.
     * Needed for VTK debug files.
     */
    template< typename T >
    inline T
    swap_byte_endian( T aValue )
    {
        T     aOutValue;
        auto *tPointer    = (char *)&aValue;
        auto *tOutPointer = (char *)&aOutValue;
        int   size        = sizeof( T );
        for ( int i = 0; i < size; i++ )
        {
            tOutPointer[ size - 1 - i ] = tPointer[ i ];
        }
        return aOutValue;
    }

}    // namespace moris
