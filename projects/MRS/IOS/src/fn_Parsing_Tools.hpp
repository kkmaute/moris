/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_Parsing_Tools.hpp
 *
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
#include "cl_Cell.hpp"
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
    inline
    Cell<std::string>
    split_string(
            const std::string  & aString,
            const std::string  & aDelim)
    {
        // create empty cell of strings
        Cell<std::string> tCellOfStrings;

        size_t start;
        size_t end = 0;

        while ((start = aString.find_first_not_of(aDelim, end)) != std::string::npos)
        {
            end = aString.find(aDelim, start);
            tCellOfStrings.push_back(aString.substr(start, end - start));
        }

        return tCellOfStrings;
    }

    // -----------------------------------------------------------------------------
    /**
     * removes leading whitespace of string
     *
     * @param s std::string
     */
    inline
    void ltrim_string(std::string& aString )
    {
        auto it = std::find_if(aString.begin(), aString.end(),
                [](char c) {
            return !std::isspace<char>(c, std::locale::classic());
        });
        aString.erase(aString.begin(), it);
    }

    // -----------------------------------------------------------------------------
    /**
     * removes trailing whitespace of string
     *
     * @param aString std::string
     */
    inline
    void rtrim_string(std::string& aString)
    {
        auto it = std::find_if(aString.rbegin(), aString.rend(),
                [](char c) {
            return !std::isspace<char>(c, std::locale::classic());
        });
        aString.erase(it.base(), aString.end());
    }

    // -----------------------------------------------------------------------------
    /**
     * removes leading and trailing whitespace of string
     *
     * @param aString std::string
     */
    inline
    void trim_string(std::string& aString)
    {
        // if aString is empty return
        if  ( aString.length() == 0 )
        {
            return;
        }

        // trim off trailing whitespaces
        rtrim_string(aString);

        // trim off leading whitespaces
        ltrim_string(aString);
    }

    // -----------------------------------------------------------------------------
    /**
     * splits a string in sub-strings defined by a set of delimiters, removes leading and trailing whitespace
     *  of sub-string, and reassembles sub-strings into single string; if not delimiter is provides input string
     *  will be trimmed.
     *
     * @param aString std::string
     * @param aDelimiter const std::string
     */
    inline
    void split_trim_string(
            std::string       & aString,
            const std::string & aDelimiter)
    {
        char tDelim[] = " ";

        // if aString is empty return
        if  ( aString.empty() )
        {
            return;
        }
        // if no delimiter is provided just trim the given string
        if ( aDelimiter.empty() )
        {
            // trim substring
            trim_string(aString);
            return;
        }

        // loop over all delimiter characters
        for (uint id=0;id<aDelimiter.length();id++)
        {
            // check if delimiter is first element of string
            bool isFirst = aString.front() ==aDelimiter[id] ? true : false;

            // check if delimiter is last element of string
            bool isLast = aString.back() == aDelimiter[id] ? true : false;

            // get char* of current delimiter
            strncpy(tDelim,&aDelimiter[id],1);

            // allocate vector for storing substrings
            std::vector<std::string> tSubStringVec;

            // get first substring
            char *token = strtok(const_cast<char*>(aString.c_str()), tDelim);

            // loop over all substrings
            while (token != nullptr)
            {
                // create substring
                std::string tSubstring(token);

                // trim substring
                trim_string(tSubstring);

                // check if substring is empty or has only white spaces
                if ( ! tSubstring.empty())
                {
                    if ( ! std::all_of(tSubstring.begin(), tSubstring.end(), isspace) )
                    {
                        tSubStringVec.push_back(tSubstring);
                    }
                }

                // get next substring substring
                token = strtok(nullptr, tDelim);
            }

            // reassemble substrings
            aString.clear();

            if (isFirst) aString = tDelim;

            if ( tSubStringVec.size() > 0 )
            {
                aString+=tSubStringVec[0];

                for (uint is=1;is<tSubStringVec.size();++is)
                {
                    aString+=tDelim;
                    aString+=tSubStringVec[is];
                }
            }
            if (isLast) aString += tDelim;
        }
    }

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
    template < typename T >
    void delimited_string_to_mat(
            const std::string & aString,
            const std::string & aDelim,
            Matrix< T >       & aMat )
    {
        // check for empty string
        if (aString.size() == 0 )
        {
            aMat.set_size(0,0);
            return;
        }

        // split string into substrings
        Cell<std::string> tSubStringVec = split_string(aString,aDelim);

        // get number of substrings
        uint tNumberOfSubStrings = tSubStringVec.size();

        if ( tNumberOfSubStrings == 0)
        {
            aMat.set_size(0,0);
            return;
        }

        // set size of matrix
        aMat.set_size(tNumberOfSubStrings,1);

        // convert strings to numerical value
        for (uint tStrgIndex=0; tStrgIndex<tNumberOfSubStrings; ++tStrgIndex)
        {
            aMat(tStrgIndex) =  stod( tSubStringVec(tStrgIndex) );
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

    template < typename T >
    void string_to_mat(
            const std::string & aString,
            Matrix< T >       & aMat )
    {
        if( aString.size() > 0 )
        {
            uint tCountRow = std::count( aString.begin(), aString.end(), ';') + 1;
            uint tCountCol = (std::count( aString.begin(), aString.end(), ',') / tCountRow) + 1;

            std::string tString( aString );

            // allocate memory
            aMat.set_size( tCountRow, tCountCol );

            // reset position
            size_t tPos = 0;

            // reset string
            tString = aString;

            // Create output matrix
            uint tRowIndex;
            uint tColIndex;
            for (tRowIndex = 0; tRowIndex < tCountRow; tRowIndex++)
            {
                for (tColIndex = 0; tColIndex < tCountCol - 1; tColIndex++)
                {
                    // find string
                    tPos = tString.find( "," );

                    // copy value into output matrix
                    if( tPos <  tString.size() )
                    {
                        aMat(tRowIndex, tColIndex) = stod(  tString.substr( 0, tPos ) );
                        tString =  tString.substr( tPos+1, tString.size() );
                    }
                }

                // Last value before next row
                tPos = tString.find( ";" );

                // copy value into output matrix
                if( tPos <  tString.size() )
                {
                    aMat(tRowIndex, tColIndex) = stod(  tString.substr( 0, tPos ) );
                    tString =  tString.substr( tPos+1, tString.size() );
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

    template < typename T >
    Matrix<T> string_to_mat( const std::string & aString)
    {
        Matrix<T> tMat;
        string_to_mat(aString, tMat);

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

    template < typename T >
    void mat_to_string(
            const Matrix< T > & aMat,
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
    /**
     * Converts a string into a cell of vectors (1-D matrices)
     * Cells are separated by ";", vector components are separated by ","
     *
     * @tparam T Type of moris::Matrix
     * @param aString String to be converted into cell of vectors
     * @return aCellMat cell of vectors converted from string
     */

    template < typename T >
    void string_to_cell_mat(
            const std::string   & aString,
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
    inline
    void string_to_cell_mat_2(
            const std::string   & aString,
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
    void string_to_cell_of_cell(
            const std::string               & aString,
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
                            tPosSubString = tStringMat.find( "," );
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
                            tPosSubString = tString.find( "," );
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

    // -----------------------------------------------------------------------------

    template < typename T >
    void string_to_cell_of_cell(
            const std::string               & aString,
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
                            tPosSubString = tStringMat.find( "," );
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
                            tString = tString.substr( tPosSubString+1, tString.size() );
                            tPosSubString = tString.find( "," );
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

    // -----------------------------------------------------------------------------

    template < typename T >
    void string_to_cell(
            const std::string             & aString,
            moris::Cell< T >              & aCell,
            moris::map< std::string, T >  & aMap )
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

                    // check that output type is member of map
                    MORIS_ERROR( aMap.key_exists(tStringMat),
                            "fn_Parsing_Tools::string_to_cell - key does not exist: %s",tString.c_str() );

                    // copy value into output matrix
                    T tComponent = aMap.find( tStringMat );

                    aCell( tCount++ ) = tComponent;

                    tString =  tString.substr( tPos+1, tString.size() );
                }
                else
                {
                    // check that output type is member of map
                    MORIS_ERROR( aMap.key_exists(tString),
                            "fn_Parsing_Tools::string_to_cell - key does not exist: %s",tString.c_str() );

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

    // -----------------------------------------------------------------------------

    template < typename T >
    void string_to_cell(
            const std::string & aString,
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
    }

    //------------------------------------------------------------------------------------------------------------------

    template < typename T >
    inline
    moris::Cell<T> string_to_cell( const std::string & aString)
    {
        moris::Cell<T> tCell;
        string_to_cell(aString, tCell);
        return tCell;
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

    /**
     * This function inverts little endian to big endian and vice versa.
     * Needed for VTK debug files.
     */
    template <typename T>  inline T swap_byte_endian(T aValue)
    {
        T aOutValue;
        auto *tPointer = (char*) &aValue;
        auto *tOutPointer = (char*)&aOutValue;
        int size = sizeof(T);
        for(int i=0; i<size; i++)
        {
            tOutPointer[size - 1 - i] = tPointer[i];
        }
        return aOutValue;
    }

}

#endif	/* MORIS_IOS_FN_PARSING_TOOLS_HPP_ */

