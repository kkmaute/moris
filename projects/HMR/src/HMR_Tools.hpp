/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_Tools.hpp
 *
 */

#pragma once

#include <string>         // std::string
#include <locale>         // std::locale, std::tolower

#include <iostream>

#include "cl_Communication_Tools.hpp" //COM/src
#include "moris_typedefs.hpp" //COR/src
#include "cl_Map.hpp" //CNT/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_unique.hpp"

#include "fn_Parsing_Tools.hpp"

namespace moris::hmr
{
// -----------------------------------------------------------------------------
// print dots for nice output

inline
std::string proc_string()
{
//        std::string tString = "              ";
    std::string tString = "#";

    if( par_size() > 1 )
    {
        uint tMyRank = par_rank();

        tString = "proc " + std::to_string( tMyRank );

        if ( tMyRank < 10 )
        {
            tString +=" ... :" ;
        }
        else if ( tMyRank < 100 )
        {
            tString +=" .. :" ;
        }
        else if ( tMyRank < 1000 )
        {
            tString +=" . :" ;
        }
        else if ( tMyRank < 10000 )
        {
            tString +="  :" ;
        }
        else
        {
            tString +=" :" ;
        }
    }

    return tString;
}

// -----------------------------------------------------------------------------

/*   template <typename T, typename U >
void unique_with_matrix( Mat< T > & aVector, Mat< U > & aMatrix )
{
    // backup input vector
    Mat< T > tVector = aVector;

    // backup input Matrix
    Mat< U > tMatrix = aMatrix;

    // make vector unique;
    aVector = unique( tVector );

    // get number of rows
    luint tRows = tMatrix.n_rows();

    // get number of cols
    luint tCols = tMatrix.n_cols();

    // get length of unique vector
    luint tLength = aVector.length();

    // reset size of output matrix
    aMatrix.set_size( tRows, tLength  );

    // create map
    map< T, T > tMap;
    for( uint j=0; j<tCols; ++j )
    {
        tMap[ tVector(j) ] = j;
    }

    // loop over all new entries
    for( luint k=0; k<tLength; ++k )
    {
        // find position in old vector
        luint j = tMap.find( aVector( k ) );

        // copy column
        aMatrix.cols( k, k ) = tMatrix.cols( j, j );
    }
} */

// -----------------------------------------------------------------------------

inline
bool string_to_bool( const std::string & aString )
{
    // locale
    std::locale loc;

    // lower string of aString
    std::string tLowerString( aString );
    for( uint i=0; i < aString.length(); ++i)
    {
        tLowerString[ i ] = std::tolower( aString[i] );
    }

    return ( tLowerString == "true" || tLowerString == "on"
          || tLowerString == "yes"  || tLowerString == "1" ) ;
}

// -----------------------------------------------------------------------------

    /**
     * returns the binominalcoefficient of n over k as real
     */
    inline
    real nchoosek( uint aN, const uint aK )
    {
        real aResult = 1.0;

        for ( uint i = 1; i <= aK; ++i )
        {
            aResult *= ( ( real ) aN + 1 - i ) / ( real( i ) );
        }

        return aResult;
    }

// -----------------------------------------------------------------------------
} /* namespace moris */
