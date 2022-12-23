/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Fortran2.hpp
 *
 */

#pragma once

#include <string.h>

//

/*
 implemented from
  http://arnholm.org/software/cppf77/cppf77.htm#Section3

  Mixed language programming
    using C++ and FORTRAN 77 by Carsten A. Arnholm
*/

namespace moris::fortran
{
    // type alias for fortran types
    using INTEGER          = int;
    using LOGICAL          = int;
    using DOUBLE_PRECISION = double;

    // fortran character class
    class CHARACTER
    {
      public:
        CHARACTER( char* aString );
        CHARACTER( char* aString, const size_t aStrLen );
        ~CHARACTER();
        CHARACTER operator()( size_t aIndex );
        void      pad( size_t aFirst, size_t aNumber = 1 );
        void      operator=( char* aStr );
                  operator char*();

      public:
        char*  mStr;    // Actual string
        size_t mLen;    // String length
    };


    inline CHARACTER::CHARACTER( char* aString )
            : mStr( aString )
            , mLen( strlen( aString ) )
    {
    }

    inline CHARACTER::CHARACTER( char* aString, const size_t lstr )
            : mStr( aString )
            , mLen( lstr )
    {
        // find position from where to start padding
        size_t slen   = strlen( mStr );                             // upper limit
        size_t actual = ( slen < mLen ) ? slen : mLen;               // actual <= len.
        for ( size_t i = actual; i < mLen; i++ ) mStr[ i ] = ' ';    // Do the padding.
    }

    inline CHARACTER::~CHARACTER()
    {
        if ( mStr[ mLen ] == '\0' ) return;    // catches string constants

        for ( int i = mLen - 1; i >= 0; i-- )
        {
            if ( mStr[ i ] == '\0' ) break;    // already zero terminated

            if ( mStr[ i ] != ' ' )
            {                           // non-blank discovered, so
                mStr[ i + 1 ] = '\0';    // zero-terminate and jump out
                break;
            }
        }
    }

    inline CHARACTER
    CHARACTER::operator()( size_t aIndex )
    {
        // Construct a temporary CHARACTER object for the array element
        // identified by "index" in order to zero-terminate that element
        size_t    tPos = aIndex * mLen;            // start pos of array element
        CHARACTER tElement( mStr + tPos, mLen );    // construct new CHARACTER.
        return tElement;                         // destructor called here.
    }

    inline void
    CHARACTER::pad( size_t aFirst, size_t aNumber )
    {

        size_t tPos = 0; 
        size_t tStop = aFirst + aNumber - 1;

        for ( size_t index = aFirst; index <= tStop; index++ )
        {
            tPos           = index * mLen;
            size_t slen   = strlen( mStr + tPos );    // upper limit
            size_t actual = ( slen < mLen ) ? slen : mLen;
            for ( size_t i = tPos + actual; i < tPos + mLen; i++ ) mStr[ i ] = ' ';    // Do the padding.
        }
    }

    inline void
    CHARACTER::operator=( char* aStr )
    {
        strncpy( mStr, aStr, mLen );                                  // this will copy a zero if str < mStr
        mStr[ mLen - 1 ] = '\0';                                     // zero terminate in case strncpy did not
        size_t slen    = strlen( mStr );                            // upper limit
        size_t actual  = ( slen < mLen ) ? slen : mLen;              // actual <= len.
        for ( size_t i = actual; i < mLen; i++ ) mStr[ i ] = ' ';    // Do the padding.
    }

    inline CHARACTER::operator char*()
    {
        return mStr;
    }
}    // namespace moris::FORTRAN
