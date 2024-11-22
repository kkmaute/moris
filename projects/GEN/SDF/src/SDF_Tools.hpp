/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * SDF_Tools.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_
#define PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_

#include <fstream>
#include <cmath>
#include <limits>

#include "moris_typedefs.hpp"

#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Communication_Tools.hpp"      // COM/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_norm.hpp"

#include "fn_Parsing_Tools.hpp"

namespace moris::sdf
{      
    //-------------------------------------------------------------------------------

    inline void
    triangle_permutation( const uint aZ, uint& aX, uint& aY )
    {
        switch ( aZ )
        {
            case 0:
                // y-z
                aX = 1;
                aY = 2;
                break;
            case 1:
                // z-x
                aX = 2;
                aY = 0;
                break;
            case 2:
                // x-y
                aX = 0;
                aY = 1;
                break;
            default:
                MORIS_ERROR( false, "triangle_permuation() expects an input less than 3" );
                break;
        }
    }

    //-------------------------------------------------------------------------------

    inline Matrix< F31RMat >
    cross( const Matrix< F31RMat >& aA, const Matrix< F31RMat >& aB )
    {
        Matrix< F31RMat > aOut( 3, 1 );

        aOut( 0 ) = aA( 1 ) * aB( 2 ) - aA( 2 ) * aB( 1 );
        aOut( 1 ) = aA( 2 ) * aB( 0 ) - aA( 0 ) * aB( 2 );
        aOut( 2 ) = aA( 0 ) * aB( 1 ) - aA( 1 ) * aB( 0 );

        return aOut;
    }

    // =============================================================================
    // LEXICAL FUNCTIONS
    // =============================================================================

    template< typename T >
    inline T
    min( const T& aA, const T& aB )
    {
        return ( aA < aB ) ? ( aA ) : ( aB );
    }

    // -----------------------------------------------------------------------------

    template< typename T >
    inline T
    max( const T& aA, const T& aB )
    {
        return ( aA > aB ) ? ( aA ) : ( aB );
    }

    // -----------------------------------------------------------------------------

    template< typename T >
    inline T
    max( const T& aA, const T& aB, const T& aC )
    {
        return max( max( aA, aB ), aC );
    }

    // -----------------------------------------------------------------------------

    template< typename T >
    inline T
    min( const T& aA, const T& aB, const T& aC )
    {
        return min( min( aA, aB ), aC );
    }

    // -----------------------------------------------------------------------------

    bool inline string_to_bool( const std::string& aString )
    {
        // locale
        std::locale loc;

        // lower string of aString
        std::string tLowerString( aString );
        for ( uint i = 0; i < aString.length(); ++i )
        {
            tLowerString[ i ] = std::tolower( aString[ i ] );
        }

        return ( tLowerString == "true"
                 || tLowerString == "on"
                 || tLowerString == "yes"
                 || tLowerString == "1" );
    }

    // =============================================================================
    // Linear Algebra
    // =============================================================================

    inline Matrix< DDRMat >
    rotation_matrix( const Matrix< DDRMat >& aAxis, const real& aAngle )
    {
        Matrix< DDRMat > aT( 3, 3 );
        real             tCos          = std::cos( aAngle );
        real             tSin          = std::sin( aAngle );
        real             tCos_minusOne = tCos - 1.0;

        aT( 0, 0 ) = tCos - aAxis( 0 ) * aAxis( 0 ) * tCos_minusOne;
        aT( 1, 0 ) = aAxis( 2 ) * tSin - aAxis( 0 ) * aAxis( 1 ) * tCos_minusOne;
        aT( 2, 0 ) = -aAxis( 1 ) * tSin - aAxis( 0 ) * aAxis( 2 ) * tCos_minusOne;
        aT( 0, 1 ) = -aAxis( 2 ) * tSin - aAxis( 0 ) * aAxis( 1 ) * tCos_minusOne;
        aT( 1, 1 ) = tCos - aAxis( 1 ) * aAxis( 1 ) * tCos_minusOne;
        aT( 2, 1 ) = aAxis( 0 ) * tSin - aAxis( 1 ) * aAxis( 2 ) * tCos_minusOne;
        aT( 0, 2 ) = aAxis( 1 ) * tSin - aAxis( 0 ) * aAxis( 2 ) * tCos_minusOne;
        aT( 1, 2 ) = -aAxis( 0 ) * tSin - aAxis( 1 ) * aAxis( 2 ) * tCos_minusOne;
        aT( 2, 2 ) = tCos - aAxis( 2 ) * aAxis( 2 ) * tCos_minusOne;
        return aT;
    }

    // Creates a 2D rotation matrix of angle aAngle
    inline Matrix< DDRMat >
    rotation_matrix( const real& aAngle )
    {
        Matrix< DDRMat > aT( 2, 2 );
        real             tCos = std::cos( aAngle );
        real             tSin = std::sin( aAngle );
        aT( 0, 0 )            = tCos;
        aT( 0, 1 )            = -tSin;
        aT( 1, 0 )            = tSin;
        aT( 1, 1 )            = tCos;
        return aT;
    }

    // =============================================================================
    // Random Stuff
    // =============================================================================
    inline uint
    random_seed()
    {
        std::ifstream file( "/dev/urandom", std::ios::binary );
        uint          tSeed;
        if ( file.is_open() )
        {
            char* tMemBlock;
            int   size = sizeof( moris::uint );
            tMemBlock  = new char[ size ];
            file.read( tMemBlock, size );
            file.close();
            tSeed = *reinterpret_cast< int* >( tMemBlock );
            delete[] tMemBlock;
        }
        else
        {
            tSeed = time( nullptr );
        }
        return tSeed;
    }

    // -----------------------------------------------------------------------------

    /**
     * @brief returns a normalized pseudorandom vector
     *
     * @return    The random vector. Must be initialized already.
     */

    inline moris::Matrix< DDRMat >
    random_axis( uint tNumDim )
    {
        moris::Matrix< DDRMat > tVector( tNumDim, 1 );

        std::srand( random_seed() );

        // generate random number for each entry of the return vector
        for ( uint i = 0; i < tNumDim; ++i )
        {
            tVector( i ) = std::rand();
        }

        // compute the norm of the vector
        real tNorm = norm( tVector );

        // Create a unit vector by dividing each entry by the norm
        for ( uint i = 0; i < tNumDim; ++i )
        {
            tVector( i ) /= tNorm;
        }

        return tVector;
    }

    // -----------------------------------------------------------------------------

    /**
     * @brief returns a pseudorandom angle between -Pi and Pi
     *
     * @return Angle: The returned angle in rad.
     *
     */
    inline real
    random_angle()
    {
        std::srand( random_seed() );

        return ( ( (real)std::rand() ) / RAND_MAX - 0.5 )
             * 4.0 * std::acos( 0.0 );
    }

    // =============================================================================
    // String tools
    // =============================================================================

    /**
     * this function removes leading and double spaces and tabs from a string
     */
    inline std::string
    clean( const std::string& aString )
    {
        // length of string
        uint tLength = aString.size();

        bool tFlag = true;    // flag telling if last char is space or tab

        std::string aResult = "";

        for ( uint k = 0; k < tLength; ++k )
        {
            if ( (int)aString[ k ] == 9 || (int)aString[ k ] == 32 )
            {
                if ( !tFlag )
                {
                    aResult = aResult + " ";
                }
                tFlag = true;
            }
            else
            {
                tFlag   = false;
                aResult = aResult + aString[ k ];
            }
        }

        // trim last space
        if ( tFlag )
        {
            aResult = aResult.substr( 0, aResult.find_last_of( ' ' ) );
        }
        return aResult;
    }

    //-------------------------------------------------------------------------------

    inline Vector< std::string >
    string_to_words( const std::string& aString )
    {
        // output cell
        Vector< std::string > aWords;

        // cleanup string
        std::string tClean = clean( aString );

        // get length of string
        int tLength = tClean.size();

        // extract words from string
        if ( tLength > 0 )
        {
            int tEnd   = -1;
            int tStart = 0;
            for ( int k = 0; k < tLength; ++k )
            {
                if ( tClean[ k ] == ' ' )
                {
                    tStart = tEnd + 1;
                    tEnd   = k;
                    aWords.push_back( tClean.substr( tStart, tEnd - tStart ) );
                }
            }

            // last word
            tStart = tEnd + 1;
            tEnd   = tLength;
            aWords.push_back( tClean.substr( tStart, tEnd - tStart ) );
        }

        return aWords;
    }

    //-------------------------------------------------------------------------------
}    // namespace moris::sdf

#endif /* PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_ */
