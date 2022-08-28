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

#include "typedefs.hpp"

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_norm.hpp"

#include "fn_Parsing_Tools.hpp"

namespace moris
{
    namespace sdf
    {
        const real   gSDFepsilon = 1e-7;
//-------------------------------------------------------------------------------

        void
        TrianglePermutation(const uint aZ, uint & aX, uint & aY)
        {
            if( aZ==0 ){  // y-z
                aX = 1;
                aY = 2;
            } else if( aZ==1 ) { // x-z
                aX = 2;
                aY = 0;
            } else { // x-y
                aX = 0;
                aY = 1;
            }
        }

//-------------------------------------------------------------------------------

        Matrix< F31RMat >
        cross( const Matrix<F31RMat> & aA,  const Matrix<F31RMat> & aB )
        {
            Matrix< F31RMat > aOut( 3, 1 );

            aOut( 0 ) = aA( 1 )*aB( 2 ) - aA( 2 )*aB( 1 );
            aOut( 1 ) = aA( 2 )*aB( 0 ) - aA( 0 )*aB( 2 );
            aOut( 2 ) = aA( 0 )*aB( 1 ) - aA( 1 )*aB( 0 );

            return aOut;
        }

// =============================================================================
// LEXICAL FUNCTIONS
// =============================================================================

        template <typename T>
        T
        min(const T& aA, const T& aB)
        {
            return (aA < aB) ? (aA) : (aB);
        }

// -----------------------------------------------------------------------------

        template <typename T>
        T
        max(const T& aA, const T& aB)
        {
            return (aA > aB) ? (aA) : (aB);
        }

// -----------------------------------------------------------------------------

        template <typename T>
        T
        max(const T& aA, const T& aB, const T& aC)
        {
            return max( max(aA, aB), aC );
        }

// -----------------------------------------------------------------------------

        template <typename T>
        T
        min(const T& aA, const T& aB, const T& aC)
        {
            return min( min(aA, aB), aC );
        }

// -----------------------------------------------------------------------------

        bool
        string_to_bool( const std::string & aString )
        {
            // locale
            std::locale loc;

            // lower string of aString
            std::string tLowerString( aString );
            for( uint i=0; i < aString.length(); ++i)
            {
                tLowerString[ i ] = std::tolower( aString[i] );
            }

            return ( tLowerString == "true"
                    || tLowerString == "on"
                            || tLowerString == "yes"
                                    || tLowerString == "1" ) ;
        }

// =============================================================================
// Linear Algebra
// =============================================================================

        Matrix< F33RMat >
        rotation_matrix( const Matrix< F31RMat > & aAxis, const real & aAngle )
        {
            Matrix< F33RMat > aT( 3, 3 );

            real tCos = std::cos( aAngle );
            real tSin = std::sin( aAngle );
            real tCos_minusOne = tCos - 1.0;

            aT( 0, 0 ) = tCos-aAxis( 0 )*aAxis( 0 )*tCos_minusOne;
            aT( 1, 0 ) = aAxis( 2 )*tSin-aAxis( 0 )*aAxis( 1 )*tCos_minusOne;
            aT( 2, 0 ) = -aAxis( 1 )*tSin-aAxis( 0 )*aAxis( 2 )*tCos_minusOne;
            aT( 0, 1 ) = -aAxis( 2 )*tSin-aAxis( 0 )*aAxis( 1 )*tCos_minusOne;
            aT( 1, 1 ) = tCos-aAxis( 1 )*aAxis( 1 )*tCos_minusOne;
            aT( 2, 1 ) = aAxis( 0 )*tSin-aAxis( 1 )*aAxis( 2 )*tCos_minusOne;
            aT( 0, 2 ) = aAxis( 1 )*tSin-aAxis( 0 )*aAxis( 2 )*tCos_minusOne;
            aT( 1, 2 ) = -aAxis( 0 )*tSin-aAxis( 1 )*aAxis( 2 )*tCos_minusOne;
            aT( 2, 2 ) = tCos-aAxis( 2 )*aAxis( 2 )*tCos_minusOne;
            return aT;
        }

// =============================================================================
// Random Stuff
// =============================================================================
        uint
        random_seed()
        {
            std::ifstream file ( "/dev/urandom", std::ios::binary );
            uint tSeed;
            if ( file.is_open() )
            {
                char * memblock;
                int size = sizeof(moris::uint);
                memblock = new char [size];
                file.read (memblock, size);
                file.close();
                tSeed = *reinterpret_cast<int*>(memblock);
                delete[] memblock;
            }
            else
            {
                tSeed = time( NULL );
            }
            return tSeed;

        }

// -----------------------------------------------------------------------------

        /**
         * @brief returns a normalized pseudorandom vector
         *
         * @return    The random vector. Must be initialized already.
         */

        moris::Matrix< F31RMat >
        random_axis()
        {
            moris::Matrix< F31RMat > aVector( 3, 1 );

            std::srand( random_seed() );

            for( uint i=0; i<3; ++i )
            {
                aVector( i ) = std::rand();
            }
            real tNorm = norm( aVector );

            for( uint i=0; i<3; ++i )
            {
                aVector( i ) /= tNorm;
            }

            return aVector;
        }

// -----------------------------------------------------------------------------

        /**
         * @brief returns a pseudorandom angle between -Pi and Pi
         *
         * @return Angle: The returned angle in rad.
         *
         */
        real
        random_angle()
        {
            std::srand( random_seed() );

            return ( ( ( real )  std::rand() )/RAND_MAX - 0.5)
                    *4.0*std::acos( 0.0 );
        }

// =============================================================================
// String tools
// =============================================================================

        /**
         * this funcitons removes leading and double spaces and tabs from a string
         */
        std::string
        clean( const std::string & aString )
        {
            // length of string
            uint tLength = aString.size();

            bool tFlag = true; // flag telling if last char is space or tab

            std::string aResult = "";

            for( uint k=0; k<tLength; ++k )
            {
                if( ( int ) aString[ k ] == 9 || ( int ) aString[ k ] == 32 )
                {
                    if( ! tFlag )
                    {
                        aResult = aResult + " ";
                    }
                    tFlag = true;
                }
                else
                {
                    tFlag = false;
                    aResult = aResult + aString[ k ];
                }
            }

            // trim last space
            if( tFlag )
            {
                aResult = aResult.substr( 0, aResult.find_last_of( " " ) );
            }
            return aResult;
        }

//-------------------------------------------------------------------------------

        moris::Cell< std::string >
        string_to_words( const std::string & aString )
        {
            // output cell
            moris::Cell< std::string > aWords;

            // cleanup string
            std::string tClean = clean( aString );

            // get length of string
            int tLength = tClean.size();

            // extract words from string
            if( tLength > 0 )
            {
                int tEnd = -1;
                int tStart = 0;
                for( int k=0; k<tLength; ++k )
                {
                    if( tClean[ k ] == ' ' )
                    {
                        tStart = tEnd + 1;
                        tEnd = k;
                        aWords.push_back( tClean.substr( tStart, tEnd-tStart ) );
                    }
                }

                // last word
                tStart = tEnd + 1;
                tEnd = tLength;
                aWords.push_back( tClean.substr( tStart, tEnd-tStart ) );
            }

            return aWords;
        }

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_ */

