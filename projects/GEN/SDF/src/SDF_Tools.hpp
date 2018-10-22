/*
 * SDF_Tools.hpp
 *
 *  Created on: Sep 30, 2018
 *      Author: messe
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

        /**
         * This function inverts little endian to big endian and vice versa.
         * Needed for VTK debug files.
         */
        template <typename T> T swap_byte_endian(T aValue)
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

// -----------------------------------------------------------------------------

        std::string
        parallelize_path( const std::string & aFilePath )
        {
            if( par_size() == 1 || aFilePath.size() == 0 )
            {
                // leave path untouched
                return aFilePath;
            }
            else
            {
                return        aFilePath.substr(0,aFilePath.find_last_of(".")) // base path
                        + "." + std::to_string( par_size() ) // rank of this processor
                + "." + std::to_string( par_rank() ) // number of procs
                +  aFilePath.substr( aFilePath.find_last_of("."), aFilePath.length() ); // file extension
            }
        }

// -----------------------------------------------------------------------------

        std::string
        proc_string()
        {
            std::string tString = "              ";

            if( par_size() > 1 )
            {
                uint tMyRank = par_rank();
                tString = "  proc " + std::to_string( tMyRank );

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

}



#endif /* PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_ */
