/*
 * SDF_Tools.hpp
 *
 *  Created on: Sep 30, 2018
 *      Author: messe
 */

#ifndef PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_
#define PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_


#include "typedefs.hpp"

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

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
// -----------------------------------------------------------------------------
    }
}



#endif /* PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_ */
