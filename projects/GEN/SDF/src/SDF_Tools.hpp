/*
 * SDF_Tools.hpp
 *
 *  Created on: Sep 30, 2018
 *      Author: messe
 */

#ifndef PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_
#define PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_

#include "typedefs.hpp"
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
            Matrix< F31RMat > aOut;

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
    }
}



#endif /* PROJECTS_GEN_SDF_SRC_SDF_TOOLS_HPP_ */
