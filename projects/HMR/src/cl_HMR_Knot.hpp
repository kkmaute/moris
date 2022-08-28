/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Knot.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_KNOT_HPP_
#define SRC_HMR_CL_HMR_KNOT_HPP_

#include <cmath>
#include "typedefs.hpp" //COR/src

namespace moris
{
    namespace hmr
    {

    // https://en.wikipedia.org/wiki/Figure-eight_knot_(mathematics)
    class Knot
     {
         real mParam0[ 4 ];
         real mParam1[ 4 ];

         // four knot parameters
         real mA;
         real mB;
         real mC;
         real mD;

         //! radius of tube that goes around the knot
         real mRmax = 1;
         real mRmin = 0.5;
         real mRmid = 0.75;

         real mXi;
         uint mNumberOfPoints = 200;
         Matrix< DDRMat > mPoints;
// ----------------------------------------------------------------------------
     public:
// ----------------------------------------------------------------------------

         Knot(  const real & aA0,
                const real & aB0,
                const real & aC0,
                const real & aD0,
                const real & aA1,
                const real & aB1,
                const real & aC1,
                const real & aD1,
                const real & aXi )
         {
             mParam0[ 0 ] = aA0;
             mParam0[ 1 ] = aB0;
             mParam0[ 2 ] = aC0;
             mParam0[ 3 ] = aD0;

             mParam1[ 0 ] = aA1;
             mParam1[ 1 ] = aB1;
             mParam1[ 2 ] = aC1;
             mParam1[ 3 ] = aD1;

             mXi   = aXi;

             mPoints.set_size( 4, mNumberOfPoints );
             real tDeltaT = 10.0/mNumberOfPoints;
             real tT = 0;

             for( uint k=0; k<mNumberOfPoints; ++k )
             {
                 mPoints( 0, k ) = tT;

                 this->calc_point_on_knot(
                         mPoints( 0, k ),
                         mPoints( 1, k ),
                         mPoints( 2, k ),
                         mPoints( 3, k ) );
                 tT += tDeltaT;
             }

         };

// ----------------------------------------------------------------------------
         ~Knot(){};
// ----------------------------------------------------------------------------

         void
         set_radius( const real & aRadius )
         {
             mRmax = aRadius;
             mRmin = 0.75*aRadius;
             mRmid = 0.5*(mRmax+mRmin);
         }

// ----------------------------------------------------------------------------

         real
         test_for_intersection(
                 const real & aX,
                 const real & aY,
                 const real & aZ ) const
         {
             // call distance calculation
             real tR = this->calculate_distance( aX, aY, aZ );
             //real tR = this->distance_to_sphere( aX, aY, aZ );

             if ( tR <= mRmax )
             {
                 return -1;
             }
             else
             {
                 return 1;
             }
         }

// ----------------------------------------------------------------------------

         real
         calculate_distance(
                 const real & aX,
                 const real & aY,
                 const real & aZ ) const
         {
             // make sure that newton does not hit a local minimum

             uint tK =  this->find_initial_param( aX, aY, aZ );

             real tDeltaT = mPoints( 0, tK+1 ) -mPoints( 0, tK-1 );
             real tT = mPoints( 0, tK );

             for( uint k=0; k<16; ++k )
             {
                 tT = this->find_better_t( tT, tDeltaT, aX, aY, aZ );
                 tDeltaT *= 0.25;
             }
             // calculate point on curve
             real tX;
             real tY;
             real tZ;
             this->calc_point_on_knot( tT, tX, tY, tZ );

             // return distance
             return this->distance( aX, aY, aZ, tX, tY, tZ );
         }

// ----------------------------------------------------------------------------

         real
         calculate_level_set(
                 const real & aX,
                 const real & aY,
                 const real & aZ )  const
         {
             real tR = this->calculate_distance( aX, aY, aZ );

             // calculate distance
             //real tR = this->distance_to_sphere( aXYZ[ 0 ], aXYZ[ 1 ], aXYZ[ 2 ] );

             if ( tR < mRmid )
             {
                 return mRmin - tR;
             }
             else
             {
                 return tR - mRmax;
             }

         }

// ----------------------------------------------------------------------------
     private:
// ----------------------------------------------------------------------------

         real distance_to_sphere(
                 const real & aX,
                 const real & aY,
                 const real & aZ ) const
         {
             real tMinDistance = 1E12;
             for ( uint k = 0; k<mNumberOfPoints; ++k )
             {
                 real tDistance = this->distance(
                         aX, aY, aZ,
                         mPoints( 1, k ), mPoints( 2, k ), mPoints( 3, k ) );

                 if ( tDistance < tMinDistance )
                 {
                     tMinDistance = tDistance;
                 }

             }
             return tMinDistance;
         }

// ----------------------------------------------------------------------------

         real
         find_better_t( const real& aT,
                        const real& aDeltaT,
                        const real & aX,
                        const real & aY,
                        const real & aZ ) const
         {
             real tT0 = aT - aDeltaT;
             real tT1 = aT + aDeltaT;

             real tStep = 6;
             real tDeltaT = (tT1-tT0)/tStep;

             real tRbest = 1E12;
             real tTbest = aT;
             real tT = tT0;

             real tX;
             real tY;
             real tZ;

             for( tT = tT0; tT<=tT1; tT=tT+tDeltaT )
             {
                 this->calc_point_on_knot( tT, tX, tY, tZ );
                 real tR = this->distance( tX, tY, tZ, aX, aY, aZ );
                 if ( tR < tRbest )
                 {
                     tTbest = tT;
                     tRbest = tR;
                 }
             }
             return tTbest;
         }

// ----------------------------------------------------------------------------

         uint find_initial_param(
                 const real & aX,
                 const real & aY,
                 const real & aZ ) const
         {
             real tMinDistance = 1E12;
             uint aK = mNumberOfPoints;
             for ( uint k = 1; k<mNumberOfPoints-1; ++k )
             {
                 real tDistance = this->distance(
                         aX, aY, aZ,
                         mPoints( 1, k ), mPoints( 2, k ), mPoints( 3, k ) );

                 if ( tDistance < tMinDistance  )
                 {
                     tMinDistance = tDistance;
                     aK = k;
                 }
             }

             return aK;
         }

// ----------------------------------------------------------------------------

         real
         bisection( const uint & aK,
                    const real & aX,
                    const real & aY,
                    const real & aZ ) const
         {
             // left value
             real tT0 = mPoints( 0, aK-1 ) - 0.5;

             // right value
             real tT1 = mPoints( 0, aK+1 ) + 0.5;

             real tT = tT1;

             real tF0  = calc_second_derivative( tT0, aX, aY, aZ );

             uint tCount = 0;
             // bisection
             while ( std::abs( tT0 - tT1 ) > 1e-4 )
             {
                 tT = 0.5*( tT0 + tT1 );
                 real tFt  = calc_second_derivative( tT, aX, aY, aZ ) ;

                 if ( tF0*tFt > 0 )
                 {
                     tT0 = tT;
                     tF0 = tFt;
                 }
                 else
                 {
                     tT1 = tT;
                 }
                 ++tCount;
                 //std::cout << tCount << " " << tT << " " << tFt << std::endl;
             }

             return tT;
         }

// ----------------------------------------------------------------------------
         /**
          * calculates a point from a given parameter
          */
         void
         calc_point_on_knot(
                 const real & aT,
                 real & aX,
                 real & aY,
                 real & aZ ) const
         {
             // knot 0
             real tCosB = std::cos( mParam0[ 1 ] * aT );
             real tCosC = std::cos( mParam0[ 2 ] * aT );
             real tSinC = std::sin( mParam0[ 2 ] * aT );
             real tSinD = std::sin( mParam0[ 3 ] * aT );

             aX = ( mParam0[ 0 ] + tCosB ) * tCosC * ( 1 - mXi );
             aY = ( mParam0[ 0 ] + tCosB ) * tSinC * ( 1 - mXi );
             aZ = tSinD * ( 1 - mXi );

             // knot 1
             tCosB = std::cos( mParam1[ 1 ] * aT );
             tCosC = std::cos( mParam1[ 2 ] * aT );
             tSinC = std::sin( mParam1[ 2 ] * aT );
             tSinD = std::sin( mParam1[ 3 ] * aT );

             aX += ( mParam1[ 0 ] + tCosB ) * tCosC * mXi;
             aY += ( mParam1[ 0 ] + tCosB ) * tSinC * mXi;
             aZ += tSinD * mXi;
         }

// ----------------------------------------------------------------------------

 /*        real newton(
                 const real& aOmega,
                 const real& aT,
                 const real& aX,
                 const real& aY,
                 const real& aZ ) const
         {
             real df0dt;
             real df1dt;
             real d2f0dt2;
             real d2f1dt2;

             // knot 0
             this->calc_derivatives(
                     mParam0, aT, aX, aY, aZ, df0dt, d2f0dt2 );

             // knot 1
             this->calc_derivatives(
                     mParam1, aT, aX, aY, aZ, df1dt, d2f1dt2 );

             // combine derivatives
             real dfdt   = df0dt * ( 1 - mXi ) + df1dt * mXi;
             real d2fdt2 = d2f0dt2 * ( 1 - mXi ) + d2f1dt2 * mXi;

             // perform Newton step

             if ( std::abs( d2fdt2 ) < 1e-4 )
             {
                 return aT;
             }
             else
             {
                 return aT - aOmega * dfdt/d2fdt2;
             }

         } */
// ----------------------------------------------------------------------------
         real
         calc_second_derivative(
                 const real& aT,
                 const real& aX,
                 const real& aY,
                 const real& aZ ) const
         {
            long double tDeltaT = 0.001;

            long double tT0 = aT-tDeltaT ;
            long double tT2 = aT+tDeltaT;

            real tX;
            real tY;
            real tZ;

            this->calc_point_on_knot( tT0, tX, tY, tZ );
            real tR0 = distance ( tX, tY, tZ, aX, aY, aZ );

            this->calc_point_on_knot( aT, tX, tY, tZ );
            real tR1 = distance ( tX, tY, tZ, aX, aY, aZ );

            this->calc_point_on_knot( tT2, tX, tY, tZ );
            real tR2 = distance ( tX, tY, tZ, aX, aY, aZ );

            return ( tR0 - 2*tR1 + tR2 )/std::pow( tDeltaT, 2 );
         }

// ----------------------------------------------------------------------------

         /**
          * the distance between two points
          */
         real distance (
                 const real & aX0,
                 const real & aY0,
                 const real & aZ0,
                 const real & aX1,
                 const real & aY1,
                 const real & aZ1 ) const
         {
             return std::sqrt(
                      std::pow( ( aX1 - aX0 ), 2 )
                    + std::pow( ( aY1 - aY0 ), 2 )
                    + std::pow( ( aZ1 - aZ0 ), 2 ) );
         }

// ----------------------------------------------------------------------------
     };

// ----------------------------------------------------------------------------
    }
}
#endif /* SRC_HMR_CL_HMR_KNOT_HPP_ */

