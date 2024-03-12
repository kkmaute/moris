/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Interpolaton.hpp
 *
 */

#ifndef SRC_TOOLS_CL_INTERPOLATON_HPP_
#define SRC_TOOLS_CL_INTERPOLATON_HPP_

// Matrix Include
#include "cl_Matrix.hpp"

#include "fn_print.hpp"
namespace moris::xtk
{
    class Interpolation
    {
      public:
        // http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch16.d/IFEM.Ch16.pdf
        /**
         * Linear interpolation of a value given a location
         *
         * @param[in] aInterpVars - Interpolation Vars (x,y,z are treated as independent interpolation variables
         * @param[in] aLclCoords  - Local coordinates to interpolate to (a point at the center of edge has {{0}}
         */

        template< typename Real_Matrix >
        static void linear_interpolation_location( const Matrix< Real_Matrix > &aInterpVars,
                const Matrix< Real_Matrix >                                    &aLocation,
                Matrix< Real_Matrix >                                          &aInterpolationResult )
        {
            typename Matrix< Real_Matrix >::Data_Type xi                    = aLocation( 0, 0 );
            size_t                                    tNumInterpolationVars = aInterpVars.n_cols();
            aInterpolationResult.set_size( 1, tNumInterpolationVars );
            for ( size_t i = 0; i < tNumInterpolationVars; i++ )
            {
                Matrix< Real_Matrix > tTmpVar = aInterpVars.get_column( i );

                aInterpolationResult( 0, i ) = ( tTmpVar( 0, 0 ) * ( 1 - xi ) + tTmpVar( 1, 0 ) * ( 1 + xi ) ) / 2;
            }
        }

        inline static Matrix< DDRMat >
        linear_interpolation_location( const Matrix< DDRMat > &aInterpVars,
                const Matrix< DDRMat >                        &aLocation )
        {
            moris::real      xi                    = aLocation( 0, 0 );
            size_t           tNumInterpolationVars = aInterpVars.n_cols();
            Matrix< DDRMat > tInterpolationResult( 1, tNumInterpolationVars );
            for ( size_t i = 0; i < tNumInterpolationVars; i++ )
            {
                tInterpolationResult( 0, i ) = ( aInterpVars( 0, i ) * ( 1 - xi ) + aInterpVars( 1, i ) * ( 1 + xi ) ) / 2;
            }

            return tInterpolationResult;
        }

        /**
         *Find a location given a value
         */
        template< typename Real, typename Real_Matrix >
        static void linear_interpolation_value( Matrix< Real_Matrix > const &aInterpVars,
                Real const                                                  &aValue,
                Matrix< Real_Matrix >                                       &aLocalCoordinate )
        {
            size_t                tNumPoints            = aInterpVars.n_rows();
            size_t                tNumInterpolationVars = aInterpVars.n_cols();
            Matrix< Real_Matrix > tTmpVar( tNumPoints, 1, 0 );
            for ( size_t i = 0; i < tNumInterpolationVars; i++ )
            {
                Matrix< Real_Matrix > tTmpVar = aInterpVars.get_column( i );

                aLocalCoordinate( 0, i ) = ( 2 * aValue - tTmpVar( 1, 0 ) - tTmpVar( 0, 0 ) ) / ( tTmpVar( 1, 0 ) - tTmpVar( 0, 0 ) );
            }
        }

        /**
         * Linear interpolation based on a local coordinate (aLclCoords) based on interpolation vars (aInterpVars)
         * Requires 1 local coordinate
         *
         * @param[in] aInterpVars - Interpolation Vars (x,y,z are treated as independent interpolation variables
         * @param[in] aLclCoords  - Local coordinates to interpolate to (a point at the center of edge has {{0}}
         */
        template< typename Matrix_View >
        static void bilinear_interpolation( const Matrix< DDRMat > &aInterpVars,
                Matrix_View                                        &aLocation,
                Matrix< DDRMat >                                   &aInterpolationResult )
        {
            // Get parametric coordinate
            moris::real xi      = aLocation( 0 );
            moris::real eta     = aLocation( 1 );
            size_t      tNumInt = aInterpVars.n_cols();
            aInterpolationResult.set_size( 1, tNumInt );
            for ( size_t i = 0; i < tNumInt; i++ )
            {
                auto tTmpVar = aInterpVars.get_column( i );

                aInterpolationResult( 0, i ) = ( tTmpVar( 0 ) * ( 1 - xi ) * ( 1 - eta )
                                                       + tTmpVar( 1 ) * ( 1 + xi ) * ( 1 - eta )
                                                       + tTmpVar( 2 ) * ( 1 + xi ) * ( 1 + eta )
                                                       + tTmpVar( 3 ) * ( 1 - xi ) * ( 1 + eta ) )
                                             / 4;
            }
        }

        /**
         * Linear interpolation based on a local coordinate (aLclCoords) based on interpolation vars (aInterpVars)
         * Requires 1 local coordinate
         *
         * @param[in] aInterpVars - Interpolation Vars (x,y,z are treated as independent interpolation variables
         * @param[in] aLclCoords  - Local coordinates to interpolate to (a point at the center of edge has {{0}}
         */
        template< typename Matrix_View >
        static void trilinear_interpolation( const Matrix< DDRMat > &aInterpVars,
                Matrix_View                                         &aLocation,
                Matrix< DDRMat >                                    &aInterpolationResult )
        {
            moris::real xi      = aLocation( 0 );
            moris::real eta     = aLocation( 1 );
            moris::real mu      = aLocation( 2 );
            size_t      tNumInt = aInterpVars.n_cols();
            aInterpolationResult.set_size( 1, tNumInt );

            for ( size_t i = 0; i < tNumInt; i++ )
            {
                auto tTmpVar = aInterpVars.get_column( i );

                aInterpolationResult( 0, i ) = ( tTmpVar( 0 ) * ( 1 - xi ) * ( 1 - eta ) * ( 1 - mu )
                                                       + tTmpVar( 1 ) * ( 1 + xi ) * ( 1 - eta ) * ( 1 - mu )
                                                       + tTmpVar( 2 ) * ( 1 + xi ) * ( 1 + eta ) * ( 1 + mu )
                                                       + tTmpVar( 3 ) * ( 1 - xi ) * ( 1 + eta ) * ( 1 - mu )
                                                       + tTmpVar( 4 ) * ( 1 - xi ) * ( 1 - eta ) * ( 1 + mu )
                                                       + tTmpVar( 5 ) * ( 1 + xi ) * ( 1 - eta ) * ( 1 + mu )
                                                       + tTmpVar( 6 ) * ( 1 + xi ) * ( 1 + eta ) * ( 1 + mu )
                                                       + tTmpVar( 7 ) * ( 1 - xi ) * ( 1 + eta ) * ( 1 + mu ) )
                                             / 8;
            }
        }
    };
}    // namespace moris::xtk

#endif /* SRC_TOOLS_CL_INTERPOLATON_HPP_ */
