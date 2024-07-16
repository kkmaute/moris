/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Interpolation.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_INTERPOLATION_HPP_
#define PROJECTS_MTK_SRC_CL_INTERPOLATION_HPP_

// Matrix Include
#include "cl_Matrix.hpp"

#include "fn_print.hpp"
namespace moris
{
    namespace mtk
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
            static void
            linear_interpolation_location( const moris::Matrix< Real_Matrix >& aInterpVars,
                    const moris::Matrix< Real_Matrix >&                        aLocation,
                    moris::Matrix< Real_Matrix >&                              aInterpolationResult )
            {
                typename moris::Matrix< Real_Matrix >::Data_Type xi                    = aLocation( 0, 0 );
                size_t                                           tNumInterpolationVars = aInterpVars.n_cols();
                aInterpolationResult.set_size( 1, tNumInterpolationVars );
                for ( size_t i = 0; i < tNumInterpolationVars; i++ )
                {
                    moris::Matrix< Real_Matrix > tTmpVar = aInterpVars.get_column( i );

                    aInterpolationResult( 0, i ) = ( tTmpVar( 0, 0 ) * ( 1 - xi ) + tTmpVar( 1, 0 ) * ( 1 + xi ) ) / 2;
                }
            }

            inline static moris::Matrix< moris::DDRMat >
            linear_interpolation_location( const moris::Matrix< moris::DDRMat >& aInterpVars,
                    const moris::Matrix< moris::DDRMat >&                        aLocation )
            {
                moris::real                    xi                    = aLocation( 0, 0 );
                size_t                         tNumInterpolationVars = aInterpVars.n_cols();
                moris::Matrix< moris::DDRMat > tInterpolationResult( 1, tNumInterpolationVars );
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
            static void
            linear_interpolation_value( moris::Matrix< Real_Matrix > const & aInterpVars,
                    Real const &                                             aValue,
                    moris::Matrix< Real_Matrix >&                            aLocalCoordinate )
            {
                size_t                       tNumPoints            = aInterpVars.n_rows();
                size_t                       tNumInterpolationVars = aInterpVars.n_cols();
                moris::Matrix< Real_Matrix > tTmpVar( tNumPoints, 1, 0 );
                for ( size_t i = 0; i < tNumInterpolationVars; i++ )
                {
                    moris::Matrix< Real_Matrix > tTmpVar = aInterpVars.get_column( i );

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
            static void
            bilinear_interpolation( const moris::Matrix< moris::DDRMat >& aInterpVars,
                    Matrix_View&                                          aLocation,
                    moris::Matrix< moris::DDRMat >&                       aInterpolationResult )
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
            static void
            trilinear_interpolation( const moris::Matrix< moris::DDRMat >& aInterpVars,
                    Matrix_View&                                           aLocation,
                    moris::Matrix< moris::DDRMat >&                        aInterpolationResult )
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
                                                           + tTmpVar( 2 ) * ( 1 + xi ) * ( 1 + eta ) * ( 1 - mu )
                                                           + tTmpVar( 3 ) * ( 1 - xi ) * ( 1 + eta ) * ( 1 - mu )
                                                           + tTmpVar( 4 ) * ( 1 - xi ) * ( 1 - eta ) * ( 1 + mu )
                                                           + tTmpVar( 5 ) * ( 1 + xi ) * ( 1 - eta ) * ( 1 + mu )
                                                           + tTmpVar( 6 ) * ( 1 + xi ) * ( 1 + eta ) * ( 1 + mu )
                                                           + tTmpVar( 7 ) * ( 1 - xi ) * ( 1 + eta ) * ( 1 + mu ) )
                                                 / 8;
                }
            }

            /**
             * Linear interpolation based on a local coordinate (aLclCoords) based on interpolation vars (aInterpVars)
             * Requires 1 local coordinate
             *
             * @param[in] aInterpVars - Interpolation Vars (x,y,z are treated as independent interpolation variables
             * @param[in] aLocations  - Local coordinates to interpolate to (a point at the center of edge has {{0}} organized in row wise manner
             * @param[in] aInterpolationResult  - physical coordinates of the aLocations
             */
            template< typename Matrix_View >
            static void
            trilinear_interpolation_multivalue( const moris::Matrix< moris::DDRMat >& aInterpVars,
                    Matrix_View&                                                      aLocations,
                    moris::Matrix< moris::DDRMat >&                                   aInterpolationResult )
            {
                for ( size_t iRow = 0; iRow < aInterpolationResult.n_rows(); iRow++ )
                {
                    moris::real xi      = aLocations( iRow, 0 );
                    moris::real eta     = aLocations( iRow, 1 );
                    moris::real mu      = aLocations( iRow, 2 );
                    size_t      tNumInt = aInterpVars.n_cols();

                    for ( size_t i = 0; i < tNumInt; i++ )
                    {
                        auto tTmpVar = aInterpVars.get_column( i );

                        aInterpolationResult( iRow, i ) = ( tTmpVar( 0 ) * ( 1 - xi ) * ( 1 - eta ) * ( 1 - mu )
                                                                  + tTmpVar( 1 ) * ( 1 + xi ) * ( 1 - eta ) * ( 1 - mu )
                                                                  + tTmpVar( 2 ) * ( 1 + xi ) * ( 1 + eta ) * ( 1 - mu )
                                                                  + tTmpVar( 3 ) * ( 1 - xi ) * ( 1 + eta ) * ( 1 - mu )
                                                                  + tTmpVar( 4 ) * ( 1 - xi ) * ( 1 - eta ) * ( 1 + mu )
                                                                  + tTmpVar( 5 ) * ( 1 + xi ) * ( 1 - eta ) * ( 1 + mu )
                                                                  + tTmpVar( 6 ) * ( 1 + xi ) * ( 1 + eta ) * ( 1 + mu )
                                                                  + tTmpVar( 7 ) * ( 1 - xi ) * ( 1 + eta ) * ( 1 + mu ) )
                                                        / 8;
                    }
                }
            }
        };
    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_CL_INTERPOLATION_HPP_ */
