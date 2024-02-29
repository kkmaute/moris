/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Interpolation.hpp
 *
 */

#pragma once

// Matrix Include
#include "cl_Matrix.hpp"

#include "fn_print.hpp"

namespace moris::gen
{
    class Interpolation
    {
        public:
            /**
             * 1D Linear interpolation of a value given a location; given values at -/+ 1
             *
             * @param[in] aInterpVars - 2xn matrix of interpolation values
             * @param[in] aLclCoords  - Local coordinates to interpolate to
             *
             * @param[out] aInterpolationResult - 1xn matrix of interpolated values
             */
            template<typename Real_Matrix>
            static void linear_interpolation_location(
                    const moris::Matrix< Real_Matrix > & aInterpVars,
                    const moris::Matrix< Real_Matrix > & aLocation,
                    moris::Matrix< Real_Matrix >       & aInterpolationResult)
            {
                // local coordinate
                typename moris::Matrix< Real_Matrix >::Data_Type xi = aLocation(0,0);

                // number of values to be interpolated
                size_t tNumInterpolationVars = aInterpVars.n_cols();

                // check that number of rows of aInterpVars is two
                MORIS_ASSERT( aInterpVars.n_rows() == 2,
                        "ERROR: linear_interpolation_location - Number of rows should be 2.\n");

                // set size of matrix to be returned
                aInterpolationResult.set_size(1, tNumInterpolationVars);

                // loop over all values to be interpolated
                for(size_t i = 0; i < tNumInterpolationVars; i++)
                {
                    moris::Matrix< Real_Matrix > tTmpVar = aInterpVars.get_column(i);

                    aInterpolationResult(0, i) = (tTmpVar(0, 0) * (1 - xi) + tTmpVar(1, 0) * (1 + xi)) / 2;
                }
            }

            //---------------------------------------------------------------------------
            /**
             * 1D Linear interpolation of a value given a location; given values at -/+ 1
             *
             * @param[in] aInterpVars - 2xn matrix of interpolation values
             * @param[in] aLclCoords  - Local coordinates to interpolate to
             *
             * @return  1xn matrix of interpolated value
             */
            inline
            static
            moris::Matrix< moris::DDRMat >
            linear_interpolation_location(
                    const moris::Matrix< moris::DDRMat > & aInterpVars,
                    const moris::Matrix< moris::DDRMat > & aLocation)
            {
                // local coordinate
                moris::real xi = aLocation(0,0);

                // number of values to be interpolated
                size_t tNumInterpolationVars = aInterpVars.n_cols();

                // check that number of rows of aInterpVars is two
                MORIS_ASSERT( aInterpVars.n_rows() == 2,
                        "ERROR: linear_interpolation_location - Number of rows should be 2.\n");

                // set size of matrix to be returned
                moris::Matrix< moris::DDRMat > tInterpolationResult(1, tNumInterpolationVars);

                // loop over all values to be interpolated
                for(size_t i = 0; i < tNumInterpolationVars; i++)
                {
                    tInterpolationResult(0, i) = (aInterpVars(0, i) * (1 - xi) + aInterpVars(1, i) * (1 + xi)) / 2;
                }

                return tInterpolationResult;
            }

            //---------------------------------------------------------------------------
            /**
             * 1D Linear interpolation: Find the local coordinate given a value
             *
             * @param[in] aInterpVars - nx2 matrix with nodal values
             * @param[in] value       - value for which location is to be found
             *
             * @param[out] aLocalCoordinate - nx1 matrix of local coordinate of points found
             */
            template<typename Real, typename Real_Matrix>
            static void linear_interpolation_value(
                    moris::Matrix< Real_Matrix > const & aInterpVars,
                    Real                         const & aValue,
                    moris::Matrix< Real_Matrix >       & aLocalCoordinate )
            {
                // number of values to be found
                size_t tNumInterpolationVars = aInterpVars.n_cols();

                // check that number of rows of aInterpVars is two
                MORIS_ASSERT( aInterpVars.n_rows() == 2,
                        "ERROR: linear_interpolation_location - Number of rows should be 2.\n");

                // loop over all points to be found
                for(size_t i = 0; i < tNumInterpolationVars; i++)
                {
                    moris::Matrix< Real_Matrix > tTmpVar = aInterpVars.get_column(i);

                    // check for denominator is not zero
                    const Real tDiffVar=tTmpVar(1) - tTmpVar(0);

                    MORIS_ASSERT( ! std::isnan(1.0/tDiffVar),
                            "ERROR: difference in nodal values for linear interpolation close to zero\n.");

                    // compute local coordinate
                    aLocalCoordinate(i) = (2.0 * aValue - tTmpVar(1) - tTmpVar(0)) / tDiffVar;
                }
            }

            //---------------------------------------------------------------------------
            /**
             * 1D Linear interpolation: Find the local coordinate given a value
             *
             * @param[in] aInterpVars - nx2 matrix with nodal values
             * @param[in] value       - value for which location is to be found
             *
             * @return nx1 matrix of local coordinate of points found
             */
            template<typename Real, typename Real_Matrix>
            static moris::Matrix< Real_Matrix >
            linear_interpolation_value(
                    moris::Matrix< Real_Matrix > const & aInterpVars,
                    Real                         const & aValue)
            {
                // create matrix to be returned
                Matrix<DDRMat> tLocalCoordinate(aInterpVars.n_cols(), 1);

                // determine local coordinates
                linear_interpolation_value(aInterpVars, aValue, tLocalCoordinate);

                return tLocalCoordinate;
            }

            //---------------------------------------------------------------------------
            /**
             * 2D linear interpolation based on a local coordinate based on nodal values
             * Requires 1 local coordinate
             *
             * @param[in] aInterpVars - 4xn matrix of nodal values given at (-1,-1),(1,-1),(1,1),(-1,1)
             * @param[in] aLclCoords  - 2x1 matrix of local coordinates to interpolate to
             *
             * @param[out] aInterpolationResult - 1xn matrix of interpolated values
             */
            template<typename Matrix_View>
            static void bilinear_interpolation(
                    const moris::Matrix< moris::DDRMat > & aInterpVars,
                    Matrix_View                          & aLocation,
                    moris::Matrix< moris::DDRMat >       & aInterpolationResult)
            {
                // parametric coordinates
                moris::real xi  = aLocation(0);
                moris::real eta = aLocation(1);

                // number of points to be interpolated
                size_t tNumInt = aInterpVars.n_cols();

                // check that number of rows of aInterpVars  four
                MORIS_ASSERT( aInterpVars.n_rows() == 4,
                        "ERROR: linear_interpolation_location - Number of rows should be 2.\n");

                // set size of matrix to be returned
                aInterpolationResult.set_size(1, tNumInt);

                // loop overall points to be interpolated
                for(size_t i = 0; i < tNumInt; i++)
                {
                    // extract column
                    auto tTmpVar = aInterpVars.get_column(i);

                    // interpolate value at point defined by local coordinate
                    aInterpolationResult(0, i) = (
                            tTmpVar(0) * (1 - xi) * (1 - eta) +
                            tTmpVar(1) * (1 + xi) * (1 - eta) +
                            tTmpVar(2) * (1 + xi) * (1 + eta) +
                            tTmpVar(3) * (1 - xi) * (1 + eta)) / 4.0;
                }
            }

            //---------------------------------------------------------------------------
            /**
             * 3D linear interpolation based on a local coordinate based on nodal values
             * Requires 1 local coordinate
             *
             * @param[in] aInterpVars - 4xn matrix of nodal values given at
             *                          (-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
             *                          (-1,-1, 1),(1,-1, 1),(1,1, 1),(-1,1, 1)
             *
             * @param[in] aLclCoords  - 3x1 matrix of local coordinates to interpolate to
             *
             * @param[out] aInterpolationResult - 1xn matrix of interpolated values
             */
            template<typename Matrix_View>
            static void trilinear_interpolation(
                    const moris::Matrix< moris::DDRMat > & aInterpVars,
                    Matrix_View                          & aLocation,
                    moris::Matrix< moris::DDRMat >       & aInterpolationResult)
            {
                // parametric coordinates
                moris::real xi  = aLocation(0);
                moris::real eta = aLocation(1);
                moris::real mu  = aLocation(2);

                // number of points to be interpolated
                size_t tNumInt = aInterpVars.n_cols();

                // set size of matrix to be returned
                aInterpolationResult.set_size(1, tNumInt);

                // loop overall points to be interpolated
                for(size_t i = 0; i < tNumInt; i++)
                {
                    // extract column
                    auto tTmpVar = aInterpVars.get_column(i);

                    // interpolate value at point defined by local coordinate
                    aInterpolationResult(0, i) = (
                                    tTmpVar(0) * (1 - xi) * (1 - eta) * (1 - mu) +
                                    tTmpVar(1) * (1 + xi) * (1 - eta) * (1 - mu) +
                                    tTmpVar(2) * (1 + xi) * (1 + eta) * (1 + mu) +
                                    tTmpVar(3) * (1 - xi) * (1 + eta) * (1 - mu) +
                                    tTmpVar(4) * (1 - xi) * (1 - eta) * (1 + mu) +
                                    tTmpVar(5) * (1 + xi) * (1 - eta) * (1 + mu) +
                                    tTmpVar(6) * (1 + xi) * (1 + eta) * (1 + mu) +
                                    tTmpVar(7) * (1 - xi) * (1 + eta) * (1 + mu))/8.0;
                }
            }
    };
}
