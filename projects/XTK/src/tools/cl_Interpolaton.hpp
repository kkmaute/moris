/*
 * cl_Interpolaton.hpp
 *
 *  Created on: Jun 21, 2017
 *      Author: ktdoble
 */

#ifndef SRC_TOOLS_CL_INTERPOLATON_HPP_
#define SRC_TOOLS_CL_INTERPOLATON_HPP_

// Matrix Include
#include "linalg/cl_XTK_Matrix.hpp"
namespace xtk
{
class Interpolation
{
public:
    //http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch16.d/IFEM.Ch16.pdf
    /**
     * Linear interpolation of a value given a location
     *
     * @param[in] aInterpVars - Interpolation Vars (x,y,z are treated as independent interpolation variables
     * @param[in] aLclCoords  - Local coordinates to interpolate to (a point at the center of edge has {{0}}
     */

    template<typename Real, typename Real_Matrix>
    static void linear_interpolation_location(const moris::Matrix<Real,Real_Matrix> & aInterpVars,
                                              const moris::Matrix<Real,Real_Matrix> & aLocation,
                                              moris::Matrix<Real,Real_Matrix> & aInterpolationResult)
    {
        Real xi = aLocation(0,0);
        size_t tNumInterpolationVars = aInterpVars.n_cols();
        aInterpolationResult.set_size(1, tNumInterpolationVars);
        for(size_t i = 0; i < tNumInterpolationVars; i++)
        {
            moris::Matrix<Real,Real_Matrix> tTmpVar = aInterpVars.get_column(i);

            aInterpolationResult(0, i) = (tTmpVar(0, 0) * (1 - xi) + tTmpVar(1, 0) * (1 + xi)) / 2;
        }
    }

    /**
     *Find a location given a value
     */
    template<typename Real, typename Real_Matrix>
    static void linear_interpolation_value(moris::Matrix<Real,Real_Matrix> const & aInterpVars,
                                           Real const & aValue,
                                           moris::Matrix<Real,Real_Matrix> & aLocalCoordinate )
    {
        size_t tNumPoints = aInterpVars.n_rows();
        size_t tNumInterpolationVars = aInterpVars.n_cols();
        moris::Matrix<Real,Real_Matrix> tTmpVar(tNumPoints, 1, 0);
        for(size_t i = 0; i < tNumInterpolationVars; i++)
        {
            moris::Matrix<Real,Real_Matrix> tTmpVar = aInterpVars.get_column(i);

            aLocalCoordinate(0, i) = (2 * aValue - tTmpVar(1, 0) - tTmpVar(0, 0)) / (tTmpVar(1, 0) - tTmpVar(0, 0));
        }
    }

    /**
     * Linear interpolation based on a local coordinate (aLclCoords) based on interpolation vars (aInterpVars)
     * Requires 1 local coordinate
     *
     * @param[in] aInterpVars - Interpolation Vars (x,y,z are treated as independent interpolation variables
     * @param[in] aLclCoords  - Local coordinates to interpolate to (a point at the center of edge has {{0}}
     */
    template<typename Real, typename Real_Matrix>
    static void bilinear_interpolation(const moris::Matrix<Real,Real_Matrix> & aInterpVars,
                                       const moris::Matrix<Real,Real_Matrix> & aLocation,
                                       moris::Matrix<Real,Real_Matrix> & aInterpolationResult)
    {
        Real xi = aLocation(0, 0);
        Real eta = aLocation(0, 1);
        size_t tNumInt = aInterpVars.n_cols();
        aInterpolationResult.set_size(1, tNumInt);
        for(size_t i = 0; i < tNumInt; i++)
        {
            moris::Matrix<Real,Real_Matrix> tTmpVar = aInterpVars.get_column(i);

            aInterpolationResult(0, i) = (tTmpVar(0, 0) * (1 - xi) * (1 - eta)
                                        + tTmpVar(1, 0) * (1 + xi) * (1 - eta)
                                        + tTmpVar(2, 0) * (1 + xi) * (1 + eta)
                                        + tTmpVar(3, 0) * (1 - xi) * (1 + eta)) / 4;
        }
    }

    /**
     * Linear interpolation based on a local coordinate (aLclCoords) based on interpolation vars (aInterpVars)
     * Requires 1 local coordinate
     *
     * @param[in] aInterpVars - Interpolation Vars (x,y,z are treated as independent interpolation variables
     * @param[in] aLclCoords  - Local coordinates to interpolate to (a point at the center of edge has {{0}}
     */
    template<typename Real, typename Real_Matrix>
    static void trilinear_interpolation(const moris::Matrix<Real,Real_Matrix> & aInterpVars,
                                        const moris::Matrix<Real,Real_Matrix> & aLocation,
                                        moris::Matrix<Real,Real_Matrix> & aInterpolationResult)
    {
        Real xi = aLocation(0, 0);
        Real eta = aLocation(0, 1);
        Real mu = aLocation(0, 2);
        size_t tNumInt = aInterpVars.n_cols();
        aInterpolationResult.set_size(1, tNumInt);

        for(size_t i = 0; i < tNumInt; i++)
        {
            moris::Matrix<Real,Real_Matrix> tTmpVar = aInterpVars.get_column(i);

            aInterpolationResult(0, i) = (tTmpVar(0, 0) * (1 - xi) * (1 - eta) * (1 - mu)
                                        + tTmpVar(1, 0) * (1 + xi) * (1 - eta) * (1 - mu)
                                        + tTmpVar(2, 0) * (1 + xi) * (1 + eta) * (1 + mu)
                                        + tTmpVar(3, 0) * (1 - xi) * (1 + eta) * (1 - mu)
                                        + tTmpVar(4, 0) * (1 - xi) * (1 - eta) * (1 + mu)
                                        + tTmpVar(5, 0) * (1 + xi) * (1 - eta) * (1 + mu)
                                        + tTmpVar(6, 0) * (1 + xi) * (1 + eta) * (1 + mu)
                                        + tTmpVar(7, 0) * (1 - xi) * (1 + eta) * (1 + mu))/ 8;
        }
    }



};
}


#endif /* SRC_TOOLS_CL_INTERPOLATON_HPP_ */
