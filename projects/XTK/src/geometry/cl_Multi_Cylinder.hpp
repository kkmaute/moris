/*
 * cl_Cylinder.hpp
 *
 *  Created on: Feb 2, 2018
 *      Author: ktdoble
 */

#ifndef SRC_GEOMETRY_CL_MULTI_CYLINDER_HPP_
#define SRC_GEOMETRY_CL_MULTI_CYLINDER_HPP_


#include <cmath>

#include "linalg/cl_XTK_Matrix.hpp"
#include "geometry/cl_Geometry.hpp"

// XTKL: Matrix Include

namespace xtk
{
class Multi_Cylinder : public Geometry
{
public:
    Multi_Cylinder(Cell<Cell<moris::real>> & mCenter,
                   Cell<moris::real>  & mRadius,
                   Cell<moris::real>  & mLength,
                   Cell<Cell<moris::real>> & mAxis   ) :
            mCenters(mCenter),
            mRadiuses(mRadius),
            mLengths(mLength),
              mAxes(mAxis)
    {
    }

    bool is_analytic() const
    {
        return true;
    }


    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 1;
    }

    moris::real evaluate_field_value_with_coordinate(moris::size_t const & aRowIndex,
                                              moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {
        moris::real lsVal = getSingleCylLSVal(mCenters(0),mAxes(0), mRadiuses(0), mLengths(0),aRowIndex, aCoordinates);

        for (size_t cInd = 1; cInd < mCenters.size(); ++cInd)
        {
            moris::real thisLsVal  =  getSingleCylLSVal(mCenters(cInd),mAxes(cInd), mRadiuses(cInd), mLengths(cInd), aRowIndex, aCoordinates);

            lsVal = std::min(thisLsVal, lsVal);
        }

        return lsVal;

    }


    moris::Matrix< moris::DDRMat > evaluate_sensitivity_dphi_dp_with_coordinate(moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {
        return moris::Matrix< moris::DDRMat >();
    }

private:

    Cell<Cell<moris::real>> mCenters;
    Cell<moris::real> mRadiuses;
    Cell<moris::real> mLengths;
    Cell<Cell<moris::real>> mAxes;

    moris::real getSingleCylLSVal(Cell<moris::real> const & aCenter,
                           Cell<moris::real> const & aAxis,
                           moris::real const & aRad,
                           moris::real const & aLength,
                           moris::size_t const & aRowIndex,
                           moris::Matrix< moris::DDRMat > const & aPointPosition) const
       {
           XTK_ASSERT(aCenter.size() == 3,"Centers need to have length 3");
           XTK_ASSERT(aAxis.size() == 3, "axis need to have length 3");
           XTK_ASSERT(aPointPosition.n_cols() == 3, "pointPosition need to have length 3");

           Cell<moris::real> relativePosition = {(aPointPosition(aRowIndex,0) - aCenter(0)),(aPointPosition(aRowIndex,1) - aCenter(1)),(aPointPosition(aRowIndex,2) - aCenter(2))};
           moris::real lsFromLeft = (relativePosition(0)*(-aAxis(0)) + relativePosition(1)*(-aAxis(1))+ relativePosition(2)*(-aAxis(2))) - aLength/2.0;
           moris::real lsFromRight = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2))) - aLength/2.0;

           moris::real axialCrd = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2)));
           Cell<moris::real> radDir = {(relativePosition(0) - aAxis(0)*axialCrd), (relativePosition(1) - aAxis(1)*axialCrd),(relativePosition(2) - aAxis(2)*axialCrd)};
           moris::real radDist = std::pow(radDir(0)*radDir(0)+radDir(1)*radDir(1)+radDir(2)*radDir(2), 0.5);
           moris::real lsFromRad = radDist - aRad;

           return std::max(std::max(lsFromLeft, lsFromRight), lsFromRad);
       }

};

}
#endif /* SRC_GEOMETRY_CL_MULTI_CYLINDER_HPP_ */
