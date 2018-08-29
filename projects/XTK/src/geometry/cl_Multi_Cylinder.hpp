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
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Multi_Cylinder : public Geometry<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:
    Multi_Cylinder(Cell<Cell<Real>> & mCenter,
                   Cell<Real>  & mRadius,
                   Cell<Real>  & mLength,
                   Cell<Cell<Real>> & mAxis   ) :
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


    void get_dphi_dp_size(Integer & aNumRows, Integer & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 1;
    }

    Real evaluate_field_value_with_coordinate(Integer const & aRowIndex,
                                              moris::Mat_New<Real,Real_Matrix> const & aCoordinates) const
    {
        Real lsVal = getSingleCylLSVal(mCenters(0),mAxes(0), mRadiuses(0), mLengths(0),aRowIndex, aCoordinates);

        for (size_t cInd = 1; cInd < mCenters.size(); ++cInd)
        {
            Real thisLsVal  =  getSingleCylLSVal(mCenters(cInd),mAxes(cInd), mRadiuses(cInd), mLengths(cInd), aRowIndex, aCoordinates);

            lsVal = std::min(thisLsVal, lsVal);
        }

        return lsVal;

    }


    moris::Mat_New<Real,Real_Matrix> evaluate_sensitivity_dphi_dp_with_coordinate(moris::Mat_New<Real, Real_Matrix> const & aCoordinates) const
    {
        return moris::Mat_New<Real,Real_Matrix>();
    }

private:

    Cell<Cell<Real>> mCenters;
    Cell<Real> mRadiuses;
    Cell<Real> mLengths;
    Cell<Cell<Real>> mAxes;

    Real getSingleCylLSVal(Cell<Real> const & aCenter,
                           Cell<Real> const & aAxis,
                           Real const & aRad,
                           Real const & aLength,
                           Integer const & aRowIndex,
                           moris::Mat_New<Real,Real_Matrix> const & aPointPosition) const
       {
           XTK_ASSERT(aCenter.size() == 3,"Centers need to have length 3");
           XTK_ASSERT(aAxis.size() == 3, "axis need to have length 3");
           XTK_ASSERT(aPointPosition.n_cols() == 3, "pointPosition need to have length 3");

           Cell<Real> relativePosition = {(aPointPosition(aRowIndex,0) - aCenter(0)),(aPointPosition(aRowIndex,1) - aCenter(1)),(aPointPosition(aRowIndex,2) - aCenter(2))};
           Real lsFromLeft = (relativePosition(0)*(-aAxis(0)) + relativePosition(1)*(-aAxis(1))+ relativePosition(2)*(-aAxis(2))) - aLength/2.0;
           Real lsFromRight = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2))) - aLength/2.0;

           Real axialCrd = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2)));
           Cell<Real> radDir = {(relativePosition(0) - aAxis(0)*axialCrd), (relativePosition(1) - aAxis(1)*axialCrd),(relativePosition(2) - aAxis(2)*axialCrd)};
           Real radDist = std::pow(radDir(0)*radDir(0)+radDir(1)*radDir(1)+radDir(2)*radDir(2), 0.5);
           Real lsFromRad = radDist - aRad;

           return std::max(std::max(lsFromLeft, lsFromRight), lsFromRad);
       }

};

}
#endif /* SRC_GEOMETRY_CL_MULTI_CYLINDER_HPP_ */
