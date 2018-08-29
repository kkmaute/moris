/*
 * cl_Analytic_Level_Set_Plane.hpp
 *
 *  Created on: Feb 1, 2018
 *      Author: ktdoble
 */

#ifndef SRC_GEOMETRY_CL_PLANE_HPP_
#define SRC_GEOMETRY_CL_PLANE_HPP_

#include "linalg/cl_XTK_Matrix.hpp"
#include "geometry/cl_Geometry.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Plane : public Geometry<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:
    Plane(Real const & aXc,
          Real const & aYc,
          Real const & aZc,
          Real const & aXn,
          Real const & aYn,
          Real const & aZn):
              mXc(aXc),
              mYc(aYc),
              mZc(aZc),
              mXn(aXn),
              mYn(aYn),
              mZn(aZn)
    {
    }

    bool is_analytic() const
    {
        return true;
    }

    void get_dphi_dp_size(Integer & aNumRows, Integer & aNumCols) const
    {
        aNumRows = 6;
        aNumCols = 3;
    }


    Real evaluate_field_value_with_coordinate(Integer const & aRowIndex,
                                              moris::Matrix<Real,Real_Matrix> const & aCoordinates) const
    {
        Real tDist = mXn*(aCoordinates(aRowIndex,0)-mXc) + mYn*(aCoordinates(aRowIndex,1)-mYc) + mZn*(aCoordinates(aRowIndex,2)-mZc);

        return tDist;
    }


    void
    get_plane_normal(moris::Matrix<Real, Real_Matrix> & aPlaneNormal)
    {
        aPlaneNormal(0,0) = mXn;
        aPlaneNormal(1,0) = mYn;
        aPlaneNormal(2,0) = mZn;
    }


private:

    Real mXc;
    Real mYc;
    Real mZc;
    Real mXn;
    Real mYn;
    Real mZn;

};
}


#endif /* SRC_GEOMETRY_CL_PLANE_HPP_ */
