/*
 * cl_Analytic_Level_Set_Plane.hpp
 *
 *  Created on: Feb 1, 2018
 *      Author: ktdoble
 */

#ifndef SRC_GEOMETRY_CL_PLANE_HPP_
#define SRC_GEOMETRY_CL_PLANE_HPP_

#include "cl_Matrix.hpp"
#include "cl_Geometry.hpp"

namespace xtk
{
class Plane : public Geometry
{
public:
    Plane(moris::real const & aXc,
          moris::real const & aYc,
          moris::real const & aZc,
          moris::real const & aXn,
          moris::real const & aYn,
          moris::real const & aZn):
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

    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 6;
        aNumCols = 3;
    }


    moris::real evaluate_field_value_with_coordinate(moris::size_t const & aRowIndex,
                                              moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {
        moris::real tDist = mXn*(aCoordinates(aRowIndex,0)-mXc) + mYn*(aCoordinates(aRowIndex,1)-mYc) + mZn*(aCoordinates(aRowIndex,2)-mZc);

        return tDist;
    }


    void
    get_plane_normal(moris::Matrix< moris::DDRMat > & aPlaneNormal)
    {
        aPlaneNormal(0,0) = mXn;
        aPlaneNormal(1,0) = mYn;
        aPlaneNormal(2,0) = mZn;
    }


private:

    moris::real mXc;
    moris::real mYc;
    moris::real mZc;
    moris::real mXn;
    moris::real mYn;
    moris::real mZn;

};
}


#endif /* SRC_GEOMETRY_CL_PLANE_HPP_ */
