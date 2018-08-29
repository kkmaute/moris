/*
 * cl_Analytic_Plane.hpp
 *
 *  Created on: Jan 23, 2018
 *      Author: doble
 */

#ifndef SRC_GEOMETRY_CL_ANALYTIC_LEVEL_SET_PLANE_HPP_
#define SRC_GEOMETRY_CL_ANALYTIC_LEVEL_SET_PLANE_HPP_



#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "geometry/cl_Geometry.hpp"
#include "assert/fn_xtk_assert.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Plane : public Geometry<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:
    Plane(Real const & aXo,
          Real const & aYo,
          Real const & aZo,
          Real const & aA ,
          Real const & aB ,
          Real const & aC
          ) :
               mXo( aXo),
               mYo( aYo),
               mZo( aZo),
               mA ( aA ),
               mB ( aB ),
               mC ( aC )
    {
    }

    bool is_analytic() const
    {
        return true;
    }


    void get_dphi_dp_size(Integer & aNumRows, Integer & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 3;
    }

    Real evaluate_field_value_with_coordinate(moris::Mat_New<Real, Real_Matrix> const & aCoordinates) const
    {
        return mA*(aCoordinates(0,0)-mXo) + mB*(aCoordinates(0,1)-mYo) + mC*(aCoordinates(0,2) - mZo);
    }


    void  get_plane_normal(moris::Mat_New<Real, Real_Matrix> & aPlaneNormal)
    {
        aPlaneNormal.resize(3,1);

        Real tLen = std::pow(std::pow(mA,2) + std::pow(mB,2) + std::pow(mC,2),0.5);

        aPlaneNormal(0,0) = mA/tLen;
        aPlaneNormal(1,0) = mB/tLen;
        aPlaneNormal(2,0) = mC/tLen;
    }

    std::shared_ptr<Matrix_Base<Real, Real_Matrix>> evaluate_sensitivity_dphi_dp_with_coordinate(moris::Mat_New<Real, Real_Matrix> const & aCoordinates) const
    {
        return NULL;
    }

private:
    Real const mXo;
    Real const mYo;
    Real const mZo;
    Real const mA ;
    Real const mB ;
    Real const mC ;



};
}



#endif /* SRC_GEOMETRY_CL_ANALYTIC_LEVEL_SET_PLANE_HPP_ */
