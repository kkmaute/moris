/*
 * cl_Analytic_Level_Set_Gyroid.hpp
 *
 *  Created on: Aug 15, 2017
 *      Author: ktdoble
 */

#ifndef SRC_GEOMETRY_CL_GYROID_HPP_
#define SRC_GEOMETRY_CL_GYROID_HPP_


#include "linalg/cl_XTK_Matrix.hpp"
#include "geometry/cl_Geometry.hpp"

// XTKL: Matrix Include

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Gyroid : public Geometry<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:
    Gyroid()
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

        Real tFuncVal = std::sin(aCoordinates(aRowIndex,0))*std::cos(aCoordinates(aRowIndex,1))+
                        std::sin(aCoordinates(aRowIndex,1))*std::cos(aCoordinates(aRowIndex,2))+
                        std::sin(aCoordinates(aRowIndex,2))*std::cos(aCoordinates(aRowIndex,0));

        return tFuncVal;
    }

private:
};
}


#endif /* SRC_GEOMETRY_CL_GYROID_HPP_ */
