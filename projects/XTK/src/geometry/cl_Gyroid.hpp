/*
 * cl_Analytic_Level_Set_Gyroid.hpp
 *
 *  Created on: Aug 15, 2017
 *      Author: ktdoble
 */

#ifndef SRC_GEOMETRY_CL_GYROID_HPP_
#define SRC_GEOMETRY_CL_GYROID_HPP_


#include "cl_Matrix.hpp"
#include "cl_Geometry.hpp"

// XTKL: Matrix Include

namespace xtk
{
class Gyroid : public Geometry
{
public:
    Gyroid()
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

        moris::real tFuncVal = std::sin(aCoordinates(aRowIndex,0))*std::cos(aCoordinates(aRowIndex,1))+
                        std::sin(aCoordinates(aRowIndex,1))*std::cos(aCoordinates(aRowIndex,2))+
                        std::sin(aCoordinates(aRowIndex,2))*std::cos(aCoordinates(aRowIndex,0));
        return tFuncVal;
    }

private:
};
}


#endif /* SRC_GEOMETRY_CL_GYROID_HPP_ */
