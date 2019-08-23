/*
 * cl_Multi_Geometry.hpp
 *
 *  Created on: Aug 9, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_GEOMETRY_CL_MULTI_GEOMETRY_HPP_
#define PROJECTS_XTK_SRC_GEOMETRY_CL_MULTI_GEOMETRY_HPP_

#include "cl_Matrix.hpp"
#include "cl_Geometry.hpp"
#include "cl_Cell.hpp"


namespace xtk
{
class Multi_Geometry : public Geometry
{
public:
    Multi_Geometry(){}

    Multi_Geometry(moris::Cell<Geometry*> const & aGeomVector) :
        mGeometries(aGeomVector)
    {
    }

    bool is_analytic() const
    {
        return true;
    }


    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 4;
    }

    moris::real evaluate_field_value_with_coordinate(moris::size_t const & aRowIndex,
                                              moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {
        moris::Matrix<moris::DDRMat> tValues(1,mGeometries.size());
        for(moris::uint i =0; i <mGeometries.size(); i++)
        {
                tValues(i) = mGeometries(i)->evaluate_field_value_with_coordinate(aRowIndex,aCoordinates);
        }


        return tValues.max();
    }


    moris::Matrix< moris::DDRMat >
    evaluate_sensitivity_dphi_dp_with_coordinate( moris::size_t const & aRowIndex,
                                                  moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {
        return moris::Matrix< moris::DDRMat >(0,0);
    }

private:
    moris::Cell<Geometry*> mGeometries;
};
}



#endif /* PROJECTS_XTK_SRC_GEOMETRY_CL_MULTI_GEOMETRY_HPP_ */
