/*
 * cl_GEN_Multi_Geometry.hpp
 *
 *  Created on: Aug 9, 2019
 *      Author: doble
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_MULTI_GEOMETRY_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_MULTI_GEOMETRY_HPP_

#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

#include "cl_GEN_Geometry_Analytic.hpp"

namespace moris
{
namespace ge
{
class Multi_Geometry : public Geometry_Analytic
{
public:
    Multi_Geometry(){}

    Multi_Geometry(moris::Cell<GEN_Geometry*> const & aGeomVector) :
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
    moris::Cell<GEN_Geometry*> mGeometries;
};
}
}

#endif /* PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_MULTI_GEOMETRY_HPP_ */
