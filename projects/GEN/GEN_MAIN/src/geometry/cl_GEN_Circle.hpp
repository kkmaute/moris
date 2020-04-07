/*
 * cl_GEN_Circle.hpp
 *
 *  Created on: Aug 7, 2019
 *      Author: ryan
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_CIRCLE_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_CIRCLE_HPP_

#include <cmath>

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

#include "cl_GEN_Enums.hpp"
#include "cl_GEN_Geometry.hpp"

namespace moris
{
namespace ge
{

class Circle : public GEN_Geometry
{
public:
    Circle()
    {

    }

    Circle(moris::real const & aRadius,
           moris::real const & aXCenter,
           moris::real const & aYCenter) :
            mRadius(aRadius),
            mXCenter(aXCenter),
            mYCenter(aYCenter)
    {

    }
    //------------------------------------------------------------------------------

    bool is_analytic() const
    {
        return true;
    }
    //------------------------------------------------------------------------------

    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 4;
    }
    //------------------------------------------------------------------------------

    moris::real evaluate_field_value_with_coordinate( moris::size_t const & aRowIndex,
                                                      moris::Matrix< moris::DDRMat > const & aCoordinates ) const
    {

        /*
         * the below form linearizes the circle equation
         */
        Matrix< DDRMat > tCoord( 1,2 );
        tCoord(0) = aCoordinates(aRowIndex, 0); tCoord(1) = aCoordinates(aRowIndex, 1);
        Matrix< DDRMat > tCenter( 1,2 );
        tCenter(0) = mXCenter; tCenter(1) = mYCenter;

        moris::real tFunctionValue = norm( tCoord-tCenter ) - mRadius;
        return tFunctionValue;
    }

    moris::real evaluate_field_value_with_single_coordinate(moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {
        return this->evaluate_field_value_with_coordinate(0,aCoordinates);
    }
    //------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat > evaluate_sensitivity_dphi_dp_with_coordinate( moris::size_t                  const & aRowIndex,
                                                                                 moris::Matrix< moris::DDRMat > const & aCoordinates ) const
    {
        moris::Matrix< moris::DDRMat > tSensitivityDxDp(3, 2, 0.0);

        moris::real sign = 0.0;

        // dx/dr = r/sqrt(r^2 - (y-yc)^2)
//        moris::real tSqrt = mRadius * mRadius - (aCoordinates(aRowIndex, 1) - mYCenter) * (aCoordinates(aRowIndex, 1) - mYCenter);
        moris::real tSqrt = mRadius * mRadius - std::pow((aCoordinates(aRowIndex, 1) - mYCenter),2);

        if(tSqrt < 0.0)
        {
            sign = -1.0;
        }
        else if(tSqrt > 0.0)
        {
            sign = 1.0;
        }
        else
        {
            std::cout << "zero denominator detected";
        }

        (tSensitivityDxDp)(0, 0) = sign * mRadius / std::sqrt(std::abs(tSqrt));

        //dy/dr = r/sqrt(r^2 - (x-xc)^2)
//        tSqrt = mRadius * mRadius - (aCoordinates(aRowIndex, 0) - mXCenter) * (aCoordinates(aRowIndex, 0) - mXCenter);
        tSqrt = mRadius * mRadius - std::pow((aCoordinates(aRowIndex, 0) - mXCenter),2);

        if(tSqrt < 0.0)
        {
            sign = -1.0;
        }
        else if(tSqrt > 0.0)
        {
            sign = 1.0;
        }
        else
            std::cout << "zero denominator detected";

        tSensitivityDxDp(0, 1) = mRadius / std::sqrt(std::abs(tSqrt));

        // fill remaining values in tSensitivity
        tSensitivityDxDp(1,0) = 1.0; // dx/dxc
        tSensitivityDxDp(1,1) = 0.0; // dy/dxc
        tSensitivityDxDp(2,0) = 0.0; // dx/dyc
        tSensitivityDxDp(2,1) = 1.0; // dy/dyc

        return tSensitivityDxDp;
    }
    //------------------------------------------------------------------------------
//    void get_pdv_values( const enum GEN_PDV aPdvType,
//                         Matrix< DDRMat > & aPdvValueMatrix )
//    {
//        switch(aPdvType)
//        {
//        case(GEN_PDV::RADIUS)   :
//               {
//                   aPdvValueMatrix.resize( 1,1 );
//                   aPdvValueMatrix( 0, 0 ) = mRadius;
//                   break;
//               }
//        case(GEN_PDV::XCOORD)   :
//               {
//                   aPdvValueMatrix.resize( 1,1 );
//                   aPdvValueMatrix( 0, 0 ) = mXCoord;
//                   break;
//               }
//        case(GEN_PDV::YCOORD)   :
//               {
//                   aPdvValueMatrix.resize( 1,1 );
//                   aPdvValueMatrix( 0, 0 ) = mYCoord;
//               }
//        default   :
//               {
//                   MORIS_ERROR( false, "Circle::get_pdv_values() - requested pdv type does not exist for this geometry " );
//                   break;
//               }
//        }
//    }

private:
    moris::real mRadius;
    moris::real mXCenter;
    moris::real mYCenter;
};

}
}

#endif /* PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_CIRCLE_HPP_ */


