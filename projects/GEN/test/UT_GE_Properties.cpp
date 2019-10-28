/*
 * UT_GE_Properties.cpp
 *
 *  Created on: Sep 25, 2019
 *      Author: sonne
 */

#include "catch.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"

#include <vector>
//------------------------------------------------------------------------------
// GE includes
#include "cl_GE_Property.hpp"
#include "cl_GE_Geometry_Library.hpp"

// LINALG includes
#include "cl_Matrix.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"

using namespace moris;
using namespace ge;

TEST_CASE("property_test_01","[GE],[properties_value_test]")
{
    // specify the parameters for the circle function
    moris::Cell< moris::real > tCircleParams = {{0.6},{0},{0} };

    // create the GE property
    Property tGEProperty( circle_function,
                          circle_function_dphi_dp,
                          tCircleParams );

    // evaluate the property at specific coordinates and check values
    Matrix< DDRMat > tCoord = { {0, 0} };
    CHECK( tGEProperty.get_field_val_at_coordinate( tCoord ) == -tCircleParams(0) );
    Matrix< DDRMat > tCheckMat(3,2);
    tCheckMat(0,0) = 1.0;            tCheckMat(0,1) = 1.0;
    tCheckMat(1,0) = 1.0;            tCheckMat(1,1) = 0.0;
    tCheckMat(2,0) = 0.0;            tCheckMat(2,1) = 1.0;
    bool tMatrixMatch = all_true( tGEProperty.get_sensitivity_dphi_dp_at_coordinate( tCoord ) == tCheckMat );
    CHECK( tMatrixMatch );

}   // end functionals_test

