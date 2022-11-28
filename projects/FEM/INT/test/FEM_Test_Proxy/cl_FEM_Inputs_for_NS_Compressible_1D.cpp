/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Inputs_for_NS_Compressible_1D.cpp
 *
 */

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_trans.hpp"

#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_FEM_Convert_Dimensions.cpp"

using namespace moris;

inline void
tConstValFunc(
        moris::Matrix< moris::DDRMat >&                aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*        aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

// heat load distribution
inline void
Func_Heat_Load_Distribution(
        moris::Matrix< moris::DDRMat >&                aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*        aFIManager )
{
    // from Heated Channel Example
    real tQfac          = 200.0; /* exponential decay factor for heat load distribution */
    real tChannelLength = 100.0; /* channel lenght from example */

    // get x-coordinate
    real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );

    // get element size in x-direction
    real tHx = 1.0;

    // get element number in x-direction
    real tElemNumber = std::floor( tX / tHx ) + 1.0;

    // get element left and right node positions
    real tLeftNodePos  = ( tElemNumber - 1.0 ) * tHx;
    real tRightNodePos = tElemNumber * tHx;

    // compute heat loads at left and right nodes
    real tQl = std::exp( -1.0 * tQfac * std::pow( 2.0 * tLeftNodePos / tChannelLength - 1.0, 2.0 ) );
    real tQr = std::exp( -1.0 * tQfac * std::pow( 2.0 * tRightNodePos / tChannelLength - 1.0, 2.0 ) );

    // compute linear interpolation
    real tAlpha = ( ( tX - tLeftNodePos ) / ( tRightNodePos - tLeftNodePos ) );
    real tQ     = tQl + tAlpha * ( tQr - tQl );
    // real tQ = std::exp( -1.0 * tQfac * std::pow( 2.0 * tX / tChannelLength - 1.0, 2.0 ) );

    // return value
    aPropMatrix = tQ * aParameters( 0 );
}

inline void
fill_data(
        moris::Matrix< moris::DDRMat >& tXHat,
        moris::Matrix< moris::DDRMat >& ttHat,
        moris::Matrix< moris::DDRMat >& tPHat,
        moris::Matrix< moris::DDRMat >& tUHat,
        moris::Matrix< moris::DDRMat >& tTHat )
{
    // set values obtained from 1D Matlab code

    // element dimensions in 1D
    real tX1 = 5.000000e+01;
    real tX2 = 5.100000e+01;
    real tX3 = 0.5 * ( tX1 + tX2 );
    // real tt1 = 3.019260e-02;
    // real tt2 = 6.038520e-02;
    real tt1 = 0.0;
    real tt2 = 3.019260e-02;

    // chose y-size
    real tY1 = 0.000000e+00;
    real tY2 = 1.000000e+00;
    real tY3 = 0.5 * ( tY1 + tY2 );

    tXHat = {
        { tX1, tY1 },
        { tX2, tY1 },
        { tX2, tY2 },
        { tX1, tY2 },
        { tX3, tY1 },
        { tX2, tY3 },
        { tX3, tY2 },
        { tX1, tY3 },
        { tX3, tY3 }
    };

    // time frame
    ttHat = { { tt1 }, { tt2 } };

    // Pressure DoF values in 1D
    Matrix< DDRMat > tPhat1D = {
        { +7.636683400005384e+04 },
        { +7.936683400004805e+04 },
        { +8.836683400005080e+04 },
        { +6.679592850378335e+04 },
        { +7.978749797591206e+04 },
        { +9.879333215766176e+04 }
    };
    // { +7.8366834e+04 },
    // { +7.8366834e+04 },
    // { +7.8366834e+04 },
    // { +7.8366834e+04 },
    // { +7.8366834e+04 },
    // { +7.8366834e+04 } };
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 } };

    // sort pressure DoFs into 2D element
    tPHat = fem::convert_DoF_vector_1D_to_2D_quadratic( tPhat1D );

    // Velocity DoF values in 1D
    Matrix< DDRMat > tUXhat1D = {
        { -1.298297133102633e-06 },
        { +3.792832711608647e-06 },
        { +7.482483137997108e-06 },
        { -2.219694454936193e-01 },
        { +9.713214928148250e-01 },
        { +3.404845514857839e-01 }
    };
    // { +0.0e+00 },
    // { +0.0e+00 },
    // { +0.0e+00 },
    // { +0.0e+00 },
    // { +0.0e+00 },
    // { +0.0e+00 } };
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 } };

    // sort velocity DoFs into 2D element
    tUHat.set_size( 18, 2, 0.0 );
    tUHat( { 0, 17 }, { 0, 0 } ) = fem::convert_DoF_vector_1D_to_2D_quadratic( tUXhat1D ).matrix_data();

    // Temperature DoF values in 1D
    Matrix< DDRMat > tThat1D = {
        { +3.430147088092934e+02 },
        { +2.830131125499078e+02 },
        { +5.630138803472934e+02 },
        { +7.564604507806324e+02 },
        { +2.962188934473666e+02 },
        { +3.763403645201755e+02 }
    };
    // { +2.73e+02 },
    // { +2.73e+02 },
    // { +2.73e+02 },
    // { +2.73e+02 },
    // { +2.73e+02 },
    // { +2.73e+02 } };
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 },
    // { +1.0e+00 } };

    // sort temperature DoFs into 2D element
    tTHat = fem::convert_DoF_vector_1D_to_2D_quadratic( tThat1D );
}
