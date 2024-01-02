/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Fields_for_NS_Compressible_UT_for_BD.cpp
 *
 */

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_trans.hpp"

#include "cl_FEM_Field_Interpolator_Manager.hpp"

using namespace moris;

//------------------------------------------------------------------------------

inline void
tConstValFunc(
        moris::Matrix< moris::DDRMat >&                aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*        aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

//------------------------------------------------------------------------------

// heat load distribution
inline void
Func_Heat_Load_Distribution(
        moris::Matrix< moris::DDRMat >&                aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*        aFIManager )
{
    // get x-coordinate
    real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );

    // constants from Heated Channel EXA
    real tQfac          = 200.0;
    real tChannelLength = 100.0;

    // get Heat Load
    real tQ = std::exp( -1.0 * tQfac * std::pow( 2.0 * tX / tChannelLength - 1.0, 2.0 ) );

    // return value
    aPropMatrix = tQ * aParameters( 0 );
}

//------------------------------------------------------------------------------

inline void
fill_data_Nitsche_element(
        moris::Matrix< moris::DDRMat >& tXHat,
        moris::Matrix< moris::DDRMat >& ttHat )
{
    // set values obtained from 1D Matlab code

    // element dimensions in 1D
    real tX1 = 0.0;
    real tX2 = 25.0;
    real tX3 = 0.5 * ( tX1 + tX2 );

    // chose y-size
    real tY1 = 0.0;
    real tY2 = 1.0;
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
    real tt1 = 3.0 * 0.0301926;
    real tt2 = 4.0 * 0.0301926;

    ttHat = { { tt1 }, { tt2 } };
}

//------------------------------------------------------------------------------

inline void
fill_Nitsche_UHat( moris::Matrix< moris::DDRMat >& tUHat )
{
    tUHat = {
        { +4.767499875461087e-03, +2.385350017824317e-09 },
        { +1.287701749849472e-01, +2.497634827421489e-09 },
        { +1.287701749451793e-01, +2.507327946509253e-09 },
        { +4.767499854282951e-03, +2.389028999797776e-09 },
        { -1.379104875687944e-01, -9.507824431176954e-10 },
        { +1.287701749651091e-01, +2.503477661075256e-09 },
        { -1.379104875414982e-01, -9.373454778273539e-10 },
        { +4.767499864861142e-03, +2.388090998103339e-09 },
        { -1.379104875551917e-01, -9.444935425401964e-10 },
        { -1.676507218088517e-01, +4.100444278400517e-09 },
        { -4.527476028069115e-01, +4.893677585820106e-09 },
        { -4.527476027851524e-01, +4.901286216323920e-09 },
        { -1.676507218013225e-01, +4.103891139367800e-09 },
        { -3.539221865043858e-01, -1.524453663194401e-09 },
        { -4.527476027959706e-01, +4.902055616175119e-09 },
        { -3.539221864543057e-01, -1.501238243891667e-09 },
        { -1.676507218052863e-01, +4.105855501378721e-09 },
        { -3.539221864793631e-01, -1.514264650866220e-09 }
    };
}

//------------------------------------------------------------------------------

inline void
fill_Nitsche_PHat( moris::Matrix< moris::DDRMat >& tPHat )
{
    tPHat = {
        { +7.871934902150265e+04 },
        { +7.881065319148856e+04 },
        { +7.881065319141212e+04 },
        { +7.871934902143817e+04 },
        { +7.822749630250456e+04 },
        { +7.881065319145036e+04 },
        { +7.822749630252989e+04 },
        { +7.871934902147045e+04 },
        { +7.822749630251728e+04 },
        { +7.892737538218674e+04 },
        { +7.919005331973224e+04 },
        { +7.919005331965104e+04 },
        { +7.892737538213342e+04 },
        { +7.817872616941837e+04 },
        { +7.919005331969161e+04 },
        { +7.817872616943391e+04 },
        { +7.892737538216007e+04 },
        { +7.817872616942617e+04 }
    };
}

//------------------------------------------------------------------------------

inline void
fill_Nitsche_THat( moris::Matrix< moris::DDRMat >& tTHat )
{
    tTHat = {
        { +2.732407703143502e+02 },
        { +2.733200232559518e+02 },
        { +2.733200232558551e+02 },
        { +2.732407703142323e+02 },
        { +2.728302777124552e+02 },
        { +2.733200232559036e+02 },
        { +2.728302777125072e+02 },
        { +2.732407703142914e+02 },
        { +2.728302777124812e+02 },
        { +2.734958407004418e+02 },
        { +2.738413345274006e+02 },
        { +2.738413345273013e+02 },
        { +2.734958407003357e+02 },
        { +2.727572086513168e+02 },
        { +2.738413345273511e+02 },
        { +2.727572086513589e+02 },
        { +2.734958407003888e+02 },
        { +2.727572086513379e+02 }
    };
}

//------------------------------------------------------------------------------
