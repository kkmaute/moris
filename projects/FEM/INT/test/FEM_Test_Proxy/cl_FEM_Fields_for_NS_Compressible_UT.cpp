/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Fields_for_NS_Compressible_UT.cpp
 *
 */

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_trans.hpp"

#include "cl_FEM_Field_Interpolator_Manager.hpp"

using namespace moris;

inline void
tConstValFunc(
        moris::Matrix< moris::DDRMat >&                aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*        aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

//------------------------------------------------------------------------------

inline void
fill_data_rectangle_element(
        moris::Matrix< moris::DDRMat >& tXHat,
        moris::Matrix< moris::DDRMat >& ttHat )
{
    // set values obtained from 1D Matlab code

    // element dimensions in 1D
    real tX1 = 0.0;
    real tX2 = 3.1;
    real tX3 = 0.5 * ( tX1 + tX2 );

    // chose y-size
    real tY1 = 0.0;
    real tY2 = 0.7;
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
    real tt1 = 1.3;
    real tt2 = 2.7;

    ttHat = { { tt1 }, { tt2 } };
}

//------------------------------------------------------------------------------

inline void
fill_smooth_UHat(
        moris::Matrix< moris::DDRMat >& tUHat,
        moris::uint                     aSpaceDim,
        moris::uint                     aInterpOrder )
{
    switch ( aSpaceDim )
    {
        case 2:
        {
            switch ( aInterpOrder )
            {
                case 1:
                    tUHat = {
                        { -3.988544000e+00, -5.796156800e+00 },
                        { -3.988544000e+00, -2.411949120e+01 },
                        { -4.126080000e-01, -1.262466240e+01 },
                        { -4.126080000e-01, -3.033833600e+00 },
                        { -4.534208000e+00, -6.478057600e+00 },
                        { -4.534208000e+00, -2.695707840e+01 },
                        { -4.690560000e-01, -1.410991680e+01 },
                        { -4.690560000e-01, -1.410991680e+01 }
                    };
                    break;
                case 2:
                    tUHat = {
                        { +1.750208351e+01, -2.377979856e+00 },
                        { +3.413302259e+01, +1.303485254e+01 },
                        { -4.562819428e+00, +2.274991354e+01 },
                        { -2.339635948e+00, -4.150322064e+00 },
                        { +2.969810550e+01, -3.302749800e+00 },
                        { +4.856563168e+00, +2.517867878e+01 },
                        { -3.969970500e+00, -5.764336200e+00 },
                        { +2.490256288e+00, -4.593407616e+00 },
                        { +4.225548000e+00, -6.379732800e+00 },
                        { +1.989650536e+01, -2.657742192e+00 },
                        { +3.880268692e+01, +1.456836461e+01 },
                        { -5.187048796e+00, +2.542637395e+01 },
                        { -2.659716436e+00, -4.638595248e+00 },
                        { +3.376103850e+01, -3.691308600e+00 },
                        { +5.520978976e+00, +2.814087629e+01 },
                        { -4.513093500e+00, -6.442493400e+00 },
                        { +2.830942816e+00, -5.133808512e+00 },
                        { +4.803636000e+00, -7.130289600e+00 }
                    };
                    break;
                case 3:
                    MORIS_ERROR( false, "Cubic not supported, yet." );
                    break;
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;
            }
            break;
        }
        case 3:
        {
            switch ( aInterpOrder )
            {
                case 1:
                    tUHat = {
                        { +7.039049472e+00, +1.325433024e+00, +4.582146800e+00 },
                        { +1.339690061e+01, -2.209055040e-01, -4.582146800e+00 },
                        { +5.249078016e+00, -9.738844800e-02, -1.908275600e+00 },
                        { +2.757990144e+00, +5.843306880e-01, +1.908275600e+00 },
                        { +2.963810304e+00, +1.590519629e+01, -1.374644040e+01 },
                        { +5.640800256e+00, -2.650866048e+00, +1.374644040e+01 },
                        { +2.210138112e+00, -1.168661376e+00, +5.724826800e+00 },
                        { +2.210138112e+00, -1.168661376e+00, +5.724826800e+00 },
                        { +7.986613824e+00, +1.570667328e+00, +5.273791600e+00 },
                        { +1.520032954e+01, -2.617778880e-01, -5.273791600e+00 },
                        { +5.955684672e+00, -1.154074560e-01, -2.196317200e+00 },
                        { +3.129258048e+00, +6.924447360e-01, +2.196317200e+00 },
                        { +3.362784768e+00, +1.884800794e+01, -1.582137480e+01 },
                        { +6.400138752e+00, -3.141334656e+00, +1.582137480e+01 },
                        { +2.507656704e+00, -1.384889472e+00, +6.588951600e+00 },
                        { +2.507656704e+00, -1.384889472e+00, +6.588951600e+00 }
                    };

                    break;
                case 2:
                    tUHat = {
                        { -5.581902643e+00, +1.334719685e-01, -2.362304806e+01 },
                        { -6.977378304e-01, -2.613042763e-01, -3.457031424e+00 },
                        { -1.445124096e-01, +1.001870856e-01, -7.271362560e-01 },
                        { -1.156099277e+00, -5.117469840e-02, -4.968764416e+00 },
                        { -2.124328827e+01, +3.757235913e+01, +2.362304806e+01 },
                        { -2.655411034e+00, -7.355715378e+01, +3.457031424e+00 },
                        { -5.499771264e-01, +2.820266460e+01, +7.271362560e-01 },
                        { -4.399817011e+00, -1.440567760e+01, +4.968764416e+00 },
                        { +1.395475661e+01, -1.253257920e-01, +2.880859520e+01 },
                        { -3.711811584e-01, -1.144995264e-02, -1.827243264e+00 },
                        { +2.890248192e+00, +4.805136000e-02, +6.059468800e+00 },
                        { -2.969449267e+00, +5.848536960e-03, -1.248616230e+01 },
                        { -1.141602808e+01, +1.481538850e+01, -1.624084554e+01 },
                        { -1.427003510e+00, -2.900477467e+01, -2.376709104e+00 },
                        { -2.955547296e-01, +1.112076650e+01, -4.999061760e-01 },
                        { -2.364437837e+00, -5.680391522e+00, -3.416025536e+00 },
                        { +5.310822067e+01, -3.527921045e+01, -2.880859520e+01 },
                        { -1.412620186e+00, -3.223161668e+00, +1.827243264e+00 },
                        { +1.099954253e+01, +1.352645784e+01, -6.059468800e+00 },
                        { -1.130096148e+01, +1.646363154e+00, +1.248616230e+01 },
                        { +1.518268877e+01, -6.095658240e-01, +1.046858120e+01 },
                        { +7.423623168e+00, -5.491584000e-03, +1.522702720e+01 },
                        { +2.825240371e+01, -1.545880896e+00, -1.522702720e+01 },
                        { -6.073075507e+00, +6.491876026e-01, -8.584236584e+00 },
                        { -7.591344384e-01, -1.270944743e+00, -1.256229744e+00 },
                        { +2.854007021e+01, -1.391116291e+01, +1.980590920e+01 },
                        { +5.911094592e+00, +5.333700960e+00, +4.165884800e+00 },
                        { -6.333312614e+00, +1.581672226e-01, -2.718879117e+01 },
                        { -7.916640768e-01, -3.096513230e-01, -3.978847488e+00 },
                        { -1.639660032e-01, +1.187239032e-01, -8.368926720e-01 },
                        { -1.311728026e+00, -6.064314480e-02, -5.718766592e+00 },
                        { -2.410296169e+01, +4.452407315e+01, +2.718879117e+01 },
                        { -3.012870211e+00, -8.716684744e+01, +3.978847488e+00 },
                        { -6.240125088e-01, +3.342077875e+01, +8.368926720e-01 },
                        { -4.992100070e+00, -1.707104526e+01, +5.718766592e+00 },
                        { +1.583328154e+01, -1.485138240e-01, +3.315706240e+01 },
                        { -4.211478528e-01, -1.356844608e-02, -2.103053568e+00 },
                        { +3.279320064e+00, +5.694192000e-02, +6.974105600e+00 },
                        { -3.369182822e+00, +6.930645120e-03, -1.437086605e+01 },
                        { -1.295280109e+01, +1.755656170e+01, -1.869229393e+01 },
                        { -1.619100137e+00, -3.437129686e+01, -2.735457648e+00 },
                        { -3.353409432e-01, +1.317835326e+01, -5.753637120e-01 },
                        { -2.682727546e+00, -6.731389073e+00, -3.931652032e+00 },
                        { +6.025740422e+01, -4.180664146e+01, -3.315706240e+01 },
                        { -1.602780595e+00, -3.819517572e+00, +2.103053568e+00 },
                        { +1.248025018e+01, +1.602915048e+01, -6.974105600e+00 },
                        { -1.282224476e+01, +1.950976601e+00, +1.437086605e+01 },
                        { +1.722651226e+01, -7.223489280e-01, +1.204874440e+01 },
                        { +8.422957056e+00, -6.507648000e-03, +1.752544640e+01 },
                        { +3.205561190e+01, -1.831902912e+00, -1.752544640e+01 },
                        { -6.890604902e+00, +7.693016083e-01, -9.879970408e+00 },
                        { -8.613256128e-01, -1.506097515e+00, -1.445849328e+00 },
                        { +3.238200274e+01, -1.648503446e+01, +2.279548040e+01 },
                        { +6.706818864e+00, +6.320553120e+00, +4.794697600e+00 }
                    };
                    break;
                case 3:
                    MORIS_ERROR( false, "Cubic not supported, yet." );
                    break;
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;
            }
            break;
        }
        default:
            MORIS_ERROR( false, "can only be 2 or 3D" );
            break;
    }
}

//------------------------------------------------------------------------------

inline void
fill_smooth_PHat(
        moris::Matrix< moris::DDRMat >& tPHat,
        moris::uint                     aSpaceDim,
        moris::uint                     aInterpOrder )
{
    switch ( aSpaceDim )
    {
        case 2:
        {
            switch ( aInterpOrder )
            {
                case 1:
                    tPHat = {
                        { +9.921836800e+01 },
                        { +8.827552000e+01 },
                        { +9.269344000e+01 },
                        { +9.951289600e+01 },
                        { +9.903289600e+01 },
                        { +8.549344000e+01 },
                        { +9.095968000e+01 },
                        { +9.095968000e+01 }
                    };
                    break;
                case 2:
                    tPHat = {
                        { +1.352743140e+02 },
                        { +1.223465549e+02 },
                        { +9.966897885e+01 },
                        { +9.947747901e+01 },
                        { +1.831070224e+02 },
                        { +1.083438348e+02 },
                        { +9.876892960e+01 },
                        { +1.131708466e+02 },
                        { +1.310307904e+02 },
                        { +1.436444901e+02 },
                        { +1.276491273e+02 },
                        { +9.959043146e+01 },
                        { +9.935349098e+01 },
                        { +2.028273328e+02 },
                        { +1.103237277e+02 },
                        { +9.847681120e+01 },
                        { +1.162961322e+02 },
                        { +1.383940288e+02 }
                    };
                    break;
                case 3:
                    MORIS_ERROR( false, "Cubic not supported, yet." );
                    break;
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;
            }
            break;
        }
        case 3:
        {
            switch ( aInterpOrder )
            {
                case 1:
                    tPHat = {
                        { +9.292708979e+01 },
                        { +9.135533197e+01 },
                        { +9.983689306e+01 },
                        { +9.986654886e+01 },
                        { +1.348204810e+02 },
                        { +1.425583657e+02 },
                        { +1.008029880e+02 },
                        { +1.008029880e+02 },
                        { +9.218812902e+01 },
                        { +9.045215770e+01 },
                        { +9.981985203e+01 },
                        { +9.985260621e+01 },
                        { +1.384584417e+02 },
                        { +1.470047621e+02 },
                        { +1.008868823e+02 },
                        { +1.008868823e+02 }
                    };
                    break;
                case 2:
                    tPHat = {
                        { +1.099386761e+02 },
                        { +1.015566601e+02 },
                        { +1.002714978e+02 },
                        { +1.017334091e+02 },
                        { +1.334134634e+02 },
                        { +1.052334340e+02 },
                        { +1.009127656e+02 },
                        { +1.058276575e+02 },
                        { +9.401284570e+01 },
                        { +1.007298094e+02 },
                        { +9.895577766e+01 },
                        { +1.046595521e+02 },
                        { +1.198316222e+02 },
                        { +1.031061577e+02 },
                        { +1.005417464e+02 },
                        { +1.034588423e+02 },
                        { +7.987140762e+01 },
                        { +1.024535922e+02 },
                        { +9.648936294e+01 },
                        { +1.156652427e+02 },
                        { +9.439899725e+01 },
                        { +9.719304090e+01 },
                        { +9.056310682e+01 },
                        { +1.092976646e+02 },
                        { +1.014562607e+02 },
                        { +8.805323965e+01 },
                        { +9.791636003e+01 },
                        { +1.109770453e+02 },
                        { +1.017192963e+02 },
                        { +1.002998632e+02 },
                        { +1.019145115e+02 },
                        { +1.369044222e+02 },
                        { +1.057802107e+02 },
                        { +1.010081292e+02 },
                        { +1.064365173e+02 },
                        { +9.338732211e+01 },
                        { +1.008060581e+02 },
                        { +9.884667981e+01 },
                        { +1.051463710e+02 },
                        { +1.219035827e+02 },
                        { +1.034306816e+02 },
                        { +1.005983468e+02 },
                        { +1.038202139e+02 },
                        { +7.776842035e+01 },
                        { +1.027099377e+02 },
                        { +9.612257997e+01 },
                        { +1.173019098e+02 },
                        { +9.381381786e+01 },
                        { +9.689977651e+01 },
                        { +8.957716275e+01 },
                        { +1.102690624e+02 },
                        { +1.016084074e+02 },
                        { +8.680507066e+01 },
                        { +9.769866630e+01 }
                    };
                    break;

                case 3:
                    MORIS_ERROR( false, "Cubic not supported, yet." );
                    break;
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;
            }
            break;
        }
        default:
            MORIS_ERROR( false, "can only be 2 or 3D" );
            break;
    }
    tPHat = trans( tPHat );
}

//------------------------------------------------------------------------------

inline void
fill_smooth_TempHat(
        moris::Matrix< moris::DDRMat >& tTempHat,
        moris::uint                     aSpaceDim,
        moris::uint                     aInterpOrder )
{
    switch ( aSpaceDim )
    {
        case 2:
        {
            switch ( aInterpOrder )
            {
                case 1:
                    tTempHat = {
                        { +1.130539528e+02 },
                        { +1.242430552e+02 },
                        { +1.165090744e+02 },
                        { +1.088895016e+02 },
                        { +1.154692216e+02 },
                        { +1.287285544e+02 },
                        { +1.195636168e+02 },
                        { +1.195636168e+02 }
                    };
                    break;
                case 2:
                    tTempHat = {
                        { +1.004433764e+02 },
                        { +1.037103604e+02 },
                        { +8.907103960e+01 },
                        { +9.869402360e+01 },
                        { +1.009334240e+02 },
                        { +9.151092640e+01 },
                        { +9.725057600e+01 },
                        { +9.898558240e+01 },
                        { +9.786438400e+01 },
                        { +1.005254108e+02 },
                        { +1.043968588e+02 },
                        { +8.704894120e+01 },
                        { +9.845238920e+01 },
                        { +1.011061280e+02 },
                        { +8.994026080e+01 },
                        { +9.674187200e+01 },
                        { +9.879789280e+01 },
                        { +9.746924800e+01 }
                    };
                    break;
                case 3:
                    MORIS_ERROR( false, "Cubic not supported, yet." );
                    break;
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;
            }
            break;
        }
        case 3:
        {
            switch ( aInterpOrder )
            {
                case 1:
                    tTempHat = {
                        { +1.223275670e+02 },
                        { +1.022532407e+02 },
                        { +1.041430555e+02 },
                        { +1.410539135e+02 },
                        { +1.094010808e+02 },
                        { +1.009487329e+02 },
                        { +1.017444444e+02 },
                        { +1.017444444e+02 },
                        { +1.266993983e+02 },
                        { +1.026944347e+02 },
                        { +1.049542831e+02 },
                        { +1.490924420e+02 },
                        { +1.112418519e+02 },
                        { +1.011344988e+02 },
                        { +1.020860140e+02 },
                        { +1.020860140e+02 }
                    };
                    break;
                case 2:
                    tTempHat = {
                        { +9.885286174e+01 },
                        { +8.388173533e+01 },
                        { +9.853673929e+01 },
                        { +9.989585961e+01 },
                        { +9.628812174e+01 },
                        { +4.784496473e+01 },
                        { +9.526522144e+01 },
                        { +9.966302541e+01 },
                        { +9.708354680e+01 },
                        { +9.310177093e+01 },
                        { +9.973523630e+01 },
                        { +9.950905245e+01 },
                        { +9.785261314e+01 },
                        { +6.982739480e+01 },
                        { +9.726084733e+01 },
                        { +9.980505427e+01 },
                        { +9.056302136e+01 },
                        { +7.767890107e+01 },
                        { +9.914328494e+01 },
                        { +9.841140550e+01 },
                        { +9.766348342e+01 },
                        { +9.875182828e+01 },
                        { +9.596120044e+01 },
                        { +9.908097014e+01 },
                        { +8.708685169e+01 },
                        { +9.454054188e+01 },
                        { +9.950437527e+01 },
                        { +9.862824726e+01 },
                        { +8.072571148e+01 },
                        { +9.825022670e+01 },
                        { +9.987546849e+01 },
                        { +9.556132040e+01 },
                        { +3.763279000e+01 },
                        { +9.433813193e+01 },
                        { +9.959704437e+01 },
                        { +9.651249303e+01 },
                        { +9.175106874e+01 },
                        { +9.968339446e+01 },
                        { +9.941292287e+01 },
                        { +9.743214578e+01 },
                        { +6.391947210e+01 },
                        { +9.672450974e+01 },
                        { +9.976688308e+01 },
                        { +8.871522135e+01 },
                        { +7.330833625e+01 },
                        { +9.897553654e+01 },
                        { +9.810035204e+01 },
                        { +9.720598367e+01 },
                        { +9.850743101e+01 },
                        { +9.517038653e+01 },
                        { +9.890102024e+01 },
                        { +8.455840307e+01 },
                        { +9.347155708e+01 },
                        { +9.940732987e+01 }
                    };
                    break;

                case 3:
                    MORIS_ERROR( false, "Cubic not supported, yet." );
                    break;
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;
            }
            break;
        }
        default:
            MORIS_ERROR( false, "can only be 2 or 3D" );
            break;
    }
    tTempHat = trans( tTempHat );
}

//------------------------------------------------------------------------------
