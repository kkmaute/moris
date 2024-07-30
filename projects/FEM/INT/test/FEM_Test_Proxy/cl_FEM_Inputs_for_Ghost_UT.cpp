/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Inputs_for_Diffusion_UT.cpp
 *
 */

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_trans.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

using namespace moris;

//------------------------------------------------------------------------------

/**
 * @brief Constant value function for unit property
 */
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
fill_X_hat(
        moris::uint                      aSpaceDim,
        moris::Matrix< moris::DDRMat > & aXHat,
        moris::mtk::Geometry_Type &      aGeometryType )
{
    switch( aSpaceDim )
    {
        // QUAD element
        case 2 :
        {
            // set geometry type
            aGeometryType = mtk::Geometry_Type::QUAD;

            // fill space coeff xHat
            aXHat = {
                    { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 }};
            break;
        }

        // HEX element
        case 3 :
        {
            // set geometry type
            aGeometryType = mtk::Geometry_Type::HEX;

            // fill space coeff xHat
            aXHat = {
                    { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 }};
            break;
        }
        default:
        {
            MORIS_ERROR( false, "fill_X_Hat() - QUAD or HEX only." );
            break;
        }
    }
}

//------------------------------------------------------------------------------

inline void
fill_u_hat_leader(
        moris::Matrix< moris::DDRMat >& aUHat,
        moris::uint                     aSpaceDim,
        moris::uint                     aInterpOrder )
{
    // initialize the vector for time-wise constant coefficient vector
    moris::Matrix< moris::DDRMat > tUHat;

    // switch between 2D and 3D
    switch ( aSpaceDim )
    {
        // QUAD elements
        case 2:
        {
            switch ( aInterpOrder )
            {
                // QUAD 4
                case 1:
                    tUHat = {
                            { 1.891288256070320e-01 }, 
                            { 8.004615534931510e-01 }, 
                            { 1.214745714812395e+00 }, 
                            { 2.870136978485235e-01 } };
                    break;

                // QUAD 9
                case 2:
                    tUHat = {
                            { 1.791257748677946e-01 }, 
                            { 4.877803160800752e-01 }, 
                            { 9.497310701550628e-01 }, 
                            { 3.487662544989165e-01 }, 
                            { 2.663465173645272e-01 }, 
                            { 6.897051952149583e-01 }, 
                            { 5.185891161035706e-01 }, 
                            { 2.532779069808584e-01 }, 
                            { 3.766051451809073e-01 } }; 
                    break;

                // QUAD 16
                case 3:
                    tUHat = {
                            { 4.636710644592379e-02 }, 
                            { 2.081211307185713e-01 }, 
                            { 1.109819586967827e+00 }, 
                            { 2.472556378443535e-01 }, 
                            { 8.458271849984333e-02 }, 
                            { 1.362546314963505e-01 }, 
                            { 4.523496496623921e-01 }, 
                            { 7.458874455200480e-01 }, 
                            { 7.265867638121548e-01 }, 
                            { 4.510428969226029e-01 }, 
                            { 1.661755462489430e-01 }, 
                            { 1.007785431698161e-01 }, 
                            { 1.838398770408118e-01 }, 
                            { 2.961483757533354e-01 }, 
                            { 4.883243651242546e-01 }, 
                            { 3.031368685109676e-01 } }; 
                    break;

                // unknown element type
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;

            } // end: switch( aInterpOrder )

            break;

        } // end: case QUAD elements

        // HEX elements
        case 3:
        {
            switch ( aInterpOrder )
            {
                // HEX 8
                case 1:
                    tUHat = {
                            { 9.267610502465551e-02 }, 
                            { 2.439025363679148e-01 }, 
                            { 7.653773237010404e-01 }, 
                            { 2.908218597932484e-01 }, 
                            { 1.424151185716165e-01 }, 
                            { 3.748043643775664e-01 }, 
                            { 1.176153251994260e+00 }, 
                            { 4.469051610947053e-01 } }; 
                    break;

                // HEX 27
                case 2:
                    tUHat = {
                            { 3.312113147023253e-04 }, 
                            { 6.334696322713911e-03 }, 
                            { 1.511088937009458e-01 }, 
                            { 7.900769475948958e-03 }, 
                            { 3.302093524490366e-03 }, 
                            { 6.315532947793853e-02 }, 
                            { 1.506517673864036e+00 }, 
                            { 7.876868502656571e-02 }, 
                            { 3.119892441439171e-03 }, 
                            { 6.301821507082660e-02 }, 
                            { 7.442243025942562e-02 }, 
                            { 3.294924460539922e-03 }, 
                            { 1.713434950294086e-03 }, 
                            { 3.277089156387216e-02 }, 
                            { 7.817222669466271e-01 }, 
                            { 4.087256066863042e-02 }, 
                            { 3.110454314413069e-02 }, 
                            { 6.282757583246801e-01 }, 
                            { 7.419729161648713e-01 }, 
                            { 3.284956836276042e-02 }, 
                            { 1.605618548342546e-01 }, 
                            { 3.103701311892221e-02 }, 
                            { 3.094312165380745e-01 }, 
                            { 1.704542833731989e-02 }, 
                            { 3.260082231929431e-01 }, 
                            { 1.613994604962295e-02 }, 
                            { 3.850049422584800e-01 } };
                    break;

                // HEX 64
                case 3:
                    tUHat = {
                            { +6.656965818039648e-02 }, 
                            { +2.847775346875384e-01 }, 
                            { +1.068177497728809e+00 }, 
                            { +2.496974031951603e-01 }, 
                            { +2.824841540827873e-01 }, 
                            { +1.208435542360656e+00 }, 
                            { +4.532743972313932e+00 }, 
                            { +1.059575212585738e+00 }, 
                            { +9.554564439898793e-02 }, 
                            { +1.634622216689139e-01 }, 
                            { +9.943399000309366e-02 }, 
                            { +1.571179566518752e-01 }, 
                            { +1.191622895413277e-01 }, 
                            { +1.869138356753473e-01 }, 
                            { +4.253674618621530e-01 }, 
                            { +6.721330043368123e-01 }, 
                            { +5.097629155814876e-01 }, 
                            { +7.995964344352301e-01 }, 
                            { +6.131335714633280e-01 }, 
                            { +3.583839837119830e-01 }, 
                            { +1.912079463003381e+00 }, 
                            { +2.999221548374744e+00 }, 
                            { +4.469681093543825e-01 }, 
                            { +7.010986786638649e-01 }, 
                            { +4.054419276301917e-01 }, 
                            { +6.936416480839619e-01 }, 
                            { +4.219418774388689e-01 }, 
                            { +6.667201588409398e-01 }, 
                            { +1.805020048515997e+00 }, 
                            { +2.852154094688128e+00 }, 
                            { +2.601793715167082e+00 }, 
                            { +1.520779875440484e+00 }, 
                            { +1.427149381518955e-01 }, 
                            { +2.255071878884365e-01 }, 
                            { +3.858041480218798e-01 }, 
                            { +2.441607987720849e-01 }, 
                            { +1.710304371915451e-01 }, 
                            { +2.926037645378185e-01 }, 
                            { +4.589681196402470e-01 }, 
                            { +2.682724136616744e-01 }, 
                            { +1.779907277710398e-01 }, 
                            { +2.791900841313841e-01 }, 
                            { +4.411547352653093e-01 }, 
                            { +2.812472822371492e-01 }, 
                            { +7.614243791746569e-01 }, 
                            { +1.203144343268560e+00 }, 
                            { +1.887210500377549e+00 }, 
                            { +1.194343880400987e+00 }, 
                            { +1.097533052854655e+00 }, 
                            { +6.415213357162445e-01 }, 
                            { +1.006268123815368e+00 }, 
                            { +1.721552292081377e+00 }, 
                            { +6.056018564729154e-01 }, 
                            { +1.036080980232837e+00 }, 
                            { +1.637135616653726e+00 }, 
                            { +9.569255566496553e-01 }, 
                            { +2.554653162833403e-01 }, 
                            { +4.370573710785387e-01 }, 
                            { +6.906045013517397e-01 }, 
                            { +4.036666786539957e-01 }, 
                            { +4.007140374050494e-01 }, 
                            { +6.855530382381665e-01 }, 
                            { +1.083258275576735e+00 }, 
                            { +6.331775558523358e-01 } }; 
                    break;

                // unknown element type
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;

            } // end: switch ( aInterpOrder )

            break;

        } // end: case HEX element

        default:
            MORIS_ERROR( false, "can only be 2D or 3D" );
            break;

    } // end: switch ( aSpaceDim )

    // repeat the coefficients to convert the coefficient vector to a time-wise linear element
    uint tNumCoeffs = tUHat.numel();
    aUHat.set_size( 2 * tNumCoeffs, 1, 0.0 );
    aUHat( { 0, tNumCoeffs - 1 }, { 0, 0 } ) = tUHat.matrix_data();
    aUHat( { tNumCoeffs, 2 * tNumCoeffs - 1 }, { 0, 0 } ) = tUHat.matrix_data();
}

//------------------------------------------------------------------------------

inline void
fill_u_hat_follower(
        moris::Matrix< moris::DDRMat >& aUHat,
        moris::uint                     aSpaceDim,
        moris::uint                     aInterpOrder )
{
    // initialize the vector for time-wise constant coefficient vector
    moris::Matrix< moris::DDRMat > tUHat;

    // switch between 2D and 3D
    switch ( aSpaceDim )
    {
        // QUAD elements
        case 2:
        {
            switch ( aInterpOrder )
            {
                // QUAD 4
                case 1:
                    tUHat = {
                            { 9.487233505307406e-02 }, 
                            { 1.724177562493723e-01 }, 
                            { 1.392912522939273e+00 }, 
                            { 7.664457910285476e-01 } }; 
                    break;

                // QUAD 9
                case 2:
                    tUHat = {
                            { 1.303471759442581e-02 }, 
                            { 8.500719406360592e-02 }, 
                            { 4.170085041779180e-01 }, 
                            { 6.394268327885244e-02 }, 
                            { 3.696396188739107e-02 }, 
                            { 1.732590347872419e-01 }, 
                            { 1.813291995453574e-01 }, 
                            { 2.656695840877514e-02 }, 
                            { 7.533880431025375e-02 } }; 
                    break;

                // QUAD 16
                case 3:
                    tUHat = {
                            { 3.848460715289934e-01 }, 
                            { 1.616122872670105e+00 }, 
                            { 4.589535336286909e+00 }, 
                            { 1.092902448311588e+00 }, 
                            { 5.773839439483581e-01 }, 
                            { 9.367582091291389e-01 }, 
                            { 1.904980209928541e+00 }, 
                            { 2.812197070858981e+00 }, 
                            { 2.660246306180847e+00 }, 
                            { 1.639679790545612e+00 }, 
                            { 6.696662818077329e-01 }, 
                            { 4.536314425896545e-01 }, 
                            { 6.805825258415442e-01 }, 
                            { 1.104189464833704e+00 }, 
                            { 1.630042329308591e+00 }, 
                            { 1.004699248671561e+00 } }; 
                    break;

                // unknown element type
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;

            } // end: switch( aInterpOrder )

            break;

        } // end: case QUAD elements

        // HEX elements
        case 3:
        {
            switch ( aInterpOrder )
            {
                // HEX 8
                case 1:
                    tUHat = {
                            { 1.378343839147868e-02 }, 
                            { 2.754323816862920e-02 }, 
                            { 1.509386926946396e-01 }, 
                            { 7.553411690047626e-02 }, 
                            { 2.363584315373315e-02 }, 
                            { 4.723115080647133e-02 }, 
                            { 2.588297030852331e-01 }, 
                            { 1.295259201012668e-01 } }; 
                    break;

                // HEX 27
                case 2:
                    tUHat = {
                            { 2.475306049631932e-02 }, 
                            { 3.184042307637808e-01 }, 
                            { 8.073815157335688e-01 }, 
                            { 6.276663929566220e-02 }, 
                            { 1.239825910908335e-01 }, 
                            { 1.594816186477111e+00 }, 
                            { 4.043994977282789e+00 }, 
                            { 3.143841778715433e-01 }, 
                            { 1.326016682045605e-01 }, 
                            { 5.337186140025686e-01 }, 
                            { 3.362396774910447e-01 }, 
                            { 4.149181406517944e-02 }, 
                            { 6.053596598082164e-02 }, 
                            { 7.786878590036176e-01 }, 
                            { 1.974528360937797e+00 }, 
                            { 1.535017918975191e-01 }, 
                            { 6.641723519163621e-01 }, 
                            { 2.673278186642266e+00 }, 
                            { 1.684149984163968e+00 }, 
                            { 2.078232959364748e-01 }, 
                            { 5.435845208378987e-01 }, 
                            { 2.222708485964533e-01 }, 
                            { 1.113305392561229e+00 }, 
                            { 1.014721813937216e-01 }, 
                            { 1.305259681540996e+00 }, 
                            { 3.242900035179521e-01 }, 
                            { 8.223061419426082e-01 } };
                    break;

                // HEX 64
                case 3:
                    tUHat = {
                            { +1.723415102815002e-01 }, 
                            { +5.763935435022657e-01 }, 
                            { +1.404280002528744e+00 }, 
                            { +4.198793328311483e-01 }, 
                            { +1.208489071085999e+00 }, 
                            { +4.041773202690752e+00 }, 
                            { +9.847059092314190e+00 }, 
                            { +2.944267948403790e+00 }, 
                            { +2.419881115028242e-01 }, 
                            { +3.622686617452231e-01 }, 
                            { +1.984964374226820e-01 }, 
                            { +2.734974617010738e-01 }, 
                            { +3.790996352747366e-01 }, 
                            { +7.022796369527619e-01 }, 
                            { +6.638682970327714e-01 }, 
                            { +9.147080748640684e-01 }, 
                            { +1.267892928172149e+00 }, 
                            { +2.348763497613022e+00 }, 
                            { +8.826029419770330e-01 }, 
                            { +5.895608472092060e-01 }, 
                            { +3.088994844670331e+00 }, 
                            { +5.722343089282838e+00 }, 
                            { +9.236086052378908e-01 }, 
                            { +1.710979002928463e+00 }, 
                            { +1.696863324490054e+00 }, 
                            { +2.540291760243716e+00 }, 
                            { +1.391892033921481e+00 }, 
                            { +1.917812446320359e+00 }, 
                            { +4.655161604967653e+00 }, 
                            { +6.414094375184315e+00 }, 
                            { +6.188967520044351e+00 }, 
                            { +4.134104658992348e+00 }, 
                            { +2.787127602252971e-01 }, 
                            { +3.840231767132378e-01 }, 
                            { +5.749024670802075e-01 }, 
                            { +4.172473516615597e-01 }, 
                            { +5.323012700869338e-01 }, 
                            { +7.968824069996742e-01 }, 
                            { +1.476219535469093e+00 }, 
                            { +9.860846804435262e-01 }, 
                            { +4.366326308000966e-01 }, 
                            { +8.088591412592574e-01 }, 
                            { +1.114483085341452e+00 }, 
                            { +6.016123904803277e-01 }, 
                            { +1.460311151181780e+00 }, 
                            { +2.012083432467447e+00 }, 
                            { +3.727371622101034e+00 }, 
                            { +2.705217018599314e+00 }, 
                            { +1.941461768841301e+00 }, 
                            { +1.296857047290633e+00 }, 
                            { +2.402419342057303e+00 }, 
                            { +3.596545444290479e+00 }, 
                            { +1.954383039549350e+00 }, 
                            { +2.925812032161925e+00 }, 
                            { +4.031317511793834e+00 }, 
                            { +2.692838256689527e+00 }, 
                            { +6.130844831012670e-01 }, 
                            { +9.178190360284113e-01 }, 
                            { +1.264613007235831e+00 }, 
                            { +8.447358154819512e-01 }, 
                            { +1.135735063162753e+00 }, 
                            { +1.700253863191518e+00 }, 
                            { +2.342687465166507e+00 }, 
                            { +1.564867667012482e+00 } }; 
                    break;

                // unknown element type
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;

            } // end: switch ( aInterpOrder )

            break;

        } // end: case HEX element

        default:
            MORIS_ERROR( false, "can only be 2D or 3D" );
            break;

    } // end: switch ( aSpaceDim )

    // repeat the coefficients to convert the coefficient vector to a time-wise linear element
    uint tNumCoeffs = tUHat.numel();
    aUHat.set_size( 2 * tNumCoeffs, 1, 0.0 );
    aUHat( { 0, tNumCoeffs - 1 }, { 0, 0 } ) = tUHat.matrix_data();
    aUHat( { tNumCoeffs, 2 * tNumCoeffs - 1 }, { 0, 0 } ) = tUHat.matrix_data();
}

//------------------------------------------------------------------------------

inline void
fill_normal(
        moris::Matrix< moris::DDRMat >& aNormal,
        moris::uint                     aSpaceDim,
        moris::uint                     aInterpOrder )
{
    // switch between 2D and 3D
    switch ( aSpaceDim )
    {
        // QUAD elements
        case 2:
        {
            switch ( aInterpOrder )
            {
                // QUAD 4
                case 1:
                    aNormal = {
                            { 9.502409374721200e-01 }, 
                            { 3.115159077031324e-01 } }; 
                    break;

                // QUAD 9
                case 2:
                    aNormal = {
                            { 2.906603354064973e-01 }, 
                            { 9.568263005485282e-01 } }; 
                    break;

                // QUAD 16
                case 3:
                    aNormal = {
                            { 1.310387997096721e-01 }, 
                            { 9.913772404945801e-01 } }; 
                    break;

                // unknown element type
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;

            } // end: switch( aInterpOrder )

            break;

        } // end: case QUAD elements

        // HEX elements
        case 3:
        {
            switch ( aInterpOrder )
            {
                // HEX 8
                case 1:
                    aNormal = {
                            { 2.738023870522189e-01 }, 
                            { 7.676889175321526e-01 }, 
                            { 5.793841374622010e-01 } }; 
                    break;

                // HEX 27
                case 2:
                    aNormal = {
                            { 7.441817958833307e-01 }, 
                            { 6.646968079454152e-01 }, 
                            { 6.611813807902152e-02 } };
                    break;

                // HEX 64
                case 3:
                    aNormal = {
                            { +8.377871164062863e-01 }, 
                            { +3.373332099529522e-01 }, 
                            { +4.293239488387261e-01 } }; 
                    break;

                // unknown element type
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;

            } // end: switch ( aInterpOrder )

            break;

        } // end: case HEX element

        default:
            MORIS_ERROR( false, "can only be 2D or 3D" );
            break;

    } // end: switch ( aSpaceDim )

}

//------------------------------------------------------------------------------

inline real
get_ref_value(
        moris::uint                     aSpaceDim,
        moris::uint                     aInterpOrder )
{
    // switch between 2D and 3D
    switch ( aSpaceDim )
    {
        // QUAD elements
        case 2:
        {
            switch ( aInterpOrder )
            {
                // QUAD 4
                case 1:
                    return 3.756416404069108e-02;
                    break;

                // QUAD 9
                case 2:
                    return 1.133495039809610e-01;
                    break;

                // QUAD 16
                case 3:
                    return 1.450355450930599e+02;
                    break;

                // unknown element type
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;

            } // end: switch( aInterpOrder )

            break;

        } // end: case QUAD elements

        // HEX elements
        case 3:
        {
            switch ( aInterpOrder )
            {
                // HEX 8
                case 1:
                    return 1.947987887813786e-01;
                    break;

                // HEX 27
                case 2:
                    return 6.142088389808256e+00;
                    break;

                // HEX 64
                case 3:
                    return 2.743572558692109e+02;
                    break;

                // unknown element type
                default:
                    MORIS_ERROR( false, "can only be 1, 2, 3" );
                    break;

            } // end: switch ( aInterpOrder )

            break;

        } // end: case HEX element

        default:
            MORIS_ERROR( false, "can only be 2D or 3D" );
            return 0.0;
            break;

    } // end: switch ( aSpaceDim )

    MORIS_ASSERT( false, "get_ref_value() - Unknown spatial dimension and polynomial order." );
    return 0.0;
    
}
