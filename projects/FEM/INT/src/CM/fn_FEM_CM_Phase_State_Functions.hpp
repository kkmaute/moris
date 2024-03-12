/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_CM_Phase_State_Functions.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_FN_FEM_CM_PHASE_STATE_FUNCTIONS_HPP_
#define PROJECTS_FEM_INT_SRC_FN_FEM_CM_PHASE_STATE_FUNCTIONS_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        inline real
        eval_phase_state_function(
                const real&         aTmelt,
                const real&         aPCconst,
                const uint&         aPCfunc,
                Field_Interpolator* aFieldInterpolator )
        {
            // get temperature
            real tTemp = aFieldInterpolator->val()( 0 );

            real tF = 0.0;

            switch ( aPCfunc )
            {
                // linear phase change model
                case 1:
                {
                    real tTupper = aTmelt + aPCconst / 2.0;
                    real tTlower = aTmelt - aPCconst / 2.0;

                    if ( tTemp < tTlower )
                    {
                        return 0.0;
                    }
                    if ( tTemp > tTupper )
                    {
                        return 1.0;
                    }
                    else
                    {
                        return ( tTemp - tTlower ) / ( tTupper - tTlower );
                    }
                    break;
                }

                // cubic phase change model
                case 2:
                {
                    real tTupper = aTmelt + aPCconst / 2.0;
                    real tTlower = aTmelt - aPCconst / 2.0;

                    if ( tTemp < tTlower )
                    {
                        return 0.0;
                    }
                    if ( tTemp > tTupper )
                    {
                        return 1.0;
                    }
                    else
                    {
                        return ( std::pow( tTemp - tTlower, 2.0 ) * ( 2.0 * tTemp + tTlower - 3.0 * tTupper ) ) / std::pow( tTlower - tTupper, 3.0 );
                    }
                    break;
                }

                // logistic function
                case 3:
                {
                    MORIS_ERROR( false, "eval_phase_state_function - logistic function option not implemented.\n" );
                    break;
                }

                // default - error
                default:
                {
                    MORIS_ERROR( false, "eval_phase_state_function - option not implemented." );
                }
            }

            return tF;
        }

        //------------------------------------------------------------------------------
        // derivative of phase state function
        inline real
        eval_dFdTemp(
                const real&         aTmelt,
                const real&         aPCconst,
                const uint&         aPCfunc,
                Field_Interpolator* aFieldInterpolator )
        {
            // get temperature
            real tTemp = aFieldInterpolator->val()( 0 );

            real tdfdT;

            switch ( aPCfunc )
            {
                // linear phase change model
                case 1:
                {
                    real tTupper = aTmelt + aPCconst / 2.0;
                    real tTlower = aTmelt - aPCconst / 2.0;

                    if ( ( tTemp > tTupper ) || ( tTemp < tTlower ) )
                    {
                        tdfdT = 0.0;
                    }
                    else
                    {
                        tdfdT = 1.0 / aPCconst;
                    }
                    break;
                }

                // cubic phase change model
                case 2:
                {
                    real tTupper = aTmelt + aPCconst / 2.0;
                    real tTlower = aTmelt - aPCconst / 2.0;

                    // cubic function: f(T)=( (tTemp - tTlower)^2 * (2*tTemp + tTlower - 3*tTupper))/(tTlower - tTupper)^3
                    if ( ( tTemp > tTupper ) || ( tTemp < tTlower ) )
                    {
                        tdfdT = 0;
                    }
                    else
                    {
                        tdfdT = 6.0 * ( tTemp - tTlower ) * ( tTupper - tTemp ) / std::pow( aPCconst, 3.0 );
                    }
                    break;
                }

                // logistic function
                case 3:
                {
                    // logistic function parameter k
                    //  real tLogisticParam = ( 2 * std::log(1/aPCconst - 1) ) / ( aTupper - 3 * aTlower );
                    real tLogisticParam = aPCconst;

                    real tExp = tLogisticParam * ( tTemp - aTmelt );
                    tdfdT     = tLogisticParam / ( std::exp( -tExp ) + 2.0 + std::exp( tExp ) );
                    break;
                }
                default:
                    MORIS_ERROR( false, "wrong option for phase change function." );
            }

            return tdfdT;
        }

        //------------------------------------------------------------------------------
        // 2nd derivative of phase state function wrt temperature
        inline real
        eval_d2FdTemp2(
                const real&         aTmelt,
                const real&         aPCconst,
                const uint&         aPCfunc,
                Field_Interpolator* aFieldInterpolator )
        {
            // get temperature
            real tTemp = aFieldInterpolator->val()( 0 );

            real td2fdT2;

            switch ( aPCfunc )
            {
                // linear phase change model
                case 1:
                {
                    td2fdT2 = 0;
                    break;
                }

                // cubic phase change model
                case 2:
                {
                    real tTupper = aTmelt + aPCconst / 2.0;
                    real tTlower = aTmelt - aPCconst / 2.0;

                    // cubic function: f(T)=( (tTemp - tTlower)^2 * (2*tTemp + tTlower - 3*tTupper))/(tTlower - tTupper)^3
                    if ( ( tTemp > tTupper ) || ( tTemp < tTlower ) )
                    {
                        td2fdT2 = 0;
                    }
                    else
                    {
                        td2fdT2 = ( 6.0 * ( tTlower + tTupper ) - ( 12.0 * tTemp ) ) / std::pow( aPCconst, 3.0 );
                    }
                    break;
                }

                // logistic function
                case 3:
                {
                    // logistic function parameter k
                    //  real tLogisticParam = ( 2.0 * std::log(1.0/aPCconst - 1.0) ) / ( aTupper - 3.0 * aTlower );
                    real tLogisticParam = aPCconst;

                    real tExp = tLogisticParam * ( tTemp - aTmelt );
                    td2fdT2   = std::pow( tLogisticParam, 2.0 ) *    //
                              ( ( std::exp( -tExp ) - 1.0 ) /        //
                                      ( std::exp( tExp ) + 3.0 + 3.0 * std::exp( -tExp ) + std::exp( -2.0 * tExp ) ) );
                    break;
                }
                default:
                    MORIS_ERROR( false, "wrong option for phase change function." );
            }

            return td2fdT2;
        }

        //------------------------------------------------------------------------------
        // derivative of the phase state function wrt temperature wrt the DoFs
        inline Matrix< DDRMat >
        eval_dFdTempdDOF(
                const real&         aTmelt,
                const real&         aPCconst,
                const uint&         aPCfunc,
                Field_Interpolator* aFieldInterpolator )
        {
            // get second derivative of phase state function
            real td2fdT2 = eval_d2FdTemp2( aTmelt, aPCconst, aPCfunc, aFieldInterpolator );

            // multiply with shape functions
            return aFieldInterpolator->N() * td2fdT2;
        }
    }    // namespace fem
}    // namespace moris

#endif /* PROJECTS_FEM_INT_SRC_FN_FEM_CM_PHASE_STATE_FUNCTIONS_HPP_ */

