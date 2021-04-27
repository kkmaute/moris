/*
 * fn_FEM_CM_Phase_State_Functions.hpp
 *
 *  Created on: May 18, 2020
 *      Author: wunsch
 */

#ifndef PROJECTS_FEM_INT_SRC_FN_FEM_CM_PHASE_STATE_FUNCTIONS_HPP_
#define PROJECTS_FEM_INT_SRC_FN_FEM_CM_PHASE_STATE_FUNCTIONS_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        // derivative of phase state function
        inline
        moris::real eval_dFdTemp(
                moris::real aTmelt,
                moris::real aPCconst,
                moris::uint aPCfunc,
                Field_Interpolator * aFieldInterpolator)
        {
            // get temperature
            real tTemp = aFieldInterpolator->val()( 0 );

            moris::real tdfdT;

            switch ( aPCfunc )
            {
                // linear phase change model
                case 1:
                {
                    real tTupper = aTmelt + aPCconst/2.0;
                    real tTlower = aTmelt - aPCconst/2.0;

                    if ( (tTemp > tTupper) || (tTemp < tTlower) )
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
                    real tTupper = aTmelt + aPCconst/2.0;
                    real tTlower = aTmelt - aPCconst/2.0;

                    // cubic function: f(T)=( (tTemp - tTlower)^2 * (2*tTemp + tTlower - 3*tTupper))/(tTlower - tTupper)^3
                    if ( (tTemp > tTupper) || (tTemp < tTlower) )
                    {
                        tdfdT = 0;
                    }
                    else
                    {
                        tdfdT = 6.0 * (tTemp - tTlower) * (tTupper - tTemp) / std::pow(aPCconst, 3.0);
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
                    tdfdT = tLogisticParam  / ( std::exp(-tExp) + 2.0 + std::exp(tExp) );
                    break;
                }
                default:
                    MORIS_ERROR(false,"wrong option for phase change function.");
            }

            // debug
            //std::cout << "dfdT = " << tdfdT << "\n" << std::flush;
            //std::cout << "Temp = " << tTemp << "\n" << std::flush;

            return tdfdT;
        }

        //------------------------------------------------------------------------------
        // 2nd derivative of phase state function wrt temperature
        inline
        moris::real eval_d2FdTemp2(
                moris::real aTmelt,
                moris::real aPCconst,
                moris::uint aPCfunc,
                Field_Interpolator * aFieldInterpolator)
        {
            // get temperature
            real tTemp = aFieldInterpolator->val()( 0 );

            moris::real td2fdT2;

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
                    real tTupper = aTmelt + aPCconst/2.0;
                    real tTlower = aTmelt - aPCconst/2.0;

                    // cubic function: f(T)=( (tTemp - tTlower)^2 * (2*tTemp + tTlower - 3*tTupper))/(tTlower - tTupper)^3
                    if ( (tTemp > tTupper) || (tTemp < tTlower) )
                    {
                        td2fdT2 = 0;
                    }
                    else
                    {
                        td2fdT2 = ( 6.0 * (tTlower + tTupper) - ( 12.0 * tTemp )  )  / std::pow(aPCconst, 3.0);
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
                    td2fdT2 = std::pow(tLogisticParam, 2) *
                            ( ( std::exp(-tExp) - 1 ) / ( std::exp(tExp) + 3 + 3*std::exp(-tExp) + std::exp(- 2.0 * tExp) ) );
                    break;
                }
                default:
                    MORIS_ERROR(false,"wrong option for phase change function.");
            }

            return td2fdT2;
        }

        //------------------------------------------------------------------------------
        // derivative of the phase state function wrt temperature wrt the DoFs
        inline
        moris::Matrix<DDRMat> eval_dFdTempdDOF(
                moris::real aTmelt,
                moris::real aPCconst,
                moris::uint aPCfunc,
                Field_Interpolator * aFieldInterpolator)
        {
            // get second derivative of phase state function
            moris::real td2fdT2 = eval_d2FdTemp2( aTmelt, aPCconst, aPCfunc, aFieldInterpolator);

            // multiply with shape functions
            return aFieldInterpolator->N() * td2fdT2;
        }
    } // namespace fem
} // namespace moris

#endif /* PROJECTS_FEM_INT_SRC_FN_FEM_CM_PHASE_STATE_FUNCTIONS_HPP_ */
