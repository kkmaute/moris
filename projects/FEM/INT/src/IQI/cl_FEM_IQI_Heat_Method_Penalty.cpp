/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Heat_Method_Penalty.cpp
 *
 */

#include "cl_FEM_IQI_Heat_Method_Penalty.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Heat_Method_Penalty::IQI_Heat_Method_Penalty()
        {
            // set FEM IQI type
            mFEMIQIType = fem::IQI_Type::H1_ERROR;

            // set the property pointer cell size
            mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "L2_Reference" ]  = static_cast< uint >( IQI_Property_Type::L2_REFERENCE_VALUE );
            mPropertyMap[ "H1S_Reference" ] = static_cast< uint >( IQI_Property_Type::H1S_REFERENCE_VALUE );
            mPropertyMap[ "Select" ]        = static_cast< uint >( IQI_Property_Type::SELECT );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Heat_Method_Penalty::initialize()
        {
            if ( !mIsInitialized )
            {
                // check for proper size of constant function parameters
                MORIS_ERROR( mParameters.size() == 7 || mParameters.size() == 8,
                        "IQI_Heat_Method_Penalty::initialize - Needs 7 or 8 constant parameters." );

                // get weights for L2 and H1 semi-norm contributions
                mPhiBound      = mParameters( 0 )( 0 );
                mPhiGradient   = mParameters( 1 )( 0 );
                mPhiGamma      = mParameters( 2 )( 0 );
                mWeightPhi1    = mParameters( 3 )( 0 );
                mWeightPhi2    = mParameters( 4 )( 0 );
                mWeightDelPhi1 = mParameters( 5 )( 0 );
                mWeightDelPhi2 = mParameters( 6 )( 0 );

                if ( mParameters.size() == 8 )
                {
                    mLevelSetSign = mParameters( 7 )( 0 );
                }

                // check mQuantityDofType is defined
                MORIS_ERROR( mQuantityDofType.size() > 0,
                        "IQI_Heat_Method_Penalty::initialize - dof_quantity parameter needs to be defined." );

                // set initialize flag to true
                mIsInitialized = true;
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI_Heat_Method_Penalty::compute_QI( Matrix< DDRMat >& aQI )
        {
            // initialize QI
            aQI.fill( 0.0 );

            // get select property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) );

            // check if Heat Penalty is used
            if ( tPropSelect != nullptr && tPropSelect->val()( 0 ) < MORIS_REAL_EPS )
            {
                return;
            }

            // initialize if needed
            this->initialize();

            // get field interpolator
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // project level set field
            real tVal = std::exp( -2.0 * mPhiGradient * mLevelSetSign * tFI->val()( 0 ) / mPhiBound );

            // check for nan, infinity
            MORIS_ASSERT( std::isfinite( tVal ),
                    "IQI_Heat_Method_Penalty::compute_QI - project level set field value is NAN or INF, exiting!" );

            // Compute phi tilde
            moris::real      tPhiTilde   = ( 2.0 / ( 1.0 + tVal ) - 1.0 ) * mPhiBound;
            Matrix< DDRMat > tPhiTildeDx = ( 4.0 * mPhiGradient * tVal ) / std::pow( 1.0 + tVal, 2 ) * mLevelSetSign * tFI->gradx( 1 );

            // Compute alpha
            moris::real tAlpha = std::exp( -mPhiGamma * std::pow( tPhiTilde / mPhiBound, 2 ) );

            // Compute weights
            moris::real tWPhi    = mWeightPhi1 * tAlpha + mWeightPhi2 * ( 1.0 - tAlpha );
            moris::real tWDelPhi = mWeightDelPhi1 * tAlpha + mWeightDelPhi2 * ( 1.0 - tAlpha );

            const std::shared_ptr< Property >& tPropL2Value =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::L2_REFERENCE_VALUE ) );

            // compute difference between dof value and reference value
            auto tL2error = tPropL2Value->val() - tPhiTilde;

            // compute L2 error
            real tL2Contribution = tWPhi * dot( tL2error, tL2error );

            const std::shared_ptr< Property >& tPropH1SValue =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::H1S_REFERENCE_VALUE ) );

            // compute difference between dof spatial gradient and reference value and flatten it
            Matrix< DDRMat > tH1Serror = vectorize( tPropH1SValue->val() - tPhiTildeDx );

            // compute H1 semi-norm error
            real tH1Contribution = tWDelPhi * dot( tH1Serror, tH1Serror );

            // return desired value
            switch ( mIQITypeIndex )
            {
                // heat penalty
                case -1:
                case 0:
                {
                    aQI( 0 ) = tL2Contribution + tH1Contribution;
                    break;
                }
                // L2 contribution of heat penalty
                case 1:
                {
                    aQI( 0 ) = tL2Contribution;
                    break;
                }
                // H1 contribution of heat penalty
                case 2:
                {
                    aQI( 0 ) = tH1Contribution;
                    break;
                }
                //  projected level set field
                case 3:
                {
                    aQI( 0 ) = tPhiTilde;
                    break;
                }
                //  norm of spatial gradients of projected level set field
                case 4:
                {
                    aQI( 0 ) = norm( tPhiTildeDx );
                    break;
                }
                //  distance to interface measure
                case 5:
                {
                    aQI( 0 ) = tAlpha;
                    break;
                }
                //  L2 weight
                case 6:
                {
                    aQI( 0 ) = tWPhi;
                    break;
                }
                //  H1 weight
                case 7:
                {
                    aQI( 0 ) = tWDelPhi;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false,
                            "IQI_Heat_Method_Penalty::compute_QI - incorrect vector index." );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( aQI ),
                    "IQI_Heat_Method_Penalty::compute_QI - IQI is NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Heat_Method_Penalty::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            Matrix< DDRMat > tQI( 1, 1 );

            this->compute_QI( tQI );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI( 0 );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Heat_Method_Penalty::compute_dQIdu( real aWStar )
        {
            // check whether to compute derivatives
            // if ( mSkipComputeDQIDU )
            if ( true )
            {
                return;
            }

            MORIS_ERROR( 0, "Intentionally not implemented - Sensitivity not used for heat method regularization." );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Heat_Method_Penalty::compute_dQIdu(
                moris::Cell< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&             adQIdu )
        {
            // check whether to compute derivatives
            if ( true )
            {
                return;
            }

            MORIS_ERROR( 0, "Intentionally not implemented - Sensitivity not used for heat method regularization." );
        }

        //------------------------------------------------------------------------------

    }    // namespace fem
}    // namespace moris
