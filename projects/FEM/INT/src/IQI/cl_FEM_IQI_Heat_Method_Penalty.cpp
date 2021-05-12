/*
 * cl_FEM_IQI_Heat_Method_Penalty.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: noel
 */
#include "cl_FEM_IQI_Heat_Method_Penalty.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

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
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "L2_Reference" ]    = static_cast< uint >( IQI_Property_Type::L2_REFERENCE_VALUE );
            mPropertyMap[ "H1S_Reference" ]   = static_cast< uint >( IQI_Property_Type::H1S_REFERENCE_VALUE );
        }

        //------------------------------------------------------------------------------

        void IQI_Heat_Method_Penalty::initialize()
        {
            if ( ! mIsInitialized )
            {

                // check for proper size of constant function parameters
                MORIS_ERROR( mParameters.size() == 6 ,
                        "IQI_Heat_Method_Penalty::initialize - Needs 6 constant parameters." );

                // get weights for L2 and H1 semi-norm contributions
                mPhiBound       = mParameters( 0 )( 0 );
                mGammaPerimReg  = mParameters( 1 )( 0 );
                mWeightPhi1     = mParameters( 2 )( 0 );
                mWeightPhi2     = mParameters( 3 )( 0 );
                mWeightDelPhi1  = mParameters( 4 )( 0 );
                mWeightDelPhi2  = mParameters( 5 )( 0 );

                // check mQuantityDofType is defined
                MORIS_ERROR( mQuantityDofType.size() > 0,
                        "IQI_Heat_Method_Penalty::initialize - dof_quantity parameter needs to be defined." );

                // set initialize flag to true
                mIsInitialized = true;
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Heat_Method_Penalty::compute_QI( Matrix< DDRMat > & aQI )
        {
            // initialize if needed
            this->initialize();

            // get field interpolator
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // initialize QI
            aQI.fill(0.0);

            // Compute phi tilde
            moris::real tPhiTilde = 2/(1+std::exp(-2*tFI->val()(0)/mPhiBound)) - 1 * mPhiBound;

            // Compute alpha
            moris::real tAlpha = std::exp((-mGammaPerimReg) * std::pow(tPhiTilde / mPhiBound,2));

            // Compute weights
            moris::real tWPhi    = mWeightPhi1 * tAlpha + mWeightPhi2 * (1-tAlpha);
            moris::real tWDelPhi = mWeightDelPhi1 * tAlpha + mWeightDelPhi2 * (1-tAlpha);


            const std::shared_ptr< Property > & tPropL2Value =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::L2_REFERENCE_VALUE ) );

            // compute difference between dof value and reference value
            auto tL2error = tFI->val() - tPropL2Value->val();

            // compute L2 error
            aQI += tWPhi * trans( tL2error ) * tL2error;


            const std::shared_ptr< Property > & tPropH1SValue =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::H1S_REFERENCE_VALUE ) );

            // compute difference between dof spatial gradient and reference value and flatten it
            Matrix<DDRMat> tH1Serror = vectorize( tFI->gradx( 1 ) - tPropH1SValue->val() );

            // compute H1 semi-norm error
            aQI += tWDelPhi * trans( tH1Serror ) * tH1Serror;

        }

        //------------------------------------------------------------------------------

        void IQI_Heat_Method_Penalty::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            Matrix<DDRMat> tQI(1,1);

            this->compute_QI( tQI );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI(0);
        }

        //------------------------------------------------------------------------------

        void IQI_Heat_Method_Penalty::compute_dQIdu( real aWStar )
        {
            // check whether to compute derivatives
            // if ( mSkipComputeDQIDU )
            if(true)
            {
                return;
            }

            MORIS_ERROR(0,"Intentionally not implemented - Sensitivitiy not used for heat method regularizations");
        }

        //------------------------------------------------------------------------------

        void IQI_Heat_Method_Penalty::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // check whether to compute derivatives
            if ( true )
            {
                return;
            }

            MORIS_ERROR(0,"Intentionally not implemented - Sensitivitiy not used for heat method regularizations");

        }

        //------------------------------------------------------------------------------

    }/* end_namespace_fem */
}/* end_namespace_moris */
