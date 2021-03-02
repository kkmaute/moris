/*
 * cl_FEM_IQI_H1_Error.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: noel
 */
#include "cl_FEM_IQI_H1_Error.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_H1_Error::IQI_H1_Error()
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

        void IQI_H1_Error::initialize()
        {
            if ( ! mIsInitialized )
            {
                // size of parameter list
                uint tParamSize = mParameters.size();

                // check for proper size of constant function parameters
                MORIS_ERROR( tParamSize == 2 || tParamSize == 3,
                        "IQI_H1_Error::IQI_H1_Error - either 2 or 3 constant parameters need to be set." );

                // get weights for L2 and H1 semi-norm contributions
                mL2Weight  = mParameters( 0 )( 0 );
                mH1SWeight = mParameters( 1 )( 0 );

                // extract parameter whether to skip computing dQIdu
                if ( tParamSize > 2 )
                {
                    mSkipComputeDQIDU = mParameters(2)(0) > 0;
                }

                // set initialize flag to true
                mIsInitialized = true;
            }
        }

        //------------------------------------------------------------------------------

        void IQI_H1_Error::compute_QI( Matrix< DDRMat > & aQI )
        {
            // initialize if needed
            this->initialize();

            // get field interpolator
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofTypes( 0 )( 0 ) );

            // initialize QI
            aQI.fill(0.0);

            // L2 contributions
            if ( mL2Weight > 0.0 )
            {
                MORIS_ASSERT (  mMasterProp( (uint) IQI_Property_Type::L2_REFERENCE_VALUE ) != nullptr,
                        "IQI_H1_Error::compute_QI - no weight for L2 contribution provided.\n");

                const std::shared_ptr< Property > & tPropL2Value =
                        mMasterProp( static_cast< uint >( IQI_Property_Type::L2_REFERENCE_VALUE ) );

                auto tL2error = tFI->val() - tPropL2Value->val();

                aQI += mL2Weight * trans( tL2error ) * tL2error;
            }

            // H1 semi-norm contribution
            if ( mH1SWeight )
            {
                MORIS_ASSERT (  mMasterProp( (uint) IQI_Property_Type::H1S_REFERENCE_VALUE ) != nullptr,
                        "IQI_H1_Error::compute_QI - no weight for H1 semi-norm contribution provided.\n");

                const std::shared_ptr< Property > & tPropH1SValue =
                        mMasterProp( static_cast< uint >( IQI_Property_Type::H1S_REFERENCE_VALUE ) );

                auto tL2error = tFI->gradx( 1 ) - tPropH1SValue->val();

                aQI += mH1SWeight * trans( tL2error ) * tL2error;
            }
        }

        //------------------------------------------------------------------------------

        void IQI_H1_Error::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            Matrix<DDRMat> tQI(1,1);

            this->compute_QI( tQI );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI(0);
        }

        //------------------------------------------------------------------------------

        void IQI_H1_Error::compute_dQIdu( real aWStar )
        {
            // check whether to compute derivatives
            if ( mSkipComputeDQIDU )
            {
                return;
            }

            MORIS_ERROR( false, "IQI_H1_Error::compute_dQIdu() - not implemented." );
        }

        //------------------------------------------------------------------------------

        void IQI_H1_Error::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // check whether to compute derivatives
            if ( mSkipComputeDQIDU )
            {
                return;
            }

            MORIS_ERROR( false, "IQI_H1_Error::compute_dQIdu() - not implemented." );
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */
