/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Jump_Dof.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Jump_Dof.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Jump_Dof::IQI_Jump_Dof()
        {
            mFEMIQIType = fem::IQI_Type::JUMP_DOF;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Jump_Dof::initialize()
        {
            if ( !mIsInitialized )
            {
                // size of parameter list
                uint tParamSize = mParameters.size();

                // extract spatial derivative order
                if ( tParamSize > 0 )
                {
                    MORIS_ERROR( mParameters( 0 ).numel() == 2,
                            "IQI_Dof::initialize - Spatial gradient definition requires exactly two coefficients.\n" );

                    mSpatialDerivativeDirection = mParameters( 0 )( 0 );
                    mSpatialDerivativeOrder     = mParameters( 0 )( 1 );
                }

                // extract time derivative order
                if ( tParamSize > 1 )
                {
                    MORIS_ERROR( mSpatialDerivativeOrder == 0,
                            "IQI_Dof::initialize - Time gradient can only be computed if spatial gradient order is zero.\n" );

                    MORIS_ERROR( mParameters( 1 ).numel() == 1,
                            "IQI_Dof::initialize - Time gradient definition requires exactly one coefficient.\n" );

                    mTimeDerivativeOrder = mParameters( 1 )( 0 );
                }

                // set initialize flag to true
                mIsInitialized = true;
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI_Jump_Dof::compute_QI( Matrix< DDRMat >& aQI )
        {
            // initialize if needed
            this->initialize();

            // evaluate QI
            this->evaluate_QI( aQI );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Jump_Dof::compute_QI( real aWStar )
        {
            // initialize if needed
            this->initialize();

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            Matrix< DDRMat > tMat;

            this->evaluate_QI( tMat );

            mSet->get_QI()( tQIIndex ) += aWStar * tMat;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Jump_Dof::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR( false, "Not Implemented for psedudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Jump_Dof::compute_dQIdu(
                moris::Cell< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&             adQIdu )
        {
            MORIS_ERROR( false, "Not Implemented for psedudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Jump_Dof::evaluate_QI( Matrix< DDRMat >& aMat )
        {
            // get field interpolator for master side
            Field_Interpolator* tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // get field interpolator for slave side
            Field_Interpolator* tFISlave =
                    mSlaveFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // check that field interpolater exists
            MORIS_ASSERT( tFIMaster != nullptr and tFISlave != nullptr,
                    "IQI_Jump_Dof::compute_QI - field interpolator does not exist." );

            // evaluate spatial derivative of dof
            if ( mSpatialDerivativeOrder > 0 )
            {
                const Matrix< DDRMat >& tSpatialGradient = tFIMaster->gradx( mSpatialDerivativeOrder ) - tFISlave->gradx( mSpatialDerivativeOrder ) ;

                aMat = { tSpatialGradient( mSpatialDerivativeDirection, mIQITypeIndex ) };
            }

            // evaluate time derivative of dof
            else if ( mTimeDerivativeOrder > 0 )
            {
                const Matrix< DDRMat >& tTemporalGradient = tFIMaster->gradt( mTimeDerivativeOrder ) - tFISlave->gradt( mTimeDerivativeOrder ) ;

                aMat = { tTemporalGradient( mIQITypeIndex ) };
            }
            else if ( mQuantityDofType.size() > 1 && mIQITypeIndex != -1 )
            {
                // Jump
                moris::real tJump = tFIMaster->val()( mIQITypeIndex ) - tFISlave->val()( mIQITypeIndex );

                // evaluate DOF value
                aMat = { tJump * tJump };
            }
            // DO NOT DELETE THIS FUNCTIONALITY AGAIN
            else
            {
                // Find the jump
                Matrix< DDRMat > tJump = tFIMaster->val() - tFISlave->val();

                // dot the product
                aMat = dot( tJump, tJump );
            }
        }

        //------------------------------------------------------------------------------

    } /* end namespace fem */
} /* end namespace moris */

