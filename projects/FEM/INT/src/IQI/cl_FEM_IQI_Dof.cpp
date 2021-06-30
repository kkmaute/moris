/*
 * cl_FEM_IQI_Dof.cpp
 *
 *  Created on: Jan 23, 2020
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Dof.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Dof::IQI_Dof(){}

        //------------------------------------------------------------------------------

        void IQI_Dof::initialize()
        {
            if ( ! mIsInitialized )
            {
                // size of parameter list
                uint tParamSize = mParameters.size();

                // extract spatial derivative order
                if ( tParamSize > 0 )
                {
                    MORIS_ERROR( mParameters( 0 ).numel() == 2,
                            "IQI_Dof::initialize - Spatial gradient definition requires exactly two coefficients.\n");

                    mSpatialDerivativeDirection = mParameters( 0 )( 0 );
                    mSpatialDerivativeOrder     = mParameters( 0 )( 1 );
                }

                // extract time derivative order
                if ( tParamSize > 1 )
                {
                    MORIS_ERROR( mSpatialDerivativeOrder == 0,
                            "IQI_Dof::initialize - Time gradient can only be computed if spatial gradient order is zero.\n");

                    MORIS_ERROR( mParameters( 1 ).numel() == 1,
                            "IQI_Dof::initialize - Time gradient definition requires exactly one coefficient.\n");

                    mTimeDerivativeOrder = mParameters( 1 )( 0 );
                }

                // set initialize flag to true
                mIsInitialized = true;
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Dof::compute_QI( real aWStar )
        {
            // initialize if needed
            this->initialize();

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // check if dof index was set (for the case of vector field)
            if( mQuantityDofType.size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Dof::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            mSet->get_QI()( tQIIndex ) += aWStar * this->evaluate_QI( mIQITypeIndex );
        }

        //------------------------------------------------------------------------------

        void IQI_Dof::compute_QI( Matrix< DDRMat > & aQI )
        {
            // initialize if needed
            this->initialize();

            // check if dof index was set (for the case of vector field)
            sint tQIIndex = 0;
            if( mQuantityDofType.size() > 1 && mIQITypeIndex != -1 )
            {
                tQIIndex = mIQITypeIndex;
            }

            // evaluate QI
            aQI = { this->evaluate_QI( tQIIndex ) };
        }

        //------------------------------------------------------------------------------

        real IQI_Dof::evaluate_QI( sint aIQITypeIndex )
        {
            // get field interpolator for a given dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // evaluate spatial derivative of dof
            if ( mSpatialDerivativeOrder > 0 )
            {
                const Matrix<DDRMat> & tSpatialGradient = tFI->gradx( mSpatialDerivativeOrder );

                return tSpatialGradient( mSpatialDerivativeDirection, aIQITypeIndex );
            }

            // evaluate time derivative of dof
            if ( mTimeDerivativeOrder > 0 )
            {
                const Matrix<DDRMat> & tTemporalGradient = tFI->gradt( mTimeDerivativeOrder );

                return tTemporalGradient( aIQITypeIndex );
            }

            // evaluate DOF value
            return tFI->val()( aIQITypeIndex );
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */
