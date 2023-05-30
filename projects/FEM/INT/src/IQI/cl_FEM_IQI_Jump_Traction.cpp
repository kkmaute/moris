/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Jump_Traction.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Jump_Traction.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Jump_Traction::IQI_Jump_Traction()
        {
            mFEMIQIType = fem::IQI_Type::JUMP_TRACTION;

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );
            mFollowerCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Jump_Traction::compute_QI( Matrix< DDRMat >& aQI )
        {
            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderElasticity =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model >& tCMFollowerElasticity =
                    mFollowerCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) );

            // evaluate average traction difference
            const Matrix< DDRMat > tTractionJump =
                    tCMLeaderElasticity->traction( mNormal ) - tCMFollowerElasticity->traction( - mNormal );

            // based on the IQI index select the norm or the individual component
            if ( mIQITypeIndex == -1 )
            {
                // compute norm if no index
                aQI = dot( tTractionJump, tTractionJump );
            }
            else
            {
                // pick the component otherwise (0,1,2)
                aQI = { tTractionJump( mIQITypeIndex ) * tTractionJump( mIQITypeIndex ) };
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI_Jump_Traction::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderElasticity =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model >& tCMFollowerElasticity =
                    mFollowerCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) );

            // evaluate average traction difference
            const Matrix< DDRMat > tTractionJump =
                    tCMLeaderElasticity->traction( mNormal ) - tCMFollowerElasticity->traction( mNormal );

            // initialize the matrix
            Matrix< DDRMat > tMat;

            // based on the IQI index select the norm or the individual component
            if ( mIQITypeIndex == -1 )
            {
                // compute norm if no index
                tMat = dot( tTractionJump, tTractionJump );
            }
            else
            {
                // pick the component otherwise (0,1,2)
                tMat = { tTractionJump( mIQITypeIndex ) * tTractionJump( mIQITypeIndex ) };
            }

            //add the contribution
            mSet->get_QI()( tQIIndex ) += aWStar * tMat;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Jump_Traction::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR( false, "Not Implemented for psedudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Jump_Traction::compute_dQIdu(
                moris::Cell< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&             adQIdu )
        {
            MORIS_ERROR( false, "Not Implemented for psedudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------

    } /* end namespace fem */
} /* end namespace moris */
