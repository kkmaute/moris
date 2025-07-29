/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Turbulence_Heat_Coefficient.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Turbulence_Heat_Coefficient.hpp"
#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Turbulence_Heat_Coefficient::IQI_Turbulence_Heat_Coefficient()
    {
        // set fem IQI type
        mFEMIQIType = fem::IQI_Type::TURBULENCE_HEAT_COEFFICIENT;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "DiffusionTurbulence" ] =
                static_cast< uint >( IQI_Constitutive_Type::DIFF_TURBULENCE );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Turbulence_Heat_Coefficient::compute_QI( Matrix< DDRMat >& aQI )
    {
        // get the CM
        const std::shared_ptr< Constitutive_Model >& tCM =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::DIFF_TURBULENCE ) );

        // cast constitutive model base class pointer to constitutive model
        CM_Diffusion_Linear_Isotropic_Turbulence* tCMPtr =
                dynamic_cast< CM_Diffusion_Linear_Isotropic_Turbulence* >( tCM.get() );

        // switch on type index
        switch ( mIQITypeIndex )
        {
            // 0 - conductivity
            case 0:
            {
                // get conductivity property from CM
                const std::shared_ptr< Property >& tPropCond = tCMPtr->get_property( "Conductivity" );

                // compute conductivity property value
                aQI = tPropCond->val();
                break;
            }
            // 1 - effective conductivity
            case 1:
            {
                aQI = tCMPtr->effective_conductivity();
                break;
            }
            // if none of the above
            default:
            {
                MORIS_ERROR( false,                                         //
                        "IQI_Turbulence_Heat_Coefficient::compute_QI - "    //
                        "Wrong mIQITypeIndex type" );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Turbulence_Heat_Coefficient::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // init QI value matrix
        Matrix< DDRMat > tQI( 1, 1 );
        this->compute_QI( tQI );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * tQI( 0 );
    }
    //------------------------------------------------------------------------------
}    // namespace moris::fem
