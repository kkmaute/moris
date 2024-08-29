/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.cpp
 *
 */

#include "cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Thermal_Energy_Diffusive_Flux::IQI_Thermal_Energy_Diffusive_Flux()
    {
        // set fem IQI type
        mFEMIQIType = fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IQI_Constitutive_Type::DIFFUSION );
    }

    //------------------------------------------------------------------------------

    void IQI_Thermal_Energy_Diffusive_Flux::compute_QI( Matrix< DDRMat > &aQI )
    {
        // get the diffusion CM
        const std::shared_ptr< Constitutive_Model > &tCMDiffusion =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

        // evaluate the QI
        // note: because diffusive flux is defined as "conductivity tensor * gradx(Temp)"
        //       the numerical diffusive flux needs to be multiplied by (-1) to obtain
        //       the physical flux
        aQI = -1.0 * dot( tCMDiffusion->flux(), mNormal );
    }

    //------------------------------------------------------------------------------

    void IQI_Thermal_Energy_Diffusive_Flux::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the diffusion CM
        const std::shared_ptr< Constitutive_Model > &tCMDiffusion =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

        // evaluate the QI
        // note: because diffusive flux is defined as "conductivity tensor * gradx(Temp)"
        //       the numerical diffusive flux needs to be multiplied by (-1) to obtain
        //       the physical flux
        mSet->get_QI()( tQIIndex ) += -1.0 * aWStar * dot( tCMDiffusion->flux(), mNormal );
    }

    //------------------------------------------------------------------------------

    void IQI_Thermal_Energy_Diffusive_Flux::compute_dQIdu( real aWStar )
    {
        // get the column index to assemble in residual
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the diffusion CM
        const std::shared_ptr< Constitutive_Model > &tCMDiffusion =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

        // get the number of leader dof type dependencies
        uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // compute dQIdu for indirect dof dependencies
        for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
        {
            // get the treated dof type
            Vector< MSI::Dof_Type > &tDofType = mRequestedLeaderGlobalDofTypes( iDof );

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // if density depends on dof type
            if ( tCMDiffusion->check_dof_dependency( tDofType ) )
            {
                // compute dQIdu
                mSet->get_residual()( tQIIndex )(
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( -1.0 * trans( tCMDiffusion->dFluxdDOF( tDofType ) ) * mNormal );
            }
        }
    }

    //------------------------------------------------------------------------------

    void IQI_Thermal_Energy_Diffusive_Flux::compute_dQIdu(
            Vector< MSI::Dof_Type > &aDofType,
            Matrix< DDRMat >        &adQIdu )
    {
        // get the diffusion CM
        const std::shared_ptr< Constitutive_Model > &tCMDiffusion =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

        // if density depends on dof type
        if ( tCMDiffusion->check_dof_dependency( aDofType ) )
        {
            // compute dQIdu
            adQIdu += -1.0 * trans( tCMDiffusion->dFluxdDOF( aDofType ) ) * mNormal;
        }
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
