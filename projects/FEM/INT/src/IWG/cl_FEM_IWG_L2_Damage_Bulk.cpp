/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_L2_Damage_Bulk.cpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_L2_DAMAGE_BULK_CPP_
#define SRC_FEM_CL_FEM_IWG_L2_DAMAGE_BULK_CPP_

#include "cl_FEM_IWG_L2_Damage_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic_Damage.hpp"

// LINALG/src
#include "fn_trans.hpp"
#include "fn_diag_vec.hpp"
#include "fn_diag_mat.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_L2_Damage_Bulk::IWG_L2_Damage_Bulk( uint aSourceType )
    {
        // set the property pointer cell size
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Weight" ] = static_cast< uint >( IWG_Property_Type::WEIGHT );
        mPropertyMap[ "Lump" ]   = static_cast< uint >( IWG_Property_Type::LUMP );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElasticDamage" ] = static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE );

        // assign source type
        mSourceType = aSourceType;

        // switch on source type
        switch ( mSourceType )
        {
            case 0:    // equivalent strain
            {
                m_get_source    = &IWG_L2_Damage_Bulk::get_source_eqStrain;
                m_get_dsourcedu = &IWG_L2_Damage_Bulk::get_dsourcedu_eqStrain;
                break;
            }
            case 1:    // history
            {
                m_get_source    = &IWG_L2_Damage_Bulk::get_source_history;
                m_get_dsourcedu = &IWG_L2_Damage_Bulk::get_dsourcedu_history;
                break;
            }
            case 2:    // smooth damage
            {
                m_get_source    = &IWG_L2_Damage_Bulk::get_source_smoothDam;
                m_get_dsourcedu = &IWG_L2_Damage_Bulk::get_dsourcedu_smoothDam;
                break;
            }
            default:
                MORIS_ERROR( false, "IWG_L2_Damage_Bulk:: unknown source." );
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG_L2_Damage_Bulk::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
    {
        // FIXME for now only implemented for the specific case of equivalent strain
        // if not for equivalent strain then default parameters and no higher order derivatives included
        if ( mSourceType == 0 )
        {
            // set parameters
            mParameters = aParameters;

            // check characteristic length provided
            MORIS_ERROR( mParameters( 0 ).numel() == 1 || mParameters( 0 ).numel() == 2,
                    "IWG_Nonlocal_Bulk::set_parameters - IWG_Nonlocal_Bulk requires a characteristic length.\n" );

            // set a characteristic length
            mCharacteristicLength = aParameters( 0 )( 0 );

            // set a required order
            mOrder = static_cast< uint >( aParameters( 1 )( 0 ) );
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG_L2_Damage_Bulk::compute_residual( real aWStar )
    {

        // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 1 );

            // get residual dof type field interpolator
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get weight coefficient
            const std::shared_ptr< Property >& tPropWeight =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT ) );

            // get lump coefficient
            const std::shared_ptr< Property >& tPropLump =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::LUMP ) );

            // if weight property specified
            real tWeight = 1.0;
            if ( tPropWeight != nullptr )
            {
                // get weight value from property
                tWeight = tPropWeight->val()( 0 );
            }

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // add source
            tRes -= aWStar * tWeight * tFI->N_trans() * this->get_source();

            // lumping parameter
            real tConsistent = 1.0;
            real tLump       = 0.0;

            // if lumping
            if ( tPropLump != nullptr )
            {
                // compute weights for lumped and consistent terms
                tLump       = tPropLump->val()( 0 );
                tConsistent = 1.0 - tLump;

                // add lumped contribution
                tRes += aWStar * tWeight * tLump * sum( tFI->N_trans() * tFI->N(), 1 ) % vectorize( tFI->get_coeff() );
            }

            // add consistent contribution
            tRes += aWStar * tWeight * tConsistent * tFI->N_trans() * tFI->val();

            // diffusion only implemented for the specific case of equivalent strain
            if ( mSourceType == 0 )
            {
                // loop over the interpolation order
                for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
                {
                    // compute diffusion coefficient
                    real tDiffusion = std::pow( mCharacteristicLength, 2.0 * iOrder ) / mOrderCoeff( iOrder - 1.0 );

                    // compute the residual
                    tRes += aWStar * tWeight * tDiffusion * trans( tFI->dnNdxn( iOrder ) ) * tFI->gradx( iOrder );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_L2_Damage_Bulk::compute_residual - Residual contains NAN or INF, exiting!" );
        }
        //------------------------------------------------------------------------------

        void IWG_L2_Damage_Bulk::compute_jacobian( real aWStar )
        {

#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for a given dof type
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get weight coefficient
            const std::shared_ptr< Property >& tPropWeight =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT ) );

            // get lump coefficient
            const std::shared_ptr< Property >& tPropLump =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::LUMP ) );

            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            // if weight property specified
            real tWeight = 1.0;
            if( tPropWeight != nullptr )
            {
                // get weight value from property
                tWeight = tPropWeight->val()( 0 );
            }

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over leader dof type dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const moris::Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // lumping parameter
                real tConsistent = 1.0;
                real tLump       = 0.0;

                // if derivative dof type is residual type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    if ( tPropLump != nullptr )
                    {
                        // compute weights for lumped and consistent terms
                        tLump       = tPropLump->val()( 0 );
                        tConsistent = 1.0 - tLump;

                        // add lumped contribution to the Jacobian
                        tJac += aWStar * tWeight * tLump * diag_mat( sum( tFI->N_trans() * tFI->N(), 1 ) );
                    }

                    // add consistent contribution to the Jacobian
                    tJac += aWStar * tWeight * tConsistent * ( tFI->N_trans() * tFI->N() );

                    // diffusion only implemented for the specific case of equivalent strain
                    if ( mSourceType == 0 )
                    {
                        // loop over the interpolation order
                        for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
                        {
                            // compute diffusion coefficient
                            real tDiffusion = std::pow( mCharacteristicLength, 2.0 * iOrder ) / mOrderCoeff( iOrder - 1.0 );

                            // compute the jacobian
                            tJac += aWStar * tWeight * tDiffusion * trans( tFI->dnNdxn( iOrder ) ) * tFI->dnNdxn( iOrder );
                        }
                    }
                }

                // if weight property specified
                if( tPropWeight != nullptr )
                {
                    // if property has dependency on the dof type
                    if ( tPropWeight->check_dof_dependency( tDofType ) )
                    {
                        MORIS_ERROR( false, "IWG_L2_Damage_Bulk::compute_jacobian - dof dependency of property not supported" );
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( tCMElasticityDamagePtr->check_dof_dependency( tDofType ) )
                {
                    // add contribution from dsourcedu
                    tJac -= aWStar * tWeight * ( tFI->N_trans() * this->get_dsourcedu( tDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_L2_Damage_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_L2_Damage_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_L2_Damage_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_L2_Damage_Bulk::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_L2_Damage_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >& IWG_L2_Damage_Bulk::get_source()
        {
            return ( this->*m_get_source )();
        }

        const Matrix< DDRMat >& IWG_L2_Damage_Bulk::get_source_eqStrain()
        {
            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            return tCMElasticityDamagePtr->equivalent_strain();
        }

        const Matrix< DDRMat >& IWG_L2_Damage_Bulk::get_source_history()
        {
            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            return tCMElasticityDamagePtr->history();
        }

        const Matrix< DDRMat >& IWG_L2_Damage_Bulk::get_source_smoothDam()
        {
            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            return tCMElasticityDamagePtr->smooth_damage();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        IWG_L2_Damage_Bulk::get_dsourcedu( const moris::Vector< MSI::Dof_Type >& aDofType )
        {
            return ( this->*m_get_dsourcedu )( aDofType );
        }

        const Matrix< DDRMat >&
        IWG_L2_Damage_Bulk::get_dsourcedu_eqStrain( const moris::Vector< MSI::Dof_Type >& aDofType )
        {
            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            return tCMElasticityDamagePtr->dEqStraindu( aDofType );
        }

        const Matrix< DDRMat >&
        IWG_L2_Damage_Bulk::get_dsourcedu_smoothDam( const moris::Vector< MSI::Dof_Type >& aDofType )
        {
            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            return tCMElasticityDamagePtr->dSmoothDamagedu( aDofType );
        }

        const Matrix< DDRMat >&
        IWG_L2_Damage_Bulk::get_dsourcedu_history( const moris::Vector< MSI::Dof_Type >& aDofType )
        {
            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            return tCMElasticityDamagePtr->dHistorydu( aDofType );
        }

        //------------------------------------------------------------------------------

}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_L2_DAMAGE_BULK_CPP_ */
