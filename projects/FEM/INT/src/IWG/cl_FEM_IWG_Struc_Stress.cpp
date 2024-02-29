/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Struc_Stress.cpp
 *
 */

#include "cl_FEM_IWG_Struc_Stress.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
// LINALG/src
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Struc_Stress::IWG_Struc_Stress( enum Stress_Type aStressType )
        {
            // set stress type
            mStressType = aStressType;

            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Struc_Stress::compute_residual( real aWStar )
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
            Field_Interpolator* tFISig = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get stress vector from Constitutive model
            const Matrix< DDRMat >& tStressTensor =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->flux();

            Matrix< DDRMat > tDummy;
            moris::real      tStressVal;
            Matrix< DDRMat > tDStressVal;

            this->eval_stress_criterion(
                    tStressVal,
                    tDStressVal,
                    tStressTensor,
                    tDummy );

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute the residual
            tRes += aWStar * (                                     //
                            tFISig->N_trans() * tFISig->val() -    //
                            tFISig->N_trans() * tStressVal );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Bulk::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Struc_Stress::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 1 );

            // get field interpolator for a given dof type
            Field_Interpolator* tFISig = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get elasticity CM
            std::shared_ptr< Constitutive_Model >& tCMElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over leader dof type dependencies
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJac += aWStar * ( tFISig->N_trans() * tFISig->N() );
                }

                // if constitutive model has dependency on the dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // get stress vector from Constitutive model
                    const Matrix< DDRMat >& tStressTensor =
                            mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->flux();

                    const Matrix< DDRMat >& tDStressTensor =
                            mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->dFluxdDOF( tDofType );

                    moris::real      tStressVal;
                    Matrix< DDRMat > tDStressVal;

                    this->eval_stress_criterion(
                            tStressVal,
                            tDStressVal,
                            tStressTensor,
                            tDStressTensor );

                    // compute the Jacobian
                    tJac += ( -1.0 ) * aWStar * ( tFISig->N_trans() * tDStressVal );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Struc_Stress::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Struc_Stress::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Struc_Stress::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Struc_Stress::eval_stress_criterion(
                moris::real&             aStressVal,
                Matrix< DDRMat >&        aDStressVal,
                Matrix< DDRMat > const & aStressTensor,
                Matrix< DDRMat > const & aDStressTensor )
        {
            // switch for different stress types
            switch ( mStressType )
            {
                case Stress_Type::VON_MISES_STRESS:
                    this->eval_Von_Mises_stress(
                            aStressVal,
                            aDStressVal,
                            aStressTensor,
                            aDStressTensor );
                    break;

                default:
                    MORIS_ERROR( false, "IWG_Struc_Stress::compute_residual - Unknown Stress Type." );
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG_Struc_Stress::eval_Von_Mises_stress(
                moris::real&             aStressVal,
                Matrix< DDRMat >&        aDStressVal,
                Matrix< DDRMat > const & aStressTensor,
                Matrix< DDRMat > const & aDStressTensor )
        {

            uint tNumStressComponents = aStressTensor.n_rows();
            uint tNumNodes            = aDStressTensor.n_cols();

            moris::real tEpsilon = std::numeric_limits< real >::epsilon() * 10;

            switch ( tNumStressComponents )
            {
                // 2D plane stress
                case 3:
                {
                    // compute contributions to von mises stress
                    real tNormalStressContribution =
                            std::pow( aStressTensor( 0 ) - aStressTensor( 1 ), 2.0 ) +    //
                            std::pow( aStressTensor( 1 ), 2.0 ) +                         //
                            std::pow( aStressTensor( 0 ), 2.0 );
                    real tShearStressContribution =
                            std::pow( aStressTensor( 2 ), 2.0 );

                    // compute Von-Mises stress value
                    aStressVal = std::sqrt( 0.5 * tNormalStressContribution + 3.0 * tShearStressContribution );

                    if ( aDStressTensor.numel() > 0 )
                    {
                        aDStressVal.set_size( 1, tNumNodes, 0.0 );

                        aDStressVal += aDStressTensor.get_row( 0 ) * ( aStressTensor( 0 ) - aStressTensor( 1 ) );
                        aDStressVal += aDStressTensor.get_row( 1 ) * ( -1.0 ) * ( aStressTensor( 0 ) - aStressTensor( 1 ) );
                        aDStressVal += aDStressTensor.get_row( 1 ) * ( aStressTensor( 1 ) );
                        aDStressVal += aDStressTensor.get_row( 0 ) * ( -1.0 ) * ( -aStressTensor( 0 ) );
                        aDStressVal += aDStressTensor.get_row( 2 ) * 6.0 * ( aStressTensor( 2 ) );

                        if ( std::abs( aStressVal ) < 2.0 * tEpsilon )
                        {
                            aDStressVal = aDStressVal * ( 0.5 / ( 2.0 * tEpsilon + aStressVal ) );
                        }
                        else
                        {
                            aDStressVal = aDStressVal * ( 0.5 / ( aStressVal ) );
                        }
                    }
                }

                break;

                // 2D plane strain and axisymmetric
                case 4:
                {
                    real tNormalStressContribution =
                            std::pow( aStressTensor( 0 ) - aStressTensor( 1 ), 2.0 ) +    //
                            std::pow( aStressTensor( 1 ) - aStressTensor( 2 ), 2.0 ) +    //
                            std::pow( aStressTensor( 2 ) - aStressTensor( 0 ), 2.0 );
                    real tShearStressContribution =
                            std::pow( aStressTensor( 3 ), 2.0 );

                    // compute Von-Mises stress value
                    aStressVal = std::sqrt( 0.5 * tNormalStressContribution + 3.0 * tShearStressContribution );

                    if ( aDStressTensor.numel() > 0 )
                    {
                        aDStressVal.set_size( 1, tNumNodes, 0.0 );

                        aDStressVal += aDStressTensor.get_row( 0 ) * ( aStressTensor( 0 ) - aStressTensor( 1 ) );
                        aDStressVal += aDStressTensor.get_row( 1 ) * ( -1.0 ) * ( aStressTensor( 0 ) - aStressTensor( 1 ) );
                        aDStressVal += aDStressTensor.get_row( 1 ) * ( aStressTensor( 1 ) - aStressTensor( 2 ) );
                        aDStressVal += aDStressTensor.get_row( 2 ) * ( -1.0 ) * ( aStressTensor( 1 ) - aStressTensor( 2 ) );
                        aDStressVal += aDStressTensor.get_row( 2 ) * ( aStressTensor( 2 ) - aStressTensor( 0 ) );
                        aDStressVal += aDStressTensor.get_row( 0 ) * ( -1.0 ) * ( aStressTensor( 2 ) - aStressTensor( 0 ) );
                        aDStressVal += aDStressTensor.get_row( 3 ) * 6.0 * ( aStressTensor( 3 ) );

                        if ( std::abs( aStressVal ) < 2 * tEpsilon )
                        {
                            aDStressVal = aDStressVal * ( 0.5 / ( 2.0 * tEpsilon + aStressVal ) );
                        }
                        else
                        {
                            aDStressVal = aDStressVal * ( 0.5 / ( aStressVal ) );
                        }
                    }
                }

                break;

                // 3D
                case 6:
                {
                    real tNormalStressContribution =
                            std::pow( aStressTensor( 0 ) - aStressTensor( 1 ), 2.0 ) +    //
                            std::pow( aStressTensor( 1 ) - aStressTensor( 2 ), 2.0 ) +    //
                            std::pow( aStressTensor( 2 ) - aStressTensor( 0 ), 2.0 );
                    real tShearStressContribution =
                            std::pow( aStressTensor( 3 ), 2.0 ) +    //
                            std::pow( aStressTensor( 4 ), 2.0 ) +    //
                            std::pow( aStressTensor( 5 ), 2.0 );

                    // compute Von-Mises stress value
                    aStressVal = std::sqrt( 0.5 * tNormalStressContribution + 3.0 * tShearStressContribution );

                    if ( aDStressTensor.numel() > 0 )
                    {
                        aDStressVal.set_size( 1, tNumNodes, 0.0 );

                        aDStressVal += aDStressTensor.get_row( 0 ) * ( aStressTensor( 0 ) - aStressTensor( 1 ) );
                        aDStressVal += aDStressTensor.get_row( 1 ) * ( -1.0 ) * ( aStressTensor( 0 ) - aStressTensor( 1 ) );
                        aDStressVal += aDStressTensor.get_row( 1 ) * ( aStressTensor( 1 ) - aStressTensor( 2 ) );
                        aDStressVal += aDStressTensor.get_row( 2 ) * ( -1.0 ) * ( aStressTensor( 1 ) - aStressTensor( 2 ) );
                        aDStressVal += aDStressTensor.get_row( 2 ) * ( aStressTensor( 2 ) - aStressTensor( 0 ) );
                        aDStressVal += aDStressTensor.get_row( 0 ) * ( -1.0 ) * ( aStressTensor( 2 ) - aStressTensor( 0 ) );
                        aDStressVal += aDStressTensor.get_row( 3 ) * 6.0 * ( aStressTensor( 3 ) );
                        aDStressVal += aDStressTensor.get_row( 4 ) * 6.0 * ( aStressTensor( 4 ) );
                        aDStressVal += aDStressTensor.get_row( 5 ) * 6.0 * ( aStressTensor( 5 ) );

                        if ( std::abs( aStressVal ) < 2 * tEpsilon )
                        {
                            aDStressVal = aDStressVal * ( 0.5 / ( 2.0 * tEpsilon + aStressVal ) );
                        }
                        else
                        {
                            aDStressVal = aDStressVal * ( 0.5 / ( aStressVal ) );
                        }
                    }
                }

                break;

                // Unknown size - error
                default:
                    MORIS_ERROR( false,
                            "IWG_Struc_Stress::get_stress_vector - CM stress vector of unknown size; 3, 4 or 6 components expected." );
            }
        }

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
