/*
 * Copyright (c) 2022 University of Colorado
 *Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Robin_Nitsche.cpp
 *
 */

// FEM/INT/src
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Diffusion_Robin_Nitsche.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Diffusion_Robin_Nitsche::IWG_Diffusion_Robin_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ]      = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            mPropertyMap[ "NeumannPenalty" ] = static_cast< uint >( IWG_Property_Type::NEUMANN_PENALTY );
            mPropertyMap[ "Traction" ]       = static_cast< uint >( IWG_Property_Type::TRACTION );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFFUSION );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "RobinNitsche" ] = static_cast< uint >( IWG_Stabilization_Type::ROBIN_NITSCHE );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Robin_Nitsche::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the master field interpolator for the residual dof type
            Field_Interpolator* tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the imposed velocity property
            const std::shared_ptr< Property >& tPropDirichlet =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the slip length property
            const std::shared_ptr< Property >& tPropNeumannPen =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::NEUMANN_PENALTY ) );

            // get the traction property
            const std::shared_ptr< Property >& tPropTraction =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the dynamic viscosity property
            const std::shared_ptr< Property >& tPropConductivity = tCMDiffusion->get_property( "Conductivity" );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::ROBIN_NITSCHE ) );

            // check that slip length is defined
            MORIS_ASSERT( tPropNeumannPen and tPropConductivity,
                    "IWG_Diffusion_Robin_Nitsche::compute_residual - Slip length not defined.\n" );

            // get the slip length
            const real tNeumannPen = tPropNeumannPen->val()( 0 );

            // compute the dirichlet jump
            Matrix< DDRMat > tJump = tPropConductivity->val() * tFITemp->val();
            if ( tPropDirichlet )
            {
                // subtract the prescribed dirichlet , by default is zero
                tJump -= tPropConductivity->val() * tPropDirichlet->val();
            }

            // compute the traction jump
            tJump += tNeumannPen * tCMDiffusion->traction( mNormal );
            if ( tPropTraction )
            {
                // subtract the prescribed traction , by default is zero
                tJump -= tNeumannPen * tPropTraction->val();
            }

            // penalty parameters
            const real tStabilityPenalty = tSPNitsche->val()( 0 );
            const real tAdjointPenalty   = tSPNitsche->val()( 1 );

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex } );

            // compute master residual
            tRes += aWStar * (                                                                            //
                            tFITemp->N_trans() * (                                                        //
                                    -tCMDiffusion->traction( mNormal )                                    //                                     //
                                    + tStabilityPenalty * tJump )                                         //
                            - mBeta * tCMDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * (    //                                                 //
                                      tAdjointPenalty * tJump ) );                                        //

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Robin_Nitsche::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Robin_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the master field interpolator for the residual dof type
            Field_Interpolator* tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the imposed velocity property
            const std::shared_ptr< Property >& tPropDirichlet =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the slip length property
            const std::shared_ptr< Property >& tPropNeumannPen =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::NEUMANN_PENALTY ) );

            // get the traction property
            const std::shared_ptr< Property >& tPropTraction =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the dynamic viscosity property
            const std::shared_ptr< Property >& tPropConductivity = tCMDiffusion->get_property( "Conductivity" );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::ROBIN_NITSCHE ) );

            // check that slip length is defined
            MORIS_ASSERT( tPropNeumannPen and tPropConductivity,
                    "IWG_Diffusion_Robin_Nitsche::compute_residual - Slip length not defined.\n" );

            // get the slip length
            const real tNeumannPen = tPropNeumannPen->val()( 0 );

            // compute the dirichlet jump
            Matrix< DDRMat > tJump = tPropConductivity->val() * tFITemp->val();
            if ( tPropDirichlet )
            {
                // subtract the prescribed dirichlet , by default is zero
                tJump -= tPropConductivity->val() * tPropDirichlet->val();
            }

            // compute the traction jump
            tJump += tNeumannPen * tCMDiffusion->traction( mNormal );
            if ( tPropTraction )
            {
                // subtract the prescribed traction , by default is zero
                tJump -= tNeumannPen * tPropTraction->val();
            }

            // penalty parameters
            const real tStabilityPenalty = tSPNitsche->val()( 0 );
            const real tAdjointPenalty   = tSPNitsche->val()( 1 );

            // get number of master dependencies
            const uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute the Jacobian for indirect dof dependencies through master
            for ( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute Jacobian direct dependencies
                    tJac += aWStar * (                                                                            //
                                    tFITemp->N_trans() * (                                                        //
                                            tStabilityPenalty * tPropConductivity->val() * tFITemp->N() )         //
                                    - mBeta * tCMDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * (    //
                                              tAdjointPenalty * tPropConductivity->val() * tFITemp->N() ) );      //
                }

                // if fluid constitutive model depends on dof type
                if ( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // compute Jacobian direct dependencies
                    tJac += aWStar * (                                                                                           //
                                    tFITemp->N_trans() * (                                                                       //
                                            -tCMDiffusion->dTractiondDOF( tDofType, mNormal ) )                                  //
                                    - mBeta * tCMDiffusion->dTestTractiondDOF( tDofType, mNormal, mResidualDofType( 0 ) ) * (    //
                                              tAdjointPenalty * tJump( 0 ) ) );

                    // compute the dependencies of the jacobian on the jump term which has traction in it
                    tJac += aWStar * (                                                                                              //
                                    tFITemp->N_trans() * (                                                                          //
                                            tStabilityPenalty * tNeumannPen * tCMDiffusion->dTractiondDOF( tDofType, mNormal ) )    //
                                    - mBeta * tCMDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * (                      //
                                              tAdjointPenalty * tNeumannPen * tCMDiffusion->dTractiondDOF( tDofType, mNormal ) ) );
                }

                // if prescribed traction depends on the dof type
                if ( tPropDirichlet )
                {
                    if ( tPropDirichlet->check_dof_dependency( tDofType ) )
                    {
                        MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_jacobian - %s.\n", "Dof dependency of prescribed traction not implemented yet" );
                    }
                }

                // if prescribed traction depends on the dof type
                if ( tPropTraction )
                {
                    if ( tPropTraction->check_dof_dependency( tDofType ) )
                    {
                        MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_jacobian - %s.\n", "Dof dependency of prescribed traction not implemented yet" );
                    }
                }

                if ( tPropNeumannPen->check_dof_dependency( tDofType ) )
                {
                    MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_jacobian - %s.\n", "Dof dependency of prescribed traction not implemented yet" );
                }

                // if stabilization parameter depends on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_jacobian - %s.\n", "Dof dependency of prescribed traction not implemented yet" );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Diffusion_Robin_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Robin_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Robin_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
