
//FEM/INT/src
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche( sint aBeta )
                {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ]  = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            mPropertyMap[ "SlipLength" ] = static_cast< uint >( IWG_Property_Type::SLIPLENGTH );
            mPropertyMap[ "Traction" ]   = static_cast< uint >( IWG_Property_Type::TRACTION );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "DirichletNitsche" ] = static_cast< uint >( IWG_Stabilization_Type::VELOCITY_SPLIPLENGTH_NITSCHE );
                }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the master field interpolator for the residual dof type
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the imposed velocity property
            const std::shared_ptr< Property > & tPropVelocity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the slip length property
            const std::shared_ptr< Property > & tPropSlipLength =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SLIPLENGTH ) );

            // get the traction property
            const std::shared_ptr< Property > & tPropTraction =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::VELOCITY_SPLIPLENGTH_NITSCHE ) );

            // get the dynamic viscosity property
            const std::shared_ptr< Property > & tPropViscosity = tCMFluid->get_property( "Viscosity" );

            // check that slip length is defined
            MORIS_ASSERT( tPropSlipLength,
                    "IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_residual - Slip length not defined.\n");

            // check that prescribed velocity parameter is defined
            MORIS_ASSERT( tPropVelocity,
                    "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_residual - Prescribed velocity not defined.\n");

            // get the dynamic viscosity value
            const real tViscosity = tPropViscosity->val()( 0 );

            // get the slip length
            const real tSplipLength = tPropSlipLength->val()( 0 );

            // get spatial dimension
            uint tSpaceDim = tFIVelocity->get_number_of_fields();

            // build projector
            const Matrix<DDRMat> tNormalProjector  = mNormal * trans( mNormal );
            const Matrix<DDRMat> tTangentProjector = moris::eye(tSpaceDim,tSpaceDim) - tNormalProjector;

            // compute velocity jump in normal direction
            const Matrix<DDRMat> tNormalVelocityJump = tNormalProjector * ( tFIVelocity->val() - tPropVelocity->val() );

            // slip condition violation
            Matrix<DDRMat> tSlipVelocityJump = tTangentProjector * ( tSplipLength * tCMFluid->traction( mNormal )
                    + tViscosity * ( tFIVelocity->val() - tPropVelocity->val() ) );

            // add contribution of prescribed traction to slip condition violation
            if ( tPropTraction )
            {
                tSlipVelocityJump -= tSplipLength * tTangentProjector * tPropTraction->val();
            }

            // penalty parameters
            const real tNormalPenalty   = tSPNitsche->val()( 0 );
            const real tTangentPenalty1 = tSPNitsche->val()( 1 );
            const real tTangentPenalty2 = tSPNitsche->val()( 2 );

            // compute master residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex } ) += aWStar * (
                            + tFIVelocity->N_trans() * (
                                    - tCMFluid->traction( mNormal )
                                    + tNormalPenalty   * tNormalVelocityJump
                                    + tTangentPenalty1 * tSlipVelocityJump )
                            - mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * (
                                    + tNormalVelocityJump
                                    + tTangentPenalty2 * tSlipVelocityJump ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the master field interpolator for residual dof type
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the imposed velocity property
            const std::shared_ptr< Property > & tPropVelocity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the slip length property
            const std::shared_ptr< Property > & tPropSlipLength =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SLIPLENGTH ) );

            // get the traction property
            const std::shared_ptr< Property > & tPropTraction =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::VELOCITY_SPLIPLENGTH_NITSCHE ) );

            // get the dynamic viscosity property
            const std::shared_ptr< Property > & tPropViscosity = tCMFluid->get_property( "Viscosity" );

            // get the dynamic viscosity value
            const real tViscosity = tPropViscosity->val()( 0 );

            // get the slip length
            const real tSplipLength = tPropSlipLength->val()( 0 );

            // get spatial dimension
            uint tSpaceDim = tFIVelocity->get_number_of_fields();

            // build projector
            const Matrix<DDRMat> tNormalProjector  = mNormal * trans( mNormal );
            const Matrix<DDRMat> tTangentProjector = moris::eye(tSpaceDim,tSpaceDim) - tNormalProjector;

            // compute velocity jump in normal direction
            const Matrix<DDRMat> tNormalVelocityJump = tNormalProjector * ( tFIVelocity->val() - tPropVelocity->val() );

            // slip condition violation
            Matrix<DDRMat> tSlipVelocityJump = tTangentProjector * ( tSplipLength * tCMFluid->traction( mNormal )
                    + tViscosity * ( tFIVelocity->val() - tPropVelocity->val() ) );

            // add contribution of prescribed traction to slip condition violation
            if ( tPropTraction )
            {
                tSlipVelocityJump -= tSplipLength * tTangentProjector * tPropTraction->val();
            }

            // penalty parameters
            const real tNormalPenalty   = tSPNitsche->val()( 0 );
            const real tTangentPenalty1 = tSPNitsche->val()( 1 );
            const real tTangentPenalty2 = tSPNitsche->val()( 2 );

            // get number of master dependencies
            const uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute the Jacobian for indirect dof dependencies through master
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute Jacobian direct dependencies
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    + tFIVelocity->N_trans() * (
                                            + tNormalPenalty   * tNormalProjector
                                            + tTangentPenalty1 * tTangentProjector * tViscosity ) * tFIVelocity->N()
                                    - mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * (
                                            + tNormalProjector
                                            + tTangentPenalty2 * tTangentProjector * tViscosity  ) * tFIVelocity->N() );
                }

                // if imposed velocity depends on dof type
                if ( tPropVelocity->check_dof_dependency( tDofType ) )
                {
                    // add contribution from property to Jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    + tFIVelocity->N_trans() * (
                                            + tNormalPenalty   * tNormalProjector
                                            + tTangentPenalty1 * tTangentProjector * tViscosity ) * tPropVelocity->dPropdDOF( tDofType )
                                    - mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * (
                                            + tNormalProjector
                                            + tTangentPenalty2 * tTangentProjector * tViscosity  ) * tPropVelocity->dPropdDOF( tDofType ) );
                }

                // if fluid constitutive model depends on dof type
                if ( tCMFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution of CM to Jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - tFIVelocity->N_trans() * tCMFluid->dTractiondDOF( tDofType, mNormal )
                                    - mBeta * tCMFluid->dTestTractiondDOF(
                                            tDofType,
                                            mNormal,
                                            tNormalVelocityJump + tTangentPenalty2 * tSlipVelocityJump,
                                            mResidualDofType( 0 ) ) );

                    // add contribution due to dependency of SlipVelocityJump on velocity and pressure
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * tSplipLength * ( (
                                    + tTangentPenalty1 * tFIVelocity->N_trans()
                                    - tTangentPenalty2 * mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) ) *
                                    tTangentProjector * tCMFluid->dTractiondDOF( tDofType, mNormal ) );
                }

                // if prescribed traction depends on the dof type
                if ( tPropTraction )
                {
                    if ( tPropTraction->check_dof_dependency( tDofType ) )
                    {
                        MORIS_ERROR(false,"IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian - %s.\n",
                                "Dof dependency of prescribed traction not implemented yet");
                    }
                }

                // if viscosity depends on the dof type
                if ( tPropViscosity->check_dof_dependency( tDofType ) )
                {
                    MORIS_ERROR(false,"IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian - %s.\n",
                            "Dof dependency of viscosity not implemented yet");
                }

                // if stabilization parameter depends on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // get derivative of penalty parameter
                    const Matrix<DDRMat> & tDSPNitsche = tSPNitsche->dSPdMasterDOF( tDofType );

                    // add contribution of SP to Jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    + tFIVelocity->N_trans() * (
                                            + tNormalVelocityJump * tDSPNitsche.get_row( 0 ) ) );

                    // add the following lines if Nitsche penalty for tangential direction depends on dofs
                    //                        + tSlipVelocityJump   * tDSPNitsche.get_row( 1 ) )
                    //                - mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) *
                    //                        tSlipVelocityJump * tDSPNitsche.get_row( 2 ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
