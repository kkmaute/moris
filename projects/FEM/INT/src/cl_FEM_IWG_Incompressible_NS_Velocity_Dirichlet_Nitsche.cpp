
//FEM/INT/src
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = IWG_Property_Type::DIRICHLET;
            mPropertyMap[ "Select" ]    = IWG_Property_Type::SELECT;
            mPropertyMap[ "Upwind" ]    = IWG_Property_Type::UPWIND;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] = IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE;
            mConstitutiveMap[ "TurbulenceFluid" ]     = IWG_Constitutive_Type::FLUID_TURBULENCE;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "DirichletNitsche" ] = IWG_Stabilization_Type::DIRICHLET_NITSCHE;
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::set_property - No slave allowed" );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // check that aConstitutiveString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::set_constitutive_model - Unknown aConstitutiveString: " ) +
                    aConstitutiveString;
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::set_property - No slave allowed" );

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::set_stabilization_parameter - Unknown aStabilizationString: " ) +
                    aStabilizationString;
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(), tErrMsg.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

//            // FIXME to remove
//            Matrix< DDRMat > tPrev = mSet->get_residual()( 0 )({ tMasterResStartIndex, tMasterResStopIndex },{ 0, 0 } );

            // get the master field interpolator for the residual dof type
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the selection matrix property
            std::shared_ptr< Property > tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFIVelocity->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get the imposed velocity property
            std::shared_ptr< Property > tPropVelocity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the upwind property
            std::shared_ptr< Property > tPropUpwind =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::UPWIND ) );

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMTurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_TURBULENCE ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFIVelocity->val() - tPropVelocity->val();

            // compute master residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) +=
                            aWStar * (
                                    - trans( tFIVelocity->N() ) * tM * tCMFluid->traction( mNormal )
                                    - mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType ) ) * tM * tVelocityJump
                                    + tSPNitsche->val()( 0 ) * trans( tFIVelocity->N() ) * tM * tVelocityJump );

            // upwind term
            if (tPropUpwind)
            {
                mSet->get_residual()( 0 )(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { 0, 0 } ) -=
                                aWStar * tPropUpwind->val()( 0 ) * trans( tFIVelocity->N() ) * mNormal * trans( tFIVelocity->val() ) * tM * tVelocityJump;
            }

            // if turbulence
            if( tCMTurbulence != nullptr )
            {
                mSet->get_residual()( 0 )(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { 0, 0 } ) +=
                                aWStar * (
                                        - trans( tFIVelocity->N() ) * tM * tCMTurbulence->traction( mNormal )
                                        - mBeta * trans( tCMTurbulence->testTraction( mNormal, mResidualDofType ) ) * tM * tVelocityJump );
            }

//            std::cout<<"4 "<<
//                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->valx()( 0 )<<" "<<
//                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->valx()( 1 )<< " "<<
//                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->valt()( 0 )<< " "<<
//                    norm( mSet->get_residual()( 0 )({ tMasterResStartIndex, tMasterResStopIndex },{ 0, 0 } ) - tPrev )<<std::endl;
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

//            // FIXME to remove
//            uint tEndIndex = mSet->get_jacobian().n_cols() - 1;
//            Matrix< DDRMat > tPrev = mSet->get_jacobian()({ tMasterResStartIndex, tMasterResStopIndex },{ 0, tEndIndex } );

            // get the master field interpolator for residual dof type
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the selection matrix property
            std::shared_ptr< Property > tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFIVelocity->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get the imposed velocity property
            std::shared_ptr< Property > tPropVelocity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the upwind property
            std::shared_ptr< Property > tPropUpwind =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::UPWIND ) );

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMTurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_TURBULENCE ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFIVelocity->val() - tPropVelocity->val();

            // get number of master dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    // compute jacobian direct dependencies
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } )
                            += aWStar * (
                                    - mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType ) ) * tM * tFIVelocity->N()
                                    + tSPNitsche->val()( 0 )  * trans( tFIVelocity->N() ) * tM * tFIVelocity->N() );
                }

                // if imposed velocity depends on dof type
                if ( tPropVelocity->check_dof_dependency( tDofType ) )
                {
                    // add contribution from property to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } )
                            -= aWStar * (
                                    - mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType ) ) * tM * tPropVelocity->dPropdDOF( tDofType )
                                    + tSPNitsche->val()( 0 ) * trans( tFIVelocity->N() ) * tM * tPropVelocity->dPropdDOF( tDofType ) );
                }

                // if fluid constitutive model depends on dof type
                if ( tCMFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution of CM to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - trans( tFIVelocity->N() ) * tM * tCMFluid->dTractiondDOF( tDofType, mNormal )
                                    - mBeta * tCMFluid->dTestTractiondDOF( tDofType, mNormal, tM * tVelocityJump, mResidualDofType ) );
                }

                // if stabilization parameter depends on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // add contribution of SP to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFIVelocity->N() ) * tM * tVelocityJump * tSPNitsche->dSPdMasterDOF( tDofType ) );
                }

                // upwind term
                if (tPropUpwind)
                {
                    // if dof type is residual dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 ) )
                    {
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        tPropUpwind->val()( 0 ) * trans( tFIVelocity->N() ) * mNormal * trans( tFIVelocity->val() ) * tM * tFIVelocity->N() +
                                        tPropUpwind->val()( 0 ) * trans( tFIVelocity->N() ) * mNormal * trans( tM * tFIVelocity->val() ) * tFIVelocity->N() );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropVelocity->check_dof_dependency( tDofType ) )
                    {
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        tPropUpwind->val()( 0 ) * trans( tFIVelocity->N() ) * mNormal * trans( tFIVelocity->val() ) * tM * tPropVelocity->dPropdDOF( tDofType ));
                    }

                    // if upwind parameter depends on the dof type
                    if ( tPropUpwind->check_dof_dependency( tDofType ) )
                    {
                        // add contribution of SP to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        trans( tFIVelocity->N() ) * mNormal * trans( tFIVelocity->val() ) * tM * tVelocityJump * tPropUpwind->dPropdDOF( tDofType ));
                    }
                }

                // if turbulence
                if( tCMTurbulence != nullptr )
                {
                    // if dof type is residual dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 ) )
                    {
                        // compute jacobian direct dependencies
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                       - mBeta * trans( tCMTurbulence->testTraction( mNormal, mResidualDofType ) ) * tM * tFIVelocity->N() );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropVelocity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution from property to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                       - mBeta * trans( tCMTurbulence->testTraction( mNormal, mResidualDofType ) ) * tM * tPropVelocity->dPropdDOF( tDofType ) );
                    }

                    // if turbulence depends on dof type
                    if( tCMTurbulence->check_dof_dependency( tDofType ) )
                    {
                        // add contribution from CM to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        - trans( tFIVelocity->N() ) * tM * tCMTurbulence->dTractiondDOF( tDofType, mNormal )
                                        - mBeta * tCMTurbulence->dTestTractiondDOF( tDofType, mNormal, tM * tVelocityJump, mResidualDofType ) );
                    }
                }
            }

//            std::cout<<"41 "<<
//                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->valx()( 0 )<<" "<<
//                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->valx()( 1 )<< " "<<
//                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->valt()( 0 )<< " "<<
//                    max(max( mSet->get_jacobian()({ tMasterResStartIndex, tMasterResStopIndex },{ 0, tEndIndex } ) - tPrev ) )<<std::endl;
//            std::cout<<"42 "<<
//                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->valx()( 0 )<<" "<<
//                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->valx()( 1 )<< " "<<
//                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->valt()( 0 )<< " "<<
//                    min( min( mSet->get_jacobian()({ tMasterResStartIndex, tMasterResStopIndex },{ 0, tEndIndex } ) - tPrev ) )<<std::endl;
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
