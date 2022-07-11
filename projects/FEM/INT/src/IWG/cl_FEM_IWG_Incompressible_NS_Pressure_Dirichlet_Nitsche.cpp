
// FEM/INT/src
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche::IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            mPropertyMap[ "Select" ]    = static_cast< uint >( IWG_Property_Type::SELECT );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the veocity dof type
            // FIXME protect dof type
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get number of fields which should equal spatial dimension
                const uint tSpaceDim = tVelocityFI->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get the imposed velocity property
            const std::shared_ptr< Property >& tPropVelocity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tVelocityFI->val() - tPropVelocity->val();

            // compute master residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) -=    //
                    aWStar * ( mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tM * tVelocityJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the veocity dof type
            // FIXME protect dof type
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get number of fields which should equal spatial dimension
                const uint tSpaceDim = tVelocityFI->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get the imposed velocity property
            const std::shared_ptr< Property >& tPropVelocity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tVelocityFI->val() - tPropVelocity->val();

            // get number of master dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master
            for ( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if dof type is velocity
                // FIXME protect dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // compute jacobian direct dependencies
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -=    //
                            aWStar * ( mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tM * tVelocityFI->N() );
                }

                // if imposed velocity property depends on dof type
                if ( tPropVelocity->check_dof_dependency( tDofType ) )
                {
                    // add contribution of CM to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) +=    //
                            aWStar * ( mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tM * tPropVelocity->dPropdDOF( tDofType ) );
                }

                // if fluid constitutive model depends on dof type
                if ( tCMFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution of CM to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -=    //
                            aWStar * ( mBeta * tCMFluid->dTestTractiondDOF( tDofType, mNormal, tM * tVelocityJump, mResidualDofType( 0 ) ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
