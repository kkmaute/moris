#include "cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Dirichlet.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::IWG_Isotropic_Struc_Linear_Pressure_Dirichlet( sint aBeta )
        {
            // sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            mPropertyMap[ "Select" ] = static_cast< uint >( IWG_Property_Type::SELECT );
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "DirichletNitsche" ] = static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here pressure), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for given dof type
            // FIXME protect dof type
            Field_Interpolator * tFIDispl =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // get the selection matrix property
            const std::shared_ptr< Property > & tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFIDispl->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get the imposed displacement property
            const std::shared_ptr< Property > & tPropDirichlet =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // get CM for elasticity
            const std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // compute displacement jump
            Matrix< DDRMat > tJump = tFIDispl->val() - tPropDirichlet->val();

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            mBeta * trans( tCMElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tM * tJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here pressure), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for given dof type
            // FIXME protect dof type
            Field_Interpolator * tFIDispl =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // get the selection matrix property
            const std::shared_ptr< Property > & tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFIDispl->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get the imposed displacement property
            const std::shared_ptr< Property > & tPropDirichlet =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get CM for elasticity
            const std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // compute jump
            Matrix< DDRMat > tJump = tFIDispl->val() - tPropDirichlet->val();

            // compute the jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for this dof type
                uint tDofDepIndex = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iDOF )( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if direct dependency on the dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::UX )
                {
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    mBeta * trans( tCMElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tM * tFIDispl->N() );
                }

                // if imposed displacement property depends on dof type
                if ( tPropDirichlet->check_dof_dependency( tDofType ) )
                {
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    mBeta * trans( tCMElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tM * tPropDirichlet->dPropdDOF( tDofType ) );
                }

                // if elasticity CM depends on dof type
                if( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    mBeta * tCMElasticity->dTestTractiondDOF( tDofType, mNormal, tM * tJump, mResidualDofType( 0 ) ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
