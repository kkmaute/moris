
#include "cl_FEM_IWG_Incompressible_NS_Pressure_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        IWG_Incompressible_NS_Pressure_Bulk::IWG_Incompressible_NS_Pressure_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Density" ] = IWG_Property_Type::DENSITY;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] = IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "IncompressibleFlow" ] = IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW;
        }

//------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the pressure FIs
            Field_Interpolator * tPressureFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX ); //FIXME this need to be protected

            // get the density and viscosity properties
            std::shared_ptr< Property > tDensityProp = mMasterProp( static_cast< uint >( IWG_Property_Type::DENSITY ) );

            // get the incompressible fluid constitutive model
            std::shared_ptr< Constitutive_Model > tIncFluidCM
            = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the incompressible flow stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tIncFlowSP
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW ) );

            // get the density value
            real tDensity = tDensityProp->val()( 0 );

            // compute the residual weak form
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
            += aWStar
            * ( trans( tPressureFI->N() ) * tVelocityFI->div() );

            // compute the residual strong form
            Matrix< DDRMat > tRM
            = tDensity * tVelocityFI->gradt( 1 )
            + tDensity * trans( tVelocityFI->val() ) * tVelocityFI->gradx( 1 )
            + tIncFluidCM->divflux();

            // add stabilization PSPG to the residual weak form
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
            += aWStar
            *  ( trans( tPressureFI->gradx( 1 ) ) * tIncFlowSP->val()( 0 ) * tRM );
        }

//------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian( real tWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian - Not implemented." );
        }

//------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

//------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_dRdp - Not implemented." );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */





