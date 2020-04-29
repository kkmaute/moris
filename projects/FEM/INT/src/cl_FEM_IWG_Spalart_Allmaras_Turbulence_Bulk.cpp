
#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        IWG_Spalart_Allmaras_Turbulence_Bulk::IWG_Spalart_Allmaras_Turbulence_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "WallDistance" ] = IWG_Property_Type::WALL_DISTANCE;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
//            mConstitutiveMap[ "IncompressibleFluid" ] = IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "SUPG" ] = IWG_Stabilization_Type::SUPG;
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity
            = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the velocity FI
            // FISME protect dof type
            Field_Interpolator * tFIVelocity
            = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropWallDistance
            = mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // compute the residual weak form
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
            += aWStar * (   trans( tFIViscosity->N() ) * tFIViscosity->gradt( 1 )
                          + trans( tFIViscosity->N() ) * trans( tFIVelocity->val() ) * tFIViscosity->gradx( 1 ) );
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity
            = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the velocity FI
            // FISME protect dof type
            Field_Interpolator * tFIVelocity
            = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropWallDistance
            = mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // get number of dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if residual dof type (here viscosity)
                if( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    // compute the jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * (   trans( tFIViscosity->N() ) * tFIViscosity->dnNdtn( 1 )
                                  + trans(tFIViscosity->N() ) * trans( tFIVelocity->val() ) * tFIViscosity->dnNdxn( 1 ) );
                }

                // if velocity dof type
                // FIXME protect dof type
                if( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // compute the jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * ( trans( tFIViscosity->N() ) * trans( tFIViscosity->gradx( 1 ) ) * tFIVelocity->N() );
                }
            }

        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_dRdp - Not implemented." );
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual_strong_form( Matrix< DDRMat > & aR )
        {
            MORIS_ASSERT( false, "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual_strong_form - not implemented.");
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_strong_form
        ( moris::Cell< MSI::Dof_Type >   aDofTypes,
          Matrix< DDRMat >             & aJ )
        {
            MORIS_ASSERT( false, "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_strong_form - not implemented.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
