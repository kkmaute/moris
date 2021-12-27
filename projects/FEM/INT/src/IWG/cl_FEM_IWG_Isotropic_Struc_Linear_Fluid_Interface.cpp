
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Fluid_Interface.hpp"

#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Linear_Fluid_Interface::IWG_Isotropic_Struc_Linear_Fluid_Interface()
        {
            // set size for the property pointer cell
             mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive  map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
            mConstitutiveMap[ "Fluid" ]       = static_cast< uint >( IWG_Constitutive_Type::FLUID );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here displacement), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for the residual dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMSlaveFluid =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex } ) -= aWStar * (
                            tFIMaster->N_trans() * tCMSlaveFluid->traction( mNormal ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here displacement), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( MSI::Dof_Type::UX, mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type (here fluid), indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( MSI::Dof_Type::VX, mtk::Master_Slave::SLAVE );

            // get field interpolator for the residual dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type(MSI::Dof_Type::UX );

            // get slave fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMSlaveFluid =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // get number of master dof dependencies
            const  uint tNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();

            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get dependent dof type
                const Cell< MSI::Dof_Type > & tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                // get the index for the dof type
                const sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                const uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                const uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                // if dependency on the dof type
                if ( tCMSlaveFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    - tFIMaster->N_trans() * tCMSlaveFluid->dTractiondDOF( tDofType, mNormal ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_jacobian_and_residual - Not implemented.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
