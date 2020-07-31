
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.hpp"
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
        IWG_Isotropic_Struc_Linear_Neumann::IWG_Isotropic_Struc_Linear_Neumann()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Traction" ] = IWG_Property_Type::TRACTION;
            mPropertyMap[ "Pressure" ] = IWG_Property_Type::PRESSURE;
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Isotropic_Struc_Linear_Neumann::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Isotropic_Struc_Linear_Neumann::set_property - No slave allowed." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }


        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here displacement), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for the residual dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get traction load property
            std::shared_ptr< Property > tPropTraction =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );

            // get traction load property
            std::shared_ptr< Property > tPropPressure =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) );

            // compute the residual
            if (tPropTraction != nullptr)
            {
                mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) -= aWStar * (
                        trans( tFI->N() ) * tPropTraction->val() );
            }

            if (tPropPressure != nullptr)
            {
                mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) -= aWStar * (
                        trans( tFI->N() ) * mNormal *  tPropPressure->val());
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here displacement), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for the residual dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get traction load property
            std::shared_ptr< Property > tPropTraction =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );

            // get traction load property
            std::shared_ptr< Property > tPropPressure =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) );

            // compute the Jacobian for indirect IWG dof dependencies through properties
            for( uint iDOF = 0; iDOF < mRequestedMasterGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if traction load depends on the dof type
                if ( tPropTraction != nullptr )
                {
                    if ( tPropTraction->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        trans( tFI->N() ) * tPropTraction->dPropdDOF( tDofType ) );
                    }
                }

                // if pressure depends on the dof type
                if ( tPropPressure != nullptr )
                {
                    if ( tPropPressure->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        trans( tFI->N() ) * mNormal * tPropPressure->dPropdDOF( tDofType ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian_and_residual - Not implemented.");
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Neumann::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
