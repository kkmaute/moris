
#include "cl_FEM_IWG_Diffusion_Radiation.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        IWG_Diffusion_Radiation::IWG_Diffusion_Radiation()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Emissivity" ] = IWG_Property_Type::EMISSIVITY;
            mPropertyMap[ "AmbientTemperature" ] = IWG_Property_Type::AMBIENT_TEMP;
            mPropertyMap[ "AbsoluteZero" ] = IWG_Property_Type::ABSOLUTE_ZERO;
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Radiation::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Diffusion_Radiation::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Diffusion_Radiation::set_property - No slave allowed." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Radiation::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get filed interpolator for residual dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get indices for SP, CM, properties
            uint tEmissivityIndex = static_cast< uint >( IWG_Property_Type::EMISSIVITY );
            uint tAmbTempIndex = static_cast< uint >( IWG_Property_Type::AMBIENT_TEMP );
            uint tZeroTempIndex = static_cast< uint >( IWG_Property_Type::ABSOLUTE_ZERO );

            // get property values
            real tT0 = mMasterProp( tZeroTempIndex )->val()( 0 );
            real tTinf = mMasterProp( tAmbTempIndex )->val()( 0 );
            real tT = tFI->val()( 0 );
            real tAlpha = mMasterProp( tEmissivityIndex )->val()( 0 );

            // compute the residual
            // N * a * (T - T_ref)
            mSet->get_residual()( 0 )( { tResStartIndex, tResStopIndex }, { 0, 0 } ) +=
                    aWStar * mStefanBoltzmannConst * tAlpha * ( std::pow( tT - tT0 , 4.0 ) - std::pow( tTinf - tT0 , 4.0 ) ) *
                    trans( tFI->N() );
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Radiation::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get indices for SP, CM, properties
            uint tEmissivityIndex = static_cast< uint >( IWG_Property_Type::EMISSIVITY );
            uint tAmbTempIndex = static_cast< uint >( IWG_Property_Type::AMBIENT_TEMP );
            uint tZeroTempIndex = static_cast< uint >( IWG_Property_Type::ABSOLUTE_ZERO );

            // get property values
            real tT0 = mMasterProp( tZeroTempIndex )->val()( 0 );
            real tTinf = mMasterProp( tAmbTempIndex )->val()( 0 );
            real tT = tFI->val()( 0 );
            real tAlpha = mMasterProp( tEmissivityIndex )->val()( 0 );

            // compute the jacobian for dof dependencies
            for( uint iDOF = 0; iDOF < mRequestedMasterGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDepDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the dof type indices for assembly
                uint tDepDofIndex   = mSet->get_dof_index_for_type( tDepDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 1 );

                // if dof type is residual dof type
                if( tDepDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    mSet->get_jacobian()(
                            { tResStartIndex, tResStopIndex },
                            { tDepStartIndex, tDepStopIndex } ) +=
                                    aWStar * tAlpha * mStefanBoltzmannConst *
                                    ( 4.0 * std::pow( tT , 3.0 ) - 12.0 * tT0 * std::pow( tT , 2.0 ) +
                                            2.0 * tT * std::pow( tT0 , 2.0 ) + 4.0 * std::pow( tT0 , 3.0 ) ) *
                                            trans( tFI->N() ) * tFI->N();
                }

                // if dependency of heat transfer coefficient on dof type
                if ( mMasterProp( tEmissivityIndex )->check_dof_dependency( tDepDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tResStartIndex, tResStopIndex },
                            { tDepStartIndex, tDepStopIndex } ) +=
                                    aWStar * mStefanBoltzmannConst *
                                    ( std::pow( tT - tT0 , 4.0 ) - std::pow( tTinf - tT0 , 4.0 ) ) *
                                    trans( tFI->N() ) * mMasterProp( tEmissivityIndex )->dPropdDOF( tDepDofType );
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Radiation::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, " IWG_Diffusion_Radiation::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Radiation::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Radiation::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
