/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Radiation.cpp
 *
 */

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
            mPropertyMap[ "Emissivity" ]         = static_cast< uint >( IWG_Property_Type::EMISSIVITY );
            mPropertyMap[ "AmbientTemperature" ] = static_cast< uint >( IWG_Property_Type::AMBIENT_TEMP );
            mPropertyMap[ "AbsoluteZero" ]       = static_cast< uint >( IWG_Property_Type::ABSOLUTE_ZERO );
            mPropertyMap[ "Thickness" ]          = static_cast< uint >( IWG_Property_Type::THICKNESS );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Radiation::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif
            // get emissivity property
            const std::shared_ptr< Property >& tPropEmissivity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::EMISSIVITY ) );

            // get ambient temperature property
            const std::shared_ptr< Property >& tPropAmbientTemp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::AMBIENT_TEMP ) );

            // get reference temperature property
            const std::shared_ptr< Property >& tPropAbsoluteZero =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::ABSOLUTE_ZERO ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get filed interpolator for residual dof type
            Field_Interpolator* tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get property values
            real tT0    = tPropAbsoluteZero->val()( 0 );
            real tTinf  = tPropAmbientTemp->val()( 0 );
            real tAlpha = tPropEmissivity->val()( 0 );

            // current temperature
            real tT = tFI->val()( 0 );

            // compute the residual
            // N * a * (T - T_ref)
            mSet->get_residual()( 0 )(
                    { tResStartIndex, tResStopIndex } ) +=                                    //
                    aWStar * (                                                                //
                            mStefanBoltzmannConst * tAlpha *                                  //
                            ( std::pow( tT - tT0, 4.0 ) - std::pow( tTinf - tT0, 4.0 ) ) *    //
                            tFI->N_trans() );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Radiation::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Radiation::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif
            // get emissivity property
            const std::shared_ptr< Property >& tPropEmissivity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::EMISSIVITY ) );

            // get ambient temperature property
            const std::shared_ptr< Property >& tPropAmbientTemp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::AMBIENT_TEMP ) );

            // get reference temperature property
            const std::shared_ptr< Property >& tPropAbsoluteZero =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::ABSOLUTE_ZERO ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator* tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get property values
            real tT0    = tPropAbsoluteZero->val()( 0 );
            real tTinf  = tPropAmbientTemp->val()( 0 );
            real tAlpha = tPropEmissivity->val()( 0 );

            // current temperature
            real tT = tFI->val()( 0 );

            // compute the jacobian for dof dependencies
            for ( uint iDOF = 0; iDOF < mRequestedMasterGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                const Cell< MSI::Dof_Type >& tDepDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the dof type indices for assembly
                uint tDepDofIndex   = mSet->get_dof_index_for_type( tDepDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tResStartIndex, tResStopIndex },
                        { tDepStartIndex, tDepStopIndex } );

                // if dof type is residual dof type
                if ( tDepDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJac += aWStar * tAlpha * mStefanBoltzmannConst *                                   //
                            ( 4.0 * std::pow( tT, 3.0 ) - 12.0 * tT0 * std::pow( tT, 2.0 )              //
                                    + 2.0 * tT * std::pow( tT0, 2.0 ) + 4.0 * std::pow( tT0, 3.0 ) )    //
                          * tFI->N_trans() * tFI->N();
                }

                // if dependency of heat transfer coefficient on dof type
                if ( tPropEmissivity->check_dof_dependency( tDepDofType ) )
                {
                    // add contribution to jacobian
                    tJac += aWStar * mStefanBoltzmannConst *                                  //
                            ( std::pow( tT - tT0, 4.0 ) - std::pow( tTinf - tT0, 4.0 ) ) *    //
                            tFI->N_trans() * tPropEmissivity->dPropdDOF( tDepDofType );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Diffusion_Neumann::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Radiation::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, " IWG_Diffusion_Radiation::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Radiation::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Radiation::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
