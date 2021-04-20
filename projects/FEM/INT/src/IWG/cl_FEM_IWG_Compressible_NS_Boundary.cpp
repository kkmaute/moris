/*
 * cl_FEM_IWG_Compressible_NS_Boundary.cpp
 *
 *  Created on: Mar 16, 2021
 *      Author: wunsch
 */

#include "cl_FEM_IWG_Compressible_NS_Boundary.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_FEM_IWG_Compressible_NS.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Compressible_NS_Boundary::IWG_Compressible_NS_Boundary()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "HeatFlux" ] = static_cast< uint >( IWG_Property_Type::HEAT_FLUX );
            mPropertyMap[ "Traction" ] = static_cast< uint >( IWG_Property_Type::TRACTION );
            mPropertyMap[ "Pressure" ] = static_cast< uint >( IWG_Property_Type::PRESSURE );

            // set size for the material model pointer cell
            mMasterMM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the material map
            mMaterialMap[ "FluidMM" ] = static_cast< uint >( IWG_Material_Type::FLUID_MM );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "FluidCM" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_CM );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Boundary::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType( 0 )  ), 
                    "IWG_Compressible_NS_Boundary::compute_residual() - Only pressure or density primitive variables supported for now." );

            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            //uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof2Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 1 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 2 ), mtk::Master_Slave::MASTER );
            //uint tMasterRes1StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 0 );
            //uint tMasterRes1StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 1 );
            uint tMasterRes2StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof2Index )( 0, 0 );
            uint tMasterRes2StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof2Index )( 0, 1 );
            uint tMasterRes3StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // get the FIs associated with each residual dof type
            //Field_Interpolator * tFIFirstDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFIThirdDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 2 ) );

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the properties
            std::shared_ptr< Property > tPropHeatFlux = mMasterProp( static_cast< uint >( IWG_Property_Type::HEAT_FLUX ) );
            std::shared_ptr< Property > tPropTraction = mMasterProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );
            std::shared_ptr< Property > tPropPressure = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) );

            // get matrix subviews for different residuals - FIXME: assuming three different Residual DoF-Types
            //auto tRes1 = mSet->get_residual()( 0 )( { tMasterRes1StartIndex, tMasterRes1StopIndex }, { 0, 0 } );
            auto tRes2 = mSet->get_residual()( 0 )( { tMasterRes2StartIndex, tMasterRes2StopIndex }, { 0, 0 } );
            auto tRes3 = mSet->get_residual()( 0 )( { tMasterRes3StartIndex, tMasterRes3StopIndex }, { 0, 0 } );

            // check if a traction is prescribed and use it, if not compute it
            if ( tPropTraction != nullptr )
            {
                tRes2 -= aWStar * ( tFIVelocity->N_trans() * tPropTraction->val() );
                tRes3 -= aWStar * ( tFIThirdDofType->N_trans() * tFIVelocity->val_trans() * tPropTraction->val() );
            }
            else
            {
                tRes2 -= aWStar * ( tFIVelocity->N_trans() * tCM->traction( mNormal, CM_Function_Type::MECHANICAL ) );
                tRes3 -= aWStar * ( tFIThirdDofType->N_trans() * tCM->traction( mNormal, CM_Function_Type::WORK ) );
            }

            // check if a heat flux is prescribed and use it
            if ( tPropHeatFlux != nullptr )
            {
                tRes3 += aWStar * ( tFIThirdDofType->N_trans() * tPropHeatFlux->val() );
            }

            // check if a pressure is prescribed and apply it
            if ( tPropPressure != nullptr )
            {
                tRes2 -= aWStar * ( tFIVelocity->N_trans() * 
                        ( tMM->pressure()( 0 ) - tPropPressure->val()( 0 ) ) * mNormal );
                tRes3 -= aWStar * ( tFIThirdDofType->N_trans() * 
                        ( tMM->pressure()( 0 ) - tPropPressure->val()( 0 ) ) * tFIVelocity->val_trans() * mNormal );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Boundary::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Boundary::compute_jacobian( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType( 0 )  ), 
                    "IWG_Compressible_NS_Boundary::compute_jacobian() - Only pressure or density primitive variables supported for now." );

            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof2Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 1 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 2 ), mtk::Master_Slave::MASTER );
            uint tMasterRes2StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof2Index )( 0, 0 );
            uint tMasterRes2StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof2Index )( 0, 1 );
            uint tMasterRes3StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // get the FIs associated with each residual dof type
            //Field_Interpolator * tFIFirstDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFIThirdDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 2 ) );

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the properties
            std::shared_ptr< Property > tPropHeatFlux = mMasterProp( static_cast< uint >( IWG_Property_Type::HEAT_FLUX ) );
            std::shared_ptr< Property > tPropTraction = mMasterProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );
            std::shared_ptr< Property > tPropPressure = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) );

            // compute the jacobian for dof dependencies
            for( uint iDOF = 0; iDOF < mRequestedMasterGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDepDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the dof type indices for assembly
                uint tDepDofIndex    = mSet->get_dof_index_for_type( tDepDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex  = mSet->get_jac_dof_assembly_map()( tMasterDof1Index )( tDepDofIndex, 0 );
                uint tDepStopIndex   = mSet->get_jac_dof_assembly_map()( tMasterDof1Index )( tDepDofIndex, 1 );

                // get matrix subviews for different residuals - FIXME: assuming three different Residual DoF-Types
                //auto tJac1 = mSet->get_jacobian()( { tMasterRes1StartIndex, tMasterRes1StopIndex }, { tDepStartIndex, tDepStopIndex } );
                auto tJac2 = mSet->get_jacobian()( { tMasterRes2StartIndex, tMasterRes2StopIndex }, { tDepStartIndex, tDepStopIndex } );
                auto tJac3 = mSet->get_jacobian()( { tMasterRes3StartIndex, tMasterRes3StopIndex }, { tDepStartIndex, tDepStopIndex } );

                // check if a traction is prescribed and use it, if not compute it
                if ( tPropTraction != nullptr )
                {
                    if ( tPropTraction->check_dof_dependency( tDepDofType ) )
                    {
                        tJac2 -= aWStar * ( tFIVelocity->N_trans() * tPropTraction->dPropdDOF( tDepDofType ) );
                        tJac3 -= aWStar * ( tFIThirdDofType->N_trans() * tFIVelocity->val_trans() * tPropTraction->dPropdDOF( tDepDofType ) );
                    }

                    if ( tDepDofType( 0 ) == mDofVelocity )
                    {
                        tJac3 -= aWStar * ( tFIThirdDofType->N_trans() * trans( tPropTraction->val() ) * tFIVelocity->N() );
                    }
                }
                else if ( tCM->check_dof_dependency( tDepDofType ) )
                {
                    tJac2 -= aWStar * ( tFIVelocity->N_trans() * tCM->dTractiondDOF( tDepDofType, mNormal, CM_Function_Type::MECHANICAL ) );
                    tJac3 -= aWStar * ( tFIThirdDofType->N_trans() * tCM->dTractiondDOF( tDepDofType, mNormal, CM_Function_Type::WORK ) );
                }

                // check if a heat flux is prescribed and use it
                if ( ( tPropHeatFlux != nullptr ) and ( tPropHeatFlux->check_dof_dependency( tDepDofType ) ) )
                {
                    tJac3 += aWStar * ( tFIThirdDofType->N_trans() * tPropHeatFlux->dPropdDOF( tDepDofType ) );
                }

                // check if a pressure is prescribed and apply it
                if ( tPropPressure != nullptr )
                {
                    if ( tPropPressure->check_dof_dependency( tDepDofType ) )
                    {
                        tJac2 -= aWStar * ( tFIVelocity->N_trans() * mNormal * tPropPressure->dPropdDOF( tDepDofType ) );
                        tJac3 -= aWStar * ( tFIThirdDofType->N_trans() * ( tFIVelocity->val_trans() * mNormal ) * tPropPressure->dPropdDOF( tDepDofType ) );
                    }
                    
                    if ( tMM->check_dof_dependency( tDepDofType ) )
                    {
                        tJac2 -= aWStar * ( tFIVelocity->N_trans() * mNormal * tMM->PressureDOF( tDepDofType ) );
                        tJac3 -= aWStar * ( tFIThirdDofType->N_trans() * ( tFIVelocity->val_trans() * mNormal ) * tMM->PressureDOF( tDepDofType ) );
                    }

                    if ( tDepDofType( 0 ) == mDofVelocity )
                    {
                        tJac3 -= aWStar * ( tFIThirdDofType->N_trans() * 
                                ( tMM->pressure()( 0 ) - tPropPressure->val()( 0 ) ) * trans( mNormal) * tFIVelocity->N() );
                    }
                }
            } // end loop over dof dependencies

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Boundary::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Boundary::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, " IWG_Compressible_NS_Boundary::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Boundary::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Boundary::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
        
    } /* namespace fem */
} /* namespace moris */
