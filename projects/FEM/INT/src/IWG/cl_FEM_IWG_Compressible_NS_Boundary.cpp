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

        void IWG_Compressible_NS_Boundary::reset_child_eval_flags()
        {
            // reset eval flags
            mTractionEval = true;
            mTractionDofEval = true;
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Boundary::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ), 
                    "IWG_Compressible_NS_Boundary::compute_residual() - Only pressure or density primitive variables supported for now." );

            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterRes1StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // compute the residual
            mSet->get_residual()( 0 )( { tMasterRes1StartIndex, tMasterRes3StopIndex }, { 0, 0 } ) -= aWStar * (
                    trans( this->W() ) * this->Traction() );

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
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ), 
                    "IWG_Compressible_NS_Boundary::compute_jacobian() - Only pressure or density primitive variables supported for now." );

            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterRes1StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // check DoF dependencies
            MORIS_ASSERT( check_dof_dependencies( mSet, mResidualDofType, mRequestedMasterGlobalDofTypes ), 
                    "IWG_Compressible_NS_Bulk::compute_jacobian - Set of DoF dependencies not suppported." );

            // get the indeces for assembly
            sint tDofFirstDepIndex     = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 0 )( 0 ), mtk::Master_Slave::MASTER );
            sint tDofThirdDepIndex     = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDep1StartIndex = mSet->get_jac_dof_assembly_map()( tMasterDof1Index )( tDofFirstDepIndex, 0 );
            uint tMasterDep3StopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDof3Index )( tDofThirdDepIndex, 1 );                

            // get subview of jacobian for += operations   
            auto tJac = mSet->get_jacobian()( { tMasterRes1StartIndex, tMasterRes3StopIndex }, { tMasterDep1StartIndex, tMasterDep3StopIndex } );  

            // compute jacobian
            tJac -= aWStar * ( trans( this->W() ) * this->dTractiondDOF() );

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
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Boundary::Traction()
        {
            // check if the A matrices have already been assembled
            if( !mTractionEval )
            {
                return mTraction;
            }  
            
            // set the eval flag
            mTractionEval = false;  

            // number of state variables and total bases
            uint tNumStateVars = this->num_space_dims() + 2; 

            // initialize traction
            mTraction.set_size( tNumStateVars, 1, 0.0 );
            auto tTraction2 = mTraction( { 1, tNumStateVars - 2 }, { 0, 0 } );
            auto tTraction3 = mTraction( { tNumStateVars - 1, tNumStateVars - 1 }, { 0, 0 } );

            // get the velocity, FIXME: this needs to be done through the material model for conservative and entropy variables
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( this->get_primary_state_var( 1 ) );

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the properties
            std::shared_ptr< Property > tPropHeatFlux = mMasterProp( static_cast< uint >( IWG_Property_Type::HEAT_FLUX ) );
            std::shared_ptr< Property > tPropTraction = mMasterProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );
            std::shared_ptr< Property > tPropPressure = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) );

            // check if a traction is prescribed and use it, if not compute it
            if ( tPropTraction != nullptr )
            {
                tTraction2 += tPropTraction->val().matrix_data();
                tTraction3 += tFIVelocity->val_trans() * tPropTraction->val();
            }
            else
            {
                tTraction2 += tCM->traction( mNormal, CM_Function_Type::MECHANICAL ).matrix_data();
                tTraction3 += tCM->traction( mNormal, CM_Function_Type::WORK ).matrix_data();
            }

            // check if a heat flux is prescribed and use it
            if ( tPropHeatFlux != nullptr )
            {
                tTraction3 -= tPropHeatFlux->val().matrix_data();
            }

            // check if a pressure is prescribed and apply it
            if ( tPropPressure != nullptr )
            {
                tTraction2 += ( tMM->pressure()( 0 ) - tPropPressure->val()( 0 ) ) * mNormal;
                tTraction3 += ( tMM->pressure()( 0 ) - tPropPressure->val()( 0 ) ) * tFIVelocity->val_trans() * mNormal;
            }

            // return 
            return mTraction;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Boundary::dTractiondDOF()
        {
            // check if the A matrices have already been assembled
            if( !mTractionDofEval )
            {
                return mTractionDOF;
            }  
            
            // set the eval flag
            mTractionDofEval = false;  

            // get the velocity, FIXME: this needs to be done through the material model for conservative and entropy variables
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( this->get_primary_state_var( 1 ) );

            // get the properties
            std::shared_ptr< Property > tPropMu = mMasterProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropKappa = mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // number of state variables and total bases
            uint tNumStateVars = this->num_space_dims() + 2;
            uint tTotNumBases = tNumStateVars * this->num_bases();

            // initialize traction
            mTractionDOF.set_size( tNumStateVars, tTotNumBases , 0.0 );

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
                uint tMasterDofIndex = mSet->get_dof_index_for_type( this->get_primary_state_var( iDOF ), mtk::Master_Slave::MASTER );
                uint tDepDofIndex     = mSet->get_dof_index_for_type( tDepDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex   = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDepDofIndex, 0 );
                uint tDepStopIndex    = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDepDofIndex, 1 );

                // get matrix subviews for different residuals - FIXME: assuming three different Residual DoF-Types
                auto tTraction2DOF = mTractionDOF( { 1, tNumStateVars - 2 }, { tDepStartIndex, tDepStopIndex } );
                auto tTraction3DOF = mTractionDOF( { tNumStateVars - 1, tNumStateVars - 1 }, { tDepStartIndex, tDepStopIndex } );

                // check if a traction is prescribed and use it, if not compute it
                if ( tPropTraction != nullptr )
                {
                    if ( tPropTraction->check_dof_dependency( tDepDofType ) )
                    {
                        tTraction2DOF += tPropTraction->dPropdDOF( tDepDofType ).matrix_data();
                        tTraction3DOF += tFIVelocity->val_trans() * tPropTraction->dPropdDOF( tDepDofType );
                    }

                    if ( tDepDofType( 0 ) == this->get_primary_state_var( 1 ) )
                    {
                        tTraction3DOF += trans( tPropTraction->val() ) * tFIVelocity->N();
                    }
                }
                else if ( tCM->check_dof_dependency( tDepDofType ) )
                {
                    tTraction2DOF += tCM->dTractiondDOF( tDepDofType, mNormal, CM_Function_Type::MECHANICAL ).matrix_data();
                    tTraction3DOF += tCM->dTractiondDOF( tDepDofType, mNormal, CM_Function_Type::WORK ).matrix_data();
                }

                // check if a heat flux is prescribed and use it
                if ( ( tPropHeatFlux != nullptr ) and ( tPropHeatFlux->check_dof_dependency( tDepDofType ) ) )
                {
                    tTraction3DOF -= tPropHeatFlux->dPropdDOF( tDepDofType ).matrix_data();
                }

                // check if a pressure is prescribed and apply it
                if ( tPropPressure != nullptr )
                {
                    if ( tPropPressure->check_dof_dependency( tDepDofType ) )
                    {
                        tTraction2DOF -= mNormal * tPropPressure->dPropdDOF( tDepDofType );
                        tTraction3DOF -= tFIVelocity->val_trans() * mNormal * tPropPressure->dPropdDOF( tDepDofType );
                    }
                    
                    if ( tMM->check_dof_dependency( tDepDofType ) )
                    {
                        tTraction2DOF += mNormal * tMM->PressureDOF( tDepDofType );
                        tTraction3DOF += tFIVelocity->val_trans() * mNormal * tMM->PressureDOF( tDepDofType );
                    }

                    if ( tDepDofType( 0 ) == this->get_primary_state_var( 1 ) )
                    {
                        tTraction3DOF += ( tMM->pressure()( 0 ) - tPropPressure->val()( 0 ) ) * trans( mNormal) * tFIVelocity->N();
                    }
                }
            } // end loop over dof dependencies

            // return 
            return mTractionDOF;
        }

        //------------------------------------------------------------------------------
        
    } /* namespace fem */
} /* namespace moris */
