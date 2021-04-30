/*
 * cl_FEM_IWG_Compressible_NS_Bulk.cpp
 *
 *  Created on: Feb 10, 2021
 *      Author: wunsch
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Bulk.hpp"
#include "fn_FEM_IWG_Compressible_NS.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

// debug - output to hdf5
#include "paths.hpp"
#include "HDF5_Tools.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Compressible_NS_Bulk::IWG_Compressible_NS_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "DynamicViscosity" ]     = static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY );
            mPropertyMap[ "ThermalConductivity" ]  = static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY );
            mPropertyMap[ "BodyForce" ]            = static_cast< uint >( IWG_Property_Type::BODY_FORCE );
            mPropertyMap[ "BodyHeatLoad" ]         = static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD );

            // set size for the material model pointer cell
            mMasterMM.resize( static_cast< uint >( IWG_Material_Type::MAX_ENUM ), nullptr );

            // populate the material map
            mMaterialMap[ "FluidMM" ] = static_cast< uint >( IWG_Material_Type::FLUID_MM );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "FluidCM" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_CM );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the stabilization parameter map
            mStabilizationMap[ "GLS" ] = static_cast< uint >( IWG_Stabilization_Type::GLS );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::reset_child_eval_flags()
        {
            // reset eval flags specific to this child IWG
            mLYEval = true;
            mLWEval = true;
            mLDofYEval = true;

            mA0invEval = true;
            mdA0invdYEval = true;

            mGEval = true;

            mMEval = true;
            mMinvEval = true;
            mdMdYEval = true;

            mTauEval = true;
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ), 
                    "IWG_Compressible_NS_Bulk::compute_jacobian() - Only pressure or density primitive variables supported for now." );


            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterRes1StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // get number of space dimensions
            uint tNumSpaceDims = this->num_space_dims();

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get subview for complete residual
            auto tRes = mSet->get_residual()( 0 )( { tMasterRes1StartIndex, tMasterRes3StopIndex }, { 0, 0 } );

            // A0 matrix contribution
            tRes += aWStar * trans( this->W() ) * this->A( 0 ) * this->dYdt();

            // loop over A-Matrices
            for ( uint iA = 1; iA < tNumSpaceDims + 1; iA++ )
            {
                // compute residual
                tRes += aWStar * trans( this->W() ) * this->A( iA ) * this->dYdx( iA - 1 );
            }

            // loop over K-matrices
            for ( uint iDim = 0; iDim < tNumSpaceDims; iDim++ )
            {
                for ( uint jDim = 0; jDim < tNumSpaceDims; jDim++ )
                {
                    // compute residual
                    tRes += aWStar * trans( this->dWdx( iDim ) ) * this->K( iDim, jDim ) * this->dYdx( jDim );
                }
            }

            // get the Stabilization Parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GLS ) );

            // add contribution of stabilization term if stabilization parameter has been set
            if ( tSP != nullptr )
            {
                // stabilization term for testing
                //tRes += aWStar * tSP->val()( 0 ) * trans( this->W() ) * this->LY();

                // debug - stabilization term without Tau for testing
                // tRes += aWStar * tSP->val()( 0 ) * trans( this->LW() ) * this->LY();

                // debug - test M
                uint tNumStateVars = tNumSpaceDims + 2;
                Matrix< DDRMat > tVR( tNumStateVars, 1, 2.3 );
                tRes += aWStar * tSP->val()( 0 ) * trans( this->W() ) * this->M() * tVR;

                // GLS stabilization term
                // tRes += aWStar * tSP->val()( 0 ) * trans( this->LW() ) * this->Tau() * this->LY();
            }           

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Bulk::compute_residual - Residual contains NAN or INF, exiting!");                                 
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_jacobian( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ), 
                    "IWG_Compressible_NS_Bulk::compute_jacobian() - Only pressure or density primitive variables supported for now." );

            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterRes1StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // get number of space dimensions
            uint tNumSpaceDims = this->num_space_dims();                   

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

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

            // add contribution from d(A0)/dDof * Y,t
            Matrix< DDRMat > tdAdY;
            eval_dAdY_VR( tMM, tCM, mMasterFIManager, mResidualDofType, this->dYdt(), 0, tdAdY );
            tJac += aWStar * trans( this->W() ) * ( tdAdY * this->W() + this->A( 0 ) * this->dWdt() ); 

            // loop over contributions from A-matrices
            for ( uint iDim = 0; iDim < tNumSpaceDims; iDim++ )
            {
                // evaluate d(Ai)/dDof * Y,i
                eval_dAdY_VR( tMM, tCM, mMasterFIManager, mResidualDofType, this->dYdx( iDim ), iDim + 1, tdAdY );

                // add contribution
                tJac += aWStar * trans( this->W() ) * ( tdAdY * this->W() + this->A( iDim + 1 ) * this->dWdx( iDim ) );

            }

            // get properties for K-matrices
            std::shared_ptr< Property > tPropMu    = mMasterProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropKappa = mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // loop over contributions from K-matrices
            Matrix< DDRMat > dKdY;
            for ( uint iDim = 0; iDim < tNumSpaceDims; iDim++ )
            {
                for ( uint jDim = 0; jDim < tNumSpaceDims; jDim++ )
                {
                    // get dKij/dY * Y,ij
                    eval_dKdY_VR( tPropMu, tPropKappa, mMasterFIManager, this->dYdx( jDim ), iDim, jDim, dKdY );
                                
                    // add contributions from K-matrices
                    tJac += aWStar * trans( this->dWdx( iDim ) ) * ( dKdY * this->W() + K( iDim, jDim ) * this->dWdx( jDim ) );
                }
            }

            // get the Stabilization Parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GLS ) );

            // add contribution of stabilization term if stabilization parameter has been set
            if ( tSP != nullptr )
            {
                // debug - stabilization term for testing
                // tJac += aWStar * tSP->val()( 0 ) * trans( this->W() ) * ( this->LW() + this->dLdDofY() );

                // debug - stabilization term without Tau for testing
                // tJac += aWStar * tSP->val()( 0 ) * ( 
                //         this->dLdDofW( this->LY() ) + 
                //         trans( this->LW() ) * ( this->LW() + this->dLdDofY() ) );

                // debug - test M
                uint tNumStateVars = tNumSpaceDims + 2;
                Matrix< DDRMat > tVR( tNumStateVars, 1, 2.3 );
                Matrix< DDRMat > tdMVRdY( tNumStateVars, tNumStateVars, 0.0 );
                for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
                {
                    tdMVRdY( { 0, tNumStateVars - 1 }, { iVar, iVar } ) = this->dMdY( iVar ) * tVR;
                }
                tJac += aWStar * tSP->val()( 0 ) * trans( this->W() ) * tdMVRdY * this->W();

                // debug - test GLS terms

                // GLS stabilization term
                // tJac += aWStar * tSP->val()( 0 ) * (
                //         trans( this->LW() ) * ( this->Tau() * this->dLdDofY() + this->dTaudY( this->LY() ) * this->W() ) +
                //         this->dLdDofW( this->Tau() * this->LY() ) );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");                     
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::compute_dRdp - Not implemented." );
        }
        
        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
