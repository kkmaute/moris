/*
 * cl_FEM_IWG_Compressible_NS_Bulk_Flux_Matrices.cpp
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

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::assemble_variable_set()
        {
            // check if the variable vectors have already been assembled
            if( !mVarVecEval )
            {      
                return;
            }         
            
            // set the eval flag
            mVarVecEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ), 
                    "IWG_Compressible_NS_Bulk::assemble_variable_set() - check for residual DoF types failed." );

            //FIXME: only density and pressure primitive variable sets supported in this function

            // get different field interpolators
            Field_Interpolator * tFIFirstVar =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFIThirdVar =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 2 ) );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // assemble field var vector and matrix based on number of spatial dimensions

            // build vector of field variables
            mY.set_size( tNumSpaceDims + 2, 1, 0.0 );
            mY( 0 ) = tFIFirstVar->val()( 0 );
            Matrix< DDRMat > tVelVec = tFIVelocity->val();
            mY( { 1, tNumSpaceDims } ) = tVelVec.matrix_data();
            mY( tNumSpaceDims + 1 ) = tFIThirdVar->val()( 0 );

            // build vector of the time derivatives of the field variables
            mdYdt.set_size( tNumSpaceDims + 2, 1, 0.0 );
            mdYdt( 0 ) = tFIFirstVar->gradt( 1 )( 0 );
            Matrix< DDRMat > tdVeldt = trans( tFIVelocity->gradt( 1 ) );
            mdYdt( { 1, tNumSpaceDims } ) = tdVeldt.matrix_data();
            mdYdt( tNumSpaceDims + 1 ) = tFIThirdVar->gradt( 1 )( 0 );

            // build matrix of spatial derivatives of field variables
            mdYdx.set_size( tNumSpaceDims + 2, tNumSpaceDims, 0.0 );
            Matrix< DDRMat > tGradP = trans( tFIFirstVar->gradx( 1 ) );
            mdYdx( { 0, 0 }, { 0, tNumSpaceDims - 1 } ) = tGradP.matrix_data();
            Matrix< DDRMat > tGradVel = trans( tFIVelocity->gradx( 1 ) );
            mdYdx( { 1, tNumSpaceDims }, { 0, tNumSpaceDims - 1 } ) = tGradVel.matrix_data();
            Matrix< DDRMat > tGradT = trans( tFIThirdVar->gradx( 1 ) );
            mdYdx( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 0, tNumSpaceDims - 1 } ) = tGradT.matrix_data();
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::assemble_variable_DOF_set()
        {  
            // check if the variable vectors have already been assembled
            if( !mVarDofVecEval )
            {
                return;
            }    
            
            // set the eval flag
            mVarDofVecEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ), 
                    "IWG_Compressible_NS_Bulk::assemble_variable_DOF_set() - check for residual DoF types failed." );

            //FIXME: only density and pressure primitive variable sets supported so far

            // get different field interpolators
            Field_Interpolator * tFIFirstVar = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            Field_Interpolator * tFIVelocity = mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFIThirdVar = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 2 ) );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();              

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

            // assemble field var vector and matrix based on number of spatial dimensions

            // initialize the cell and matrices 
            mYDOF.set_size( tNumSpaceDims + 2, ( tNumSpaceDims + 2 ) * tNumBases, 0.0 );
            mdYdtDOF.set_size( tNumSpaceDims + 2, ( tNumSpaceDims + 2 ) * tNumBases, 0.0 );
            mdYdxDOF.assign( tNumSpaceDims, mYDOF );

            // build vector of field variables
            mYDOF( { 0, 0 }, { 0, tNumBases - 1 } ) = tFIFirstVar->N().matrix_data();
            mYDOF( { 1, tNumSpaceDims }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) = tFIVelocity->N().matrix_data();
            mYDOF( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = tFIThirdVar->N().matrix_data();

            Matrix< DDRMat > tdNveldt;
            this->compute_dnNdtn( tdNveldt );
            mdYdtDOF( { 0, 0 }, { 0, tNumBases - 1 } ) = tFIFirstVar->dnNdtn( 1 ).matrix_data();
            mdYdtDOF( { 1, tNumSpaceDims }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) = tdNveldt.matrix_data();
            mdYdtDOF( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tFIThirdVar->dnNdtn( 1 ).matrix_data();

            for ( uint iDim = 0; iDim < tNumSpaceDims; iDim++ )
            {
                mdYdxDOF( iDim )( { 0, 0 }, { 0, tNumBases - 1 } ) = tFIFirstVar->dnNdxn( 1 )( { iDim, iDim },{ 0, tNumBases - 1 } );
                mdYdxDOF( iDim )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tFIVelocity->dnNdxn( 1 )( { iDim, iDim },{ 0, tNumBases - 1 } );
                mdYdxDOF( iDim )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tFIVelocity->dnNdxn( 1 )( { iDim, iDim },{ 0, tNumBases - 1 } );
                if ( tNumSpaceDims == 3 )
                {
                    mdYdxDOF( iDim )( { 3, 3 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = tFIVelocity->dnNdxn( 1 )( { iDim, iDim },{ 0, tNumBases - 1 } );
                }
                mdYdxDOF( iDim )( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tFIThirdVar->dnNdxn( 1 )( { iDim, iDim },{ 0, tNumBases - 1 } );
            } 
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::eval_A_matrices()
        {
            // check if the variable vectors have already been assembled
            if( !mFluxAMatEval )
            {
                return;
            }  
            
            // set the eval flag
            mFluxAMatEval = false;          

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // evaluate A-matrices and store them
            eval_A( tMM, tCM, mMasterFIManager, mResidualDofType, mA );  
        }

       //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::eval_A_DOF_matrices()
        {
            // check if the A-Dof flux matrices have already been assembled
            if( !mFluxADofMatEval )
            {
                return;
            }
            
            // set the eval flag
            mFluxADofMatEval = false;

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

            // create standard empty matrix
            const Matrix< DDRMat > tEmptyADof( tNumSpaceDims + 2, ( tNumSpaceDims + 2 ) * tNumBases, 0.0 );

            // initialize derivatives of A-matrices
            mADOF.resize( tNumSpaceDims + 1 );
            for ( uint iMat = 0; iMat < tNumSpaceDims + 1; iMat++)
            {
                mADOF( iMat ).assign( tNumSpaceDims + 2, tEmptyADof );
            }      

            // evaluate the derivatives for each of the matrices and store them
            eval_A0_DOF( tMM, tCM, mMasterFIManager, mResidualDofType, mADOF( 0 ) ); 
            eval_A1_DOF( tMM, tCM, mMasterFIManager, mResidualDofType, mADOF( 1 ) );  
            eval_A2_DOF( tMM, tCM, mMasterFIManager, mResidualDofType, mADOF( 2 ) ); 
            if ( tNumSpaceDims == 3 )
            {
                eval_A3_DOF( tMM, tCM, mMasterFIManager, mResidualDofType, mADOF( 3 ) );  
            } 
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::eval_K_matrices()
        {
            // check if K matrices have already been evaluated
            if ( !mFluxKMatEval )
            {
                return;
            }

            // set the eval flag
            mFluxKMatEval = false;            

            // get the viscosity
            std::shared_ptr< Property > tPropDynamicViscosity = mMasterProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropThermalConductivity = mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // eval K matrices and store them
            eval_K( tPropDynamicViscosity, tPropThermalConductivity, mMasterFIManager, mK );
        }                    

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
