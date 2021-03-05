/*
 * cl_FEM_IWG_Compressible_NS_Bulk_Flux_Matrices.cpp
 *
 *  Created on: Feb 10, 2021
 *      Author: wunsch
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Bulk.hpp"

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
            MORIS_ASSERT( this->check_residual_dof_types(), "IWG_Compressible_NS_Bulk::assemble_variable_set() - check for residual DoF types failed." );

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
            MORIS_ASSERT( this->check_residual_dof_types(), "IWG_Compressible_NS_Bulk::assemble_variable_DOF_set() - check for residual DoF types failed." );

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

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            real tRho = tMM->density()( 0 );
            real tUx = tFIVelocity->val()( 0 );
            real tUy = tFIVelocity->val()( 1 );
            real tAlphaP = tMM->AlphaP()( 0 );
            real tBetaT = tMM->BetaT()( 0 );
            real tEtot = tCM->Energy()( 0 );

            // constants computed as in M. Polner's 2005 PhD thesis
            real tE1p = tBetaT * tEtot;
            real tE2p = tE1p + 1.0;
            real tE3p = tEtot + tMM->pressure()( 0 );
            real tE4p = -1.0 * tAlphaP * tEtot + tRho * tMM->Cp()( 0 );         

            // reset A matrices
            Matrix< DDRMat > tEmptyA( tNumSpaceDims + 2, tNumSpaceDims + 2, 0.0 );
            mA.assign( tNumSpaceDims + 2, tEmptyA );
            //mA.resize( tNumSpaceDims + 2 );          

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // evaluate A_0
                    mA( 0 ) = { 
                        {       tRho * tBetaT,        0.0,        0.0,       -tRho * tAlphaP }, 
                        { tUx * tRho * tBetaT,       tRho,        0.0, tUx * -tRho * tAlphaP }, 
                        { tUy * tRho * tBetaT,        0.0,       tRho, tUy * -tRho * tAlphaP }, 
                        {                tE1p, tUx * tRho, tUy * tRho,                  tE4p } };

                    // evaluate A_1
                    mA( 1 ) = { 
                        {             tRho * tBetaT * tUx,                    tRho,              0.0,       -tRho * tAlphaP * tUx }, 
                        { 1.0 + tUx * tRho * tBetaT * tUx,        2.0 * tUx * tRho,              0.0, tUx * -tRho * tAlphaP * tUx }, 
                        {       tUy * tRho * tBetaT * tUx,              tUy * tRho,       tUx * tRho, tUy * -tRho * tAlphaP * tUx },  
                        {                      tE2p * tUx, tE3p + tUx * tUx * tRho, tUy * tUx * tRho,                  tE4p * tUx } };

                    // evaluate A_2
                    mA( 2 ) = { 
                        {             tRho * tBetaT * tUy,              0.0,                    tRho,       -tRho * tAlphaP * tUy }, 
                        {       tUx * tRho * tBetaT * tUy,       tUy * tRho,              tUx * tRho, tUx * -tRho * tAlphaP * tUy }, 
                        { 1.0 + tUy * tRho * tBetaT * tUy,              0.0,        2.0 * tUy * tRho, tUy * -tRho * tAlphaP * tUy }, 
                        {                      tE2p * tUy, tUx * tUy * tRho, tE3p + tUy * tUy * tRho,                  tE4p * tUy } };

                    break;
                }

                // for 3D
                case 3 :
                {
                    // get the z-velocity
                    real tUz = tFIVelocity->val()( 2 );

                    // evaluate A_0
                    mA( 0 ) = { 
                        {       tRho * tBetaT,        0.0,        0.0,        0.0,       -tRho * tAlphaP }, 
                        { tUx * tRho * tBetaT,       tRho,        0.0,        0.0, tUx * -tRho * tAlphaP }, 
                        { tUy * tRho * tBetaT,        0.0,       tRho,        0.0, tUy * -tRho * tAlphaP }, 
                        { tUz * tRho * tBetaT,        0.0,        0.0,       tRho, tUz * -tRho * tAlphaP }, 
                        {                tE1p, tUx * tRho, tUy * tRho, tUz * tRho,                  tE4p } };      

                    // evaluate A_1
                    mA( 1 ) = { 
                        {             tRho * tBetaT * tUx,                    tRho,              0.0,              0.0,       -tRho * tAlphaP * tUx }, 
                        { 1.0 + tUx * tRho * tBetaT * tUx,        2.0 * tUx * tRho,              0.0,              0.0, tUx * -tRho * tAlphaP * tUx }, 
                        {       tUy * tRho * tBetaT * tUx,              tUy * tRho,       tUx * tRho,              0.0, tUy * -tRho * tAlphaP * tUx }, 
                        {       tUz * tRho * tBetaT * tUx,              tUz * tRho,              0.0,       tUx * tRho, tUz * -tRho * tAlphaP * tUx }, 
                        {                      tE2p * tUx, tE3p + tUx * tUx * tRho, tUy * tUx * tRho, tUz * tUx * tRho,                  tE4p * tUx } };

                    // evaluate A_2
                    mA( 2 ) = { 
                        {             tRho * tBetaT * tUy,              0.0,                    tRho,              0.0,       -tRho * tAlphaP * tUy }, 
                        {       tUx * tRho * tBetaT * tUy,       tUy * tRho,              tUx * tRho,              0.0, tUx * -tRho * tAlphaP * tUy }, 
                        { 1.0 + tUy * tRho * tBetaT * tUy,              0.0,        2.0 * tUy * tRho,              0.0, tUy * -tRho * tAlphaP * tUy }, 
                        {       tUz * tRho * tBetaT * tUy,              0.0,              tUz * tRho,       tUy * tRho, tUz * -tRho * tAlphaP * tUy }, 
                        {                      tE2p * tUy, tUx * tUy * tRho, tE3p + tUy * tUy * tRho, tUz * tUy * tRho,                  tE4p * tUy } };


                    // evaluate A_3
                    mA( 3 ) = { 
                        {             tRho * tBetaT * tUz,              0.0,              0.0,                    tRho,       -tRho * tAlphaP * tUz }, 
                        {       tUx * tRho * tBetaT * tUz,       tUz * tRho,              0.0,              tUx * tRho, tUx * -tRho * tAlphaP * tUz }, 
                        {       tUy * tRho * tBetaT * tUz,              0.0,       tUz * tRho,              tUy * tRho, tUy * -tRho * tAlphaP * tUz }, 
                        { 1.0 + tUz * tRho * tBetaT * tUz,              0.0,              0.0,        2.0 * tUz * tRho, tUz * -tRho * tAlphaP * tUz }, 
                        {                      tE2p * tUz, tUx * tUz * tRho, tUy * tUz * tRho, tE3p + tUz * tUz * tRho,                  tE4p * tUz } };                        

                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::assemble_A_matrices() - Number of space dimensions must be 2 or 3" );
                };
            }
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

            // evaluate the derivatives for each of the matrices
            this->eval_A0_DOF();
            this->eval_A1_DOF();
            this->eval_A2_DOF();
            if ( tNumSpaceDims == 3 )
            {
                this->eval_A3_DOF();
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

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            real tKa = tPropThermalConductivity->val()( 0 );
            real tMu = tPropDynamicViscosity->val()( 0 );
            real tLa = -2.0 * tMu / 3.0;
            real tCh = 2.0 * tMu + tLa;
            real tUx = tFIVelocity->val()( 0 );
            real tUy = tFIVelocity->val()( 1 );

            // set number of K matrices
            mK.resize( tNumSpaceDims );
            for ( uint iDim = 0; iDim < tNumSpaceDims ; iDim++)
            {
                mK( iDim ).resize( tNumSpaceDims );
            }            

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // evaluate K_11
                    mK( 0 )( 0 ) = { 
                        { 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       tCh,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tMu, 0.0 }, 
                        { 0.0, tUx * tCh, tUy * tMu, tKa } };

                    // evaluate K_12
                    mK( 0 )( 1 ) = { 
                        { 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tLa, 0.0 }, 
                        { 0.0,       tMu,       0.0, 0.0 }, 
                        { 0.0, tUy * tMu, tUx * tLa, 0.0 } };

                    // =======================

                    // evaluate K_21
                    mK( 1 )( 0 ) = { 
                        { 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tMu, 0.0 }, 
                        { 0.0,       tLa,       0.0, 0.0 }, 
                        { 0.0, tUy * tLa, tUx * tMu, 0.0 } };

                        // evaluate K_22
                    mK( 1 )( 1 ) = { 
                        { 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       tMu,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tCh, 0.0 }, 
                        { 0.0, tUx * tMu, tUy * tCh, tKa } };

                    break;
                }

                // for 3D
                case 3 :
                {
                    // get the z-velocity
                    real tUz = tFIVelocity->val()( 2 );

                    // evaluate K_11
                    mK( 0 )( 0 ) = { 
                        { 0.0,       0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       tCh,       0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tMu,       0.0, 0.0 }, 
                        { 0.0,       0.0,       0.0,       tMu, 0.0 }, 
                        { 0.0, tUx * tCh, tUy * tMu, tUz * tMu, tKa } };

                    // evaluate K_12
                    mK( 0 )( 1 ) = { 
                        { 0.0,       0.0,       0.0, 0.0, 0.0 }, 
                        { 0.0,       0.0,       tLa, 0.0, 0.0 }, 
                        { 0.0,       tMu,       0.0, 0.0, 0.0 }, 
                        { 0.0,       0.0,       0.0, 0.0, 0.0 }, 
                        { 0.0, tUy * tMu, tUx * tLa, 0.0, 0.0 } };

                    // evaluate K_13
                    mK( 0 )( 2 ) = { 
                        { 0.0,       0.0, 0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0, 0.0,       tLa, 0.0 }, 
                        { 0.0,       0.0, 0.0,       0.0, 0.0 }, 
                        { 0.0,       tMu, 0.0,       0.0, 0.0 }, 
                        { 0.0, tUz * tMu, 0.0, tUx * tLa, 0.0 } };

                    // =======================

                    // evaluate K_21
                    mK( 1 )( 0 ) = { 
                        { 0.0,       0.0,       0.0, 0.0, 0.0 }, 
                        { 0.0,       0.0,       tMu, 0.0, 0.0 }, 
                        { 0.0,       tLa,       0.0, 0.0, 0.0 }, 
                        { 0.0,       0.0,       0.0, 0.0, 0.0 }, 
                        { 0.0, tUy * tLa, tUx * tMu, 0.0, 0.0 } };

                        // evaluate K_22
                    mK( 1 )( 1 ) = { 
                        { 0.0,       0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       tMu,       0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tCh,       0.0, 0.0 }, 
                        { 0.0,       0.0,       0.0,       tMu, 0.0 }, 
                        { 0.0, tUx * tMu, tUy * tCh, tUz * tMu, tKa } };

                        // evaluate K_23
                    mK( 1 )( 2 ) = { 
                        { 0.0, 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0, 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0, 0.0,       0.0,       tLa, 0.0 }, 
                        { 0.0, 0.0,       tMu,       0.0, 0.0 }, 
                        { 0.0, 0.0, tUy * tMu, tUy * tLa, 0.0 } };

                    // =======================

                    // evaluate K_31
                    mK( 2 )( 0 ) = { 
                        { 0.0,       0.0, 0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0, 0.0,       tMu, 0.0 }, 
                        { 0.0,       0.0, 0.0,       0.0, 0.0 }, 
                        { 0.0,       tLa, 0.0,       0.0, 0.0 }, 
                        { 0.0, tUz * tLa, 0.0, tUx * tMu, 0.0 } };

                    // evaluate K_32
                    mK( 2 )( 1 ) = { 
                        { 0.0, 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0, 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0, 0.0,       0.0,       tMu, 0.0 }, 
                        { 0.0, 0.0,       tLa,       0.0, 0.0 }, 
                        { 0.0, 0.0, tUz * tLa, tUy * tMu, 0.0 } };

                    // evaluate K_33
                    mK( 2 )( 2 ) = { 
                        { 0.0,       0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       tMu,       0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tMu,       0.0, 0.0 }, 
                        { 0.0,       0.0,       0.0,       tCh, 0.0 }, 
                        { 0.0, tUx * tMu, tUy * tMu, tUz * tCh, tKa } };

                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::assemble_A_matrices() - Number of space dimensions must be 2 or 3" );
                };
            }
        }                    

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
