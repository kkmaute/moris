/*
 * cl_FEM_IWG_Compressible_NS_Bulk_Flux_Matrix_Derivs.cpp
 *
 *  Created on: Feb 23, 2021
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

        void IWG_Compressible_NS_Bulk::eval_A0_DOF()
        {
            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

            // get commonly used values
            real tRho    = tMM->density()( 0 );
            real tAlphaP = tMM->AlphaP()( 0 );
            real tBetaT  = tMM->BetaT()( 0 );
            real tCv     = tMM->Cv()( 0 );
            real tEtot   = tCM->Energy()( 0 );
            real tUx     = tFIVelocity->val()( 0 );
            real tUy     = tFIVelocity->val()( 1 );
            real tUz     = 0.0;
            if ( tNumSpaceDims == 3 )
            {
                tUz = tFIVelocity->val()( 2 );
            }      

            // divide the N-vector for the velocity
            Matrix< DDRMat > tNUx = tFIVelocity->N()( { 0, 0 }, { 0, tNumBases - 1 } );
            Matrix< DDRMat > tNUy = tFIVelocity->N()( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } );
            Matrix< DDRMat > tNUz;
            if ( tNumSpaceDims == 3 )
            {
                tNUz = tFIVelocity->N()( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } );
            } 

            // variable index for the third variable depending on spatial dimension
            uint tThirdVarIndex = tNumSpaceDims + 1;

            // =======================
            // Assemble A0 derivatives
            // =======================

            // derivative matrix for FIRST ROW OF A0
        
            mADOF( 0 )( 0 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } );

            mADOF( 0 )( 0 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } );                  

            mADOF( 0 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    -1.0 * tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) - tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } );

            mADOF( 0 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    -1.0 * tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) - tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } );                

            // derivative matrix for SECOND ROW OF A0
            
            mADOF( 0 )( 1 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUx * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 0 )( 1 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    tRho * tBetaT * tNUx;                    

            mADOF( 0 )( 1 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUx * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 0 )( 1 )( { 1, 1 }, { 0, tNumBases - 1 } ) = 
                    tMM->DensityDOF( { mResidualDofType( 0 ) } ).matrix_data();

            mADOF( 0 )( 1 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tMM->DensityDOF( { mResidualDofType( 2 ) } ).matrix_data();


            mADOF( 0 )( 1 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    tUx * ( -1.0 * tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) - tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 0 )( 1 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    -1.0 * tRho * tAlphaP * tNUx;  

            mADOF( 0 )( 1 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUx * ( -1.0 * tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) - tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) );


            // derivative matrix for THIRD ROW OF A0
            
            mADOF( 0 )( 2 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUy * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 0 )( 2 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    tRho * tBetaT * tNUy;                    

            mADOF( 0 )( 2 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUy * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 0 )( 2 )( { 2, 2 }, { 0, tNumBases - 1 } ) = 
                    tMM->DensityDOF( { mResidualDofType( 0 ) } ).matrix_data();

            mADOF( 0 )( 2 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tMM->DensityDOF( { mResidualDofType( 2 ) } ).matrix_data();


            mADOF( 0 )( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    tUy * ( -1.0 * tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) - tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 0 )( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    -1.0 * tRho * tAlphaP * tNUy;  

            mADOF( 0 )( 2 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUy * ( -1.0 * tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) - tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) );

            
            // derivative matrix for FOURTH ROW OF A0 IF 3D
            if ( tNumSpaceDims == 3 )
            {
                mADOF( 0 )( 3 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                        tUz * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

                mADOF( 0 )( 3 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                        tRho * tBetaT * tNUz;                    

                mADOF( 0 )( 3 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tUz * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


                mADOF( 0 )( 3 )( { 3, 3 }, { 0, tNumBases - 1 } ) = 
                        tMM->DensityDOF( { mResidualDofType( 0 ) } ).matrix_data();

                mADOF( 0 )( 3 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tMM->DensityDOF( { mResidualDofType( 2 ) } ).matrix_data();


                mADOF( 0 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                        tUz * ( -1.0 * tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) - tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

                mADOF( 0 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                        -1.0 * tRho * tAlphaP * tNUz;  

                mADOF( 0 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tUz * ( -1.0 * tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) - tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) );  
            }

            // derivative matrix for LAST ROW OF A0
            
            mADOF( 0 )( tThirdVarIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tEtot * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) + tBetaT * tCM->dEnergydDOF( { mResidualDofType( 0 ) } );

            mADOF( 0 )( tThirdVarIndex )( { 0, 0 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) = 
                    tBetaT * tCM->dEnergydDOF( { mResidualDofType( 1 ) } );                    

            mADOF( 0 )( tThirdVarIndex )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tEtot * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) + tBetaT * tCM->dEnergydDOF( { mResidualDofType( 2 ) } );


            mADOF( 0 )( tThirdVarIndex )( { 1, tNumSpaceDims }, { 0, tNumBases - 1 } ) = 
                    tFIVelocity->val() * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 0 )( tThirdVarIndex )( { 1, tNumSpaceDims }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) = 
                    tRho * tFIVelocity->N();

            mADOF( 0 )( tThirdVarIndex )( { 1, tNumSpaceDims }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tFIVelocity->val() * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            mADOF( 0 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    tCv * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->CpDOF( { mResidualDofType( 0 ) } ) - 
                    tEtot * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) - tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 0 ) } );

            mADOF( 0 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) = 
                    -1.0 * tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 1 ) } );                    

            mADOF( 0 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tCv * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->CpDOF( { mResidualDofType( 2 ) } ) - 
                    tEtot * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) - tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 2 ) } );                              
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::eval_A1_DOF()
        {
            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

            // get commonly used values
            real tRho    = tMM->density()( 0 );
            real tAlphaP = tMM->AlphaP()( 0 );
            real tBetaT  = tMM->BetaT()( 0 );
            real tCv     = tMM->Cv()( 0 );
            real tEtot   = tCM->Energy()( 0 );
            real tUx     = tFIVelocity->val()( 0 );
            real tUy     = tFIVelocity->val()( 1 );
            real tUz     = 0.0;
            if ( tNumSpaceDims == 3 )
            {
                tUz = tFIVelocity->val()( 2 );
            }      

            // divide the N-vector for the velocity
            Matrix< DDRMat > tNUx = tFIVelocity->N()( { 0, 0 }, { 0, tNumBases - 1 } );
            Matrix< DDRMat > tNUy = tFIVelocity->N()( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } );
            Matrix< DDRMat > tNUz;
            if ( tNumSpaceDims == 3 )
            {
                tNUz = tFIVelocity->N()( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } );
            } 

            // variable index for the third variable depending on spatial dimension
            uint tThirdVarIndex = tNumSpaceDims + 1;

            // =======================
            // Assemble A1 derivatives
            // =======================

            // derivative matrix for FIRST ROW OF A1
        
            mADOF( 1 )( 0 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUx * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 1 )( 0 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    tRho * tBetaT * tNUx;

            mADOF( 1 )( 0 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUx * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 1 )( 0 )( { 1, 1 }, { 0, tNumBases - 1 } ) = 
                    tMM->DensityDOF( { mResidualDofType( 0 ) } ).matrix_data(); 

            mADOF( 1 )( 0 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tMM->DensityDOF( { mResidualDofType( 2 ) } ).matrix_data();   


            mADOF( 1 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    -1.0 * tUx * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 1 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    -1.0 * tRho * tAlphaP * tNUx;

            mADOF( 1 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    -1.0 * tUx * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) );                


            // derivative matrix for SECOND ROW OF A1
            
            mADOF( 1 )( 1 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUx * tUx * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 1 )( 1 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    2.0 * tRho * tBetaT * tUx * tNUx;                    

            mADOF( 1 )( 1 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUx * tUx * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 1 )( 1 )( { 1, 1 }, { 0, tNumBases - 1 } ) = 
                    2.0 * tUx * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 1 )( 1 )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    2.0 * tRho * tNUx;

            mADOF( 1 )( 1 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    2.0 * tUx * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            mADOF( 1 )( 1 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    -1.0 * tUx * tUx * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 1 )( 1 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    -2.0 * tRho * tAlphaP * tUx * tNUx;  

            mADOF( 1 )( 1 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    -1.0 * tUx * tUx * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) );


            // derivative matrix for THIRD ROW OF A1
            
            mADOF( 1 )( 2 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUx * tUy * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 1 )( 2 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    tRho * tBetaT * tUy * tNUx;   

            mADOF( 1 )( 2 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    tRho * tBetaT * tUx * tNUy;                    

            mADOF( 1 )( 2 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUx * tUy * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 1 )( 2 )( { 1, 1 }, { 0, tNumBases - 1 } ) = 
                    tUy * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 1 )( 2 )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    tRho * tNUy;

            mADOF( 1 )( 2 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUy * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            mADOF( 1 )( 2 )( { 2, 2 }, { 0, tNumBases - 1 } ) = 
                    tUx * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 1 )( 2 )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    tRho * tNUx;

            mADOF( 1 )( 2 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUx * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            mADOF( 1 )( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    -1.0 * tUx * tUy * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 1 )( 2 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    -1.0 * tRho * tAlphaP * tUy * tNUx;

            mADOF( 1 )( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    -1.0 * tRho * tAlphaP * tUx * tNUy;  

            mADOF( 1 )( 2 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    -1.0 * tUx * tUy * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) );

            
            // derivative matrix for FOURTH ROW OF A1 IF 3D
            if ( tNumSpaceDims == 3 )
            {
                mADOF( 1 )( 3 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                        tUx * tUz * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

                mADOF( 1 )( 3 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                        tRho * tBetaT * tUz * tNUx;   

                mADOF( 1 )( 3 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                        tRho * tBetaT * tUx * tNUz;                    

                mADOF( 1 )( 3 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tUx * tUz * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


                mADOF( 1 )( 3 )( { 1, 1 }, { 0, tNumBases - 1 } ) = 
                        tUz * tMM->DensityDOF( { mResidualDofType( 0 ) } );

                mADOF( 1 )( 3 )( { 1, 1 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                        tRho * tNUz;

                mADOF( 1 )( 3 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tUz * tMM->DensityDOF( { mResidualDofType( 2 ) } );


                mADOF( 1 )( 3 )( { 3, 3 }, { 0, tNumBases - 1 } ) = 
                        tUx * tMM->DensityDOF( { mResidualDofType( 0 ) } );

                mADOF( 1 )( 3 )( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                        tRho * tNUx;

                mADOF( 1 )( 3 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tUx * tMM->DensityDOF( { mResidualDofType( 2 ) } );


                mADOF( 1 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                        -1.0 * tUx * tUz * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

                mADOF( 1 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) = 
                        -1.0 * tRho * tAlphaP * tUz * tNUx;

                mADOF( 1 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                        -1.0 * tRho * tAlphaP * tUx * tNUz;  

                mADOF( 1 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        -1.0 * tUx * tUz * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) ); 
            }

            // derivative matrix for LAST ROW OF A1
            
            mADOF( 1 )( tThirdVarIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUx * ( tEtot * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) + tBetaT * tCM->dEnergydDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 1 )( tThirdVarIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    ( tBetaT * tEtot + 1.0 ) * tNUx;
                    
            mADOF( 1 )( tThirdVarIndex )( { 0, 0 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) += 
                    tUx * tBetaT * tCM->dEnergydDOF( { mResidualDofType( 1 ) } );                    

            mADOF( 1 )( tThirdVarIndex )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUx * ( tEtot * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) + tBetaT * tCM->dEnergydDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 1 )( tThirdVarIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) = 
                   tCM->dEnergydDOF( { mResidualDofType( 0 ) } ) + tMM->PressureDOF( { mResidualDofType( 0 ) } ) + tUx * tUx * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 1 )( tThirdVarIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                   2.0 * tRho * tUx * tNUx;

            mADOF( 1 )( tThirdVarIndex )( { 1, 1 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) += 
                    tCM->dEnergydDOF( { mResidualDofType( 1 ) } ).matrix_data();

            mADOF( 1 )( tThirdVarIndex )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                   tCM->dEnergydDOF( { mResidualDofType( 2 ) } ) + tMM->PressureDOF( { mResidualDofType( 2 ) } ) + tUx * tUx * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            mADOF( 1 )( tThirdVarIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) = 
                   tUx * tUy * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 1 )( tThirdVarIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                   tRho * tUy * tNUx;

            mADOF( 1 )( tThirdVarIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                   tRho * tUx * tNUy;

            mADOF( 1 )( tThirdVarIndex )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                   tUx * tUy * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            if ( tNumSpaceDims == 3 ) // third velocity component only for 3D
            {
                mADOF( 1 )( tThirdVarIndex )( { 3, 3 }, { 0, tNumBases - 1 } ) = 
                    tUx * tUz * tMM->DensityDOF( { mResidualDofType( 0 ) } );

                mADOF( 1 )( tThirdVarIndex )( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    tRho * tUz * tNUx;

                mADOF( 1 )( tThirdVarIndex )( { 3, 3 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                    tRho * tUx * tNUz;

                mADOF( 1 )( tThirdVarIndex )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUx * tUz * tMM->DensityDOF( { mResidualDofType( 2 ) } );
            }

            mADOF( 1 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    tUx * ( tCv * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->CpDOF( { mResidualDofType( 0 ) } ) - 
                    tEtot * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) - tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 0 ) } ) ) ;

            mADOF( 1 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) = 
                    ( tRho * tCv - tAlphaP * tEtot ) * tNUx;  
            
            mADOF( 1 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) -= 
                    tUx * tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 1 ) } );                    

            mADOF( 1 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUx * ( tCv * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->CpDOF( { mResidualDofType( 2 ) } ) - 
                    tEtot * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) - tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 2 ) } ) );                              
        }  

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::eval_A2_DOF()
        {
            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

            // get commonly used values
            real tRho    = tMM->density()( 0 );
            real tAlphaP = tMM->AlphaP()( 0 );
            real tBetaT  = tMM->BetaT()( 0 );
            real tCv     = tMM->Cv()( 0 );
            real tEtot   = tCM->Energy()( 0 );
            real tUx     = tFIVelocity->val()( 0 );
            real tUy     = tFIVelocity->val()( 1 );
            real tUz     = 0.0;
            if ( tNumSpaceDims == 3 )
            {
                tUz = tFIVelocity->val()( 2 );
            }      

            // divide the N-vector for the velocity
            Matrix< DDRMat > tNUx = tFIVelocity->N()( { 0, 0 }, { 0, tNumBases - 1 } );
            Matrix< DDRMat > tNUy = tFIVelocity->N()( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } );
            Matrix< DDRMat > tNUz;
            if ( tNumSpaceDims == 3 )
            {
                tNUz = tFIVelocity->N()( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } );
            } 

            // variable index for the third variable depending on spatial dimension
            uint tThirdVarIndex = tNumSpaceDims + 1;

            // =======================
            // Assemble A2 derivatives
            // =======================

            // derivative matrix for FIRST ROW OF A2
        
            mADOF( 2 )( 0 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUy * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 2 )( 0 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    tRho * tBetaT * tNUy;

            mADOF( 2 )( 0 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUy * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 2 )( 0 )( { 2, 2 }, { 0, tNumBases - 1 } ) = 
                    tMM->DensityDOF( { mResidualDofType( 0 ) } ).matrix_data(); 

            mADOF( 2 )( 0 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tMM->DensityDOF( { mResidualDofType( 2 ) } ).matrix_data();   


            mADOF( 2 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    -1.0 * tUy * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 2 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    -1.0 * tRho * tAlphaP * tNUy;

            mADOF( 2 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    -1.0 * tUy * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) ); 

            // derivative matrix for SECOND ROW OF A2 (is equivalent ot third row of A1)
            
            mADOF( 2 )( 1 ) = mADOF( 1 )( 2 );

            // derivative matrix for THIRD ROW OF A2
            
            mADOF( 2 )( 2 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUy * tUy * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 2 )( 2 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    2.0 * tRho * tBetaT * tUy * tNUy;                    

            mADOF( 2 )( 2 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUy * tUy * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 2 )( 2 )( { 2, 2 }, { 0, tNumBases - 1 } ) = 
                    2.0 * tUy * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 2 )( 2 )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    2.0 * tRho * tNUy;

            mADOF( 2 )( 2 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    2.0 * tUy * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            mADOF( 2 )( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    -1.0 * tUy * tUy * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 2 )( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    -2.0 * tRho * tAlphaP * tUy * tNUy;  

            mADOF( 2 )( 2 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    -1.0 * tUy * tUy * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) );

            // derivative matrix for FOURTH ROW OF A2
            if ( tNumSpaceDims == 3 )
            {
                mADOF( 2 )( 3 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                        tUy * tUz * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

                mADOF( 2 )( 3 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                        tRho * tBetaT * tUz * tNUy;   

                mADOF( 2 )( 3 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                        tRho * tBetaT * tUy * tNUz;                    

                mADOF( 2 )( 3 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tUy * tUz * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


                mADOF( 2 )( 3 )( { 2, 2 }, { 0, tNumBases - 1 } ) = 
                        tUz * tMM->DensityDOF( { mResidualDofType( 0 ) } );

                mADOF( 2 )( 3 )( { 2, 2 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                        tRho * tNUz;

                mADOF( 2 )( 3 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tUz * tMM->DensityDOF( { mResidualDofType( 2 ) } );


                mADOF( 2 )( 3 )( { 3, 3 }, { 0, tNumBases - 1 } ) = 
                        tUy * tMM->DensityDOF( { mResidualDofType( 0 ) } );

                mADOF( 2 )( 3 )( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                        tRho * tNUy;

                mADOF( 2 )( 3 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tUy * tMM->DensityDOF( { mResidualDofType( 2 ) } );


                mADOF( 2 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                        -1.0 * tUy * tUz * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

                mADOF( 2 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                        -1.0 * tRho * tAlphaP * tUz * tNUy;

                mADOF( 2 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                        -1.0 * tRho * tAlphaP * tUy * tNUz;  

                mADOF( 2 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        -1.0 * tUy * tUz * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) ); 
            }

            // derivative matrix for LAST ROW OF A2
            
            mADOF( 2 )( tThirdVarIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUy * ( tEtot * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) + tBetaT * tCM->dEnergydDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 2 )( tThirdVarIndex )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    ( tBetaT * tEtot + 1.0 ) * tNUy;
                    
            mADOF( 2 )( tThirdVarIndex )( { 0, 0 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) += 
                    tUy * tBetaT * tCM->dEnergydDOF( { mResidualDofType( 1 ) } );                    

            mADOF( 2 )( tThirdVarIndex )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUy * ( tEtot * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) + tBetaT * tCM->dEnergydDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 2 )( tThirdVarIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) = 
                   tUx * tUy * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 2 )( tThirdVarIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                   tRho * tUy * tNUx;

            mADOF( 2 )( tThirdVarIndex )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                   tRho * tUx * tNUy;

            mADOF( 2 )( tThirdVarIndex )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                   tUx * tUy * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            mADOF( 2 )( tThirdVarIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) = 
                   tCM->dEnergydDOF( { mResidualDofType( 0 ) } ) + tMM->PressureDOF( { mResidualDofType( 0 ) } ) + tUy * tUy * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 2 )( tThirdVarIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                   2.0 * tRho * tUy * tNUy;

            mADOF( 2 )( tThirdVarIndex )( { 2, 2 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) += 
                    tCM->dEnergydDOF( { mResidualDofType( 1 ) } ).matrix_data();

            mADOF( 2 )( tThirdVarIndex )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                   tCM->dEnergydDOF( { mResidualDofType( 2 ) } ) + tMM->PressureDOF( { mResidualDofType( 2 ) } ) + tUy * tUy * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            if ( tNumSpaceDims == 3 ) // third velocity component only for 3D
            {
                mADOF( 2 )( tThirdVarIndex )( { 3, 3 }, { 0, tNumBases - 1 } ) = 
                    tUy * tUz * tMM->DensityDOF( { mResidualDofType( 0 ) } );

                mADOF( 2 )( tThirdVarIndex )( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    tRho * tUz * tNUy;

                mADOF( 2 )( tThirdVarIndex )( { 3, 3 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                    tRho * tUy * tNUz;

                mADOF( 2 )( tThirdVarIndex )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUy * tUz * tMM->DensityDOF( { mResidualDofType( 2 ) } );
            }

            mADOF( 2 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    tUy * ( tCv * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->CpDOF( { mResidualDofType( 0 ) } ) - 
                    tEtot * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) - tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 0 ) } ) ) ;

            mADOF( 2 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                    ( tRho * tCv - tAlphaP * tEtot ) * tNUy;  
            
            mADOF( 2 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) -= 
                    tUy * tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 1 ) } );                    

            mADOF( 2 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUy * ( tCv * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->CpDOF( { mResidualDofType( 2 ) } ) - 
                    tEtot * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) - tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 2 ) } ) );   

        }   

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Bulk::eval_A3_DOF()
        {
            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

            // get commonly used values
            real tRho    = tMM->density()( 0 );
            real tAlphaP = tMM->AlphaP()( 0 );
            real tBetaT  = tMM->BetaT()( 0 );
            real tCv     = tMM->Cv()( 0 );
            real tEtot   = tCM->Energy()( 0 );
            real tUx     = tFIVelocity->val()( 0 );
            real tUy     = tFIVelocity->val()( 1 );
            real tUz     = tFIVelocity->val()( 2 );   

            // divide the N-vector for the velocity
            Matrix< DDRMat > tNUx = tFIVelocity->N()( { 0, 0 }, { 0, tNumBases - 1 } );
            Matrix< DDRMat > tNUy = tFIVelocity->N()( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } );
            Matrix< DDRMat > tNUz = tFIVelocity->N()( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } );

            // variable index for the third variable depending on spatial dimension
            uint tThirdVarIndex = tNumSpaceDims + 1;

            // =======================
            // Assemble A3 derivatives
            // =======================

            // derivative matrix for FIRST ROW OF A3
        
            mADOF( 3 )( 0 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUz * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 3 )( 0 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                    tRho * tBetaT * tNUz;

            mADOF( 3 )( 0 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUz * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 3 )( 0 )( { 3, 3 }, { 0, tNumBases - 1 } ) = 
                    tMM->DensityDOF( { mResidualDofType( 0 ) } ).matrix_data(); 

            mADOF( 3 )( 0 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tMM->DensityDOF( { mResidualDofType( 2 ) } ).matrix_data();   


            mADOF( 3 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    -1.0 * tUz * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 3 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                    -1.0 * tRho * tAlphaP * tNUz;

            mADOF( 3 )( 0 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    -1.0 * tUz * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) ); 


            // derivative matrix for SECOND ROW OF A3 (identical to fourth row of A1)

            mADOF( 3 )( 1 ) = mADOF( 1 )( 3 );


            // derivative matrix for THIRD ROW OF A3 (identical to fourth row of A2)

            mADOF( 3 )( 2 ) = mADOF( 2 )( 3 );


            // derivative matrix for FOURTH ROW OF A3

            mADOF( 3 )( 3 )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUz * tUz * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 3 )( 3 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                    2.0 * tRho * tBetaT * tUz * tNUz;                    

            mADOF( 3 )( 3 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUz * tUz * ( tBetaT * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 3 )( 3 )( { 3, 3 }, { 0, tNumBases - 1 } ) = 
                    2.0 * tUz * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 3 )( 3 )( { 3, 3 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                    2.0 * tRho * tNUz;

            mADOF( 3 )( 3 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    2.0 * tUz * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            mADOF( 3 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    -1.0 * tUz * tUz * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 3 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                    -2.0 * tRho * tAlphaP * tUz * tNUz;  

            mADOF( 3 )( 3 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    -1.0 * tUz * tUz * ( tAlphaP * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) );


            // derivative matrix for LAST ROW OF A3

            mADOF( 3 )( tThirdVarIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                    tUz * ( tEtot * tMM->BetaTDOF( { mResidualDofType( 0 ) } ) + tBetaT * tCM->dEnergydDOF( { mResidualDofType( 0 ) } ) );

            mADOF( 3 )( tThirdVarIndex )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                    ( tBetaT * tEtot + 1.0 ) * tNUz;
                    
            mADOF( 3 )( tThirdVarIndex )( { 0, 0 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) += 
                    tUz * tBetaT * tCM->dEnergydDOF( { mResidualDofType( 1 ) } );                    

            mADOF( 3 )( tThirdVarIndex )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUz * ( tEtot * tMM->BetaTDOF( { mResidualDofType( 2 ) } ) + tBetaT * tCM->dEnergydDOF( { mResidualDofType( 2 ) } ) );


            mADOF( 3 )( tThirdVarIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) = 
                   tUx * tUz * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 3 )( tThirdVarIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = 
                   tRho * tUz * tNUx;

            mADOF( 3 )( tThirdVarIndex )( { 1, 1 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                   tRho * tUx * tNUz;

            mADOF( 3 )( tThirdVarIndex )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                   tUx * tUz * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            mADOF( 3 )( tThirdVarIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) = 
                tUy * tUz * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 3 )( tThirdVarIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 
                tRho * tUz * tNUy;

            mADOF( 3 )( tThirdVarIndex )( { 2, 2 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                tRho * tUy * tNUz;

            mADOF( 3 )( tThirdVarIndex )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                tUy * tUz * tMM->DensityDOF( { mResidualDofType( 2 ) } );


            mADOF( 3 )( tThirdVarIndex )( { 3, 3 }, { 0, tNumBases - 1 } ) = 
                   tCM->dEnergydDOF( { mResidualDofType( 0 ) } ) + tMM->PressureDOF( { mResidualDofType( 0 ) } ) + tUz * tUz * tMM->DensityDOF( { mResidualDofType( 0 ) } );

            mADOF( 3 )( tThirdVarIndex )( { 3, 3 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                   2.0 * tRho * tUz * tNUz;

            mADOF( 3 )( tThirdVarIndex )( { 3, 3 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) += 
                    tCM->dEnergydDOF( { mResidualDofType( 1 ) } ).matrix_data();

            mADOF( 3 )( tThirdVarIndex )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                   tCM->dEnergydDOF( { mResidualDofType( 2 ) } ) + tMM->PressureDOF( { mResidualDofType( 2 ) } ) + tUz * tUz * tMM->DensityDOF( { mResidualDofType( 2 ) } );



            mADOF( 3 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) = 
                    tUz * ( tCv * tMM->DensityDOF( { mResidualDofType( 0 ) } ) + tRho * tMM->CpDOF( { mResidualDofType( 0 ) } ) - 
                    tEtot * tMM->AlphaPDOF( { mResidualDofType( 0 ) } ) - tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 0 ) } ) ) ;

            mADOF( 3 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) = 
                    ( tRho * tCv - tAlphaP * tEtot ) * tNUz;  
            
            mADOF( 3 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) -= 
                    tUz * tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 1 ) } );                    

            mADOF( 3 )( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                    tUz * ( tCv * tMM->DensityDOF( { mResidualDofType( 2 ) } ) + tRho * tMM->CpDOF( { mResidualDofType( 2 ) } ) - 
                    tEtot * tMM->AlphaPDOF( { mResidualDofType( 2 ) } ) - tAlphaP * tCM->dEnergydDOF( { mResidualDofType( 2 ) } ) );  

        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
