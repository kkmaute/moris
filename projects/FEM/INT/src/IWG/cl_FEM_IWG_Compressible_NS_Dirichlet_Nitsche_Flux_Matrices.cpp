
//FEM/INT/src
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Dirichlet_Nitsche.hpp"
#include "fn_FEM_IWG_Compressible_NS.hpp"

//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Dirichlet_Nitsche::jump()
        {
            // check if the variable vectors have already been assembled
            if( !mJumpEval )
            {      
                return mJump;
            }         
            
            // set the eval flag
            mJumpEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType( 0 )  ), 
                    "IWG_Compressible_NS_Dirichlet_Nitsche::eval_jump() - check for residual DoF types failed." );

            //FIXME: only density and pressure primitive variable sets supported in this function

            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get the properties
            std::shared_ptr< Property > tPropPrescDof1 = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_1 ) );
            std::shared_ptr< Property > tPropPrescVel  = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_VELOCITY ) );
            std::shared_ptr< Property > tPropSelectVel = mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT_VELOCITY ) );
            std::shared_ptr< Property > tPropPrescDof3 = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_3 ) );

            // assemble field var vector and matrix based on number of spatial dimensions
            mJump.set_size( tNumSpaceDims + 2, 1, 0.0 );

            // set a default selection matrix if needed
            Matrix< DDRMat > tSelectVelMat;
            if ( tPropSelectVel == nullptr )
            {
                // set selection matrix as identity
                eye( tNumSpaceDims, tNumSpaceDims, tSelectVelMat );
            }
            else
            {
                tSelectVelMat = tPropSelectVel->val();
            }

            // for each variable check if something is prescribed 
            // first dof type
            if ( tPropPrescDof1 != nullptr )
            {
                // get field interpolator
                Field_Interpolator * tFIFirstDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

                // compute jump term for velocity, if prescribed
                mJump( 0 ) = tFIFirstDofType->val()( 0 ) - tPropPrescDof1->val()( 0 );
            }

            if ( tPropPrescVel != nullptr )
            {
                // compute jump term for velocity, if prescribed
                mJump( { 1, tNumSpaceDims } ) = tSelectVelMat * ( tFIVelocity->val() - tPropPrescVel->val() );
            }

            if ( tPropPrescDof3 != nullptr )
            {
                // get field interpolator
                Field_Interpolator * tFIThirdDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 2 ) );

                // compute jump term for third Dof, if prescribed
                mJump( tNumSpaceDims ) = tFIThirdDofType->val()( 0 ) - tPropPrescDof3->val()( 0 );
            }

            // return jump
            return mJump;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Dirichlet_Nitsche::dJumpdDOF()
        {
                        // check if the variable vectors have already been assembled
            if( !mJumpDofEval )
            {      
                return mdJumpdDOF;
            }         
            
            // set the eval flag
            mJumpDofEval = false;

            //FIXME: only density and pressure primitive variable sets supported in this function

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType( 0 )  ), 
                    "IWG_Compressible_NS_Dirichlet_Nitsche::mJumpDofEval() - check for residual DoF types failed." );

            // check Dof dependencies
            MORIS_ASSERT( check_dof_dependencies( mSet, mResidualDofType( 0 ), mRequestedMasterGlobalDofTypes ),
                    "IWG_Compressible_NS_Dirichlet_Nitsche::eval_dJumpdDOF() - List of Dof Dependencies not supported. See error messages above." );
            
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

            // get the properties
            std::shared_ptr< Property > tPropPrescDof1 = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_1 ) );
            std::shared_ptr< Property > tPropPrescVel  = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_VELOCITY ) );
            std::shared_ptr< Property > tPropSelectVel = mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT_VELOCITY ) );
            std::shared_ptr< Property > tPropPrescDof3 = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_3 ) );

            // assemble field var vector and matrix based on number of spatial dimensions
            mdJumpdDOF.set_size( tNumSpaceDims + 2, ( tNumSpaceDims + 2 ) * tNumBases, 0.0 );

            // set a default selection matrix if needed
            Matrix< DDRMat > tSelectVelMat;
            if ( tPropSelectVel == nullptr )
            {
                // set selection matrix as identity
                eye( tNumSpaceDims, tNumSpaceDims, tSelectVelMat );
            }
            else
            {
                tSelectVelMat = tPropSelectVel->val();
            }

            // for each variable check if something is prescribed, and add
            // first dof type
            if ( tPropPrescDof1 != nullptr )
            {
                // get field interpolator
                Field_Interpolator * tFIFirstVar =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

                // direct dependency of the state variable
                mdJumpdDOF( { 0, 0 }, { 0, tNumBases - 1 } ) = 
                        tFIFirstVar->N().matrix_data();
            }

            // second (velocity) DoF
            if ( tPropPrescVel != nullptr )
            {
                // direct dependency of the state variable
                mdJumpdDOF( { 1, tNumSpaceDims }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) = 
                        tSelectVelMat * tFIVelocity->N().matrix_data();
            }

            // third (temperature) DoF
            if ( tPropPrescDof3 != nullptr )
            {
                // get field interpolator
                Field_Interpolator * tFIThirdVar =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 2 ) );

                // direct dependency of the state variable
                mdJumpdDOF( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) = 
                        tFIThirdVar->N().matrix_data();
            }

            // loop over dof dependencies and add contributions from property derivatives
            for (uint iDof = 0; iDof < mRequestedMasterGlobalDofTypes.size(); iDof++)
            {   
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDof );

                // get index
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( iDof ), mtk::Master_Slave::MASTER );
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // first DoF type
                if ( tPropPrescDof1 != nullptr )
                {
                    if ( tPropPrescDof1->check_dof_dependency( tDofType ) )
                    {
                        // add contribution
                        mdJumpdDOF( { 0, 0 }, { tMasterDepStartIndex, tMasterDepStopIndex } ) -= 
                                tPropPrescDof1->dPropdDOF( tDofType ).matrix_data();
                    }
                }

                // second (velocity) DoF
                if ( tPropPrescVel != nullptr )
                {
                    if ( tPropPrescVel->check_dof_dependency( tDofType ) )
                    {
                        // direct dependency of the state variable
                        mdJumpdDOF( { 1, tNumSpaceDims }, { tMasterDepStartIndex, tMasterDepStopIndex } ) -= 
                                tPropPrescVel->dPropdDOF( tDofType ).matrix_data();
                    }
                }

                // third (temperature) DoF
                if ( tPropPrescDof3 != nullptr )
                {
                    if ( tPropPrescDof3->check_dof_dependency( tDofType ) )
                    {
                        // add contribution
                        mdJumpdDOF( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { tMasterDepStartIndex, tMasterDepStopIndex } ) -= 
                                tPropPrescDof3->dPropdDOF( tDofType ).matrix_data();
                    }
                }
            } // end: loop over dof dependencies

            // return matrix
            return mdJumpdDOF;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Dirichlet_Nitsche::A_Matrix( uint aAindex )
        {
            // check if the A matrices have already been assembled
            if( !mFluxAMatEval )
            {
                return mA( aAindex );
            }  
            
            // set the eval flag
            mFluxAMatEval = false;          

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // evaluate A-matrices and store them
            eval_A( tMM, tCM, mMasterFIManager, mResidualDofType( 0 ), mA );  

            // return 
            return mA( aAindex );
        }

       //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Dirichlet_Nitsche::eval_A_DOF_matrices()
        {
            // check if the A-Dof flux matrices have already been assembled
            if( !mFluxAMatDofEval )
            {
                return;
            }
            
            // set the eval flag
            mFluxAMatDofEval = false;

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
            eval_A0_DOF( tMM, tCM, mMasterFIManager, mResidualDofType( 0 ), mADOF( 0 ) ); 
            eval_A1_DOF( tMM, tCM, mMasterFIManager, mResidualDofType( 0 ), mADOF( 1 ) );  
            eval_A2_DOF( tMM, tCM, mMasterFIManager, mResidualDofType( 0 ), mADOF( 2 ) ); 
            if ( tNumSpaceDims == 3 )
            {
                eval_A3_DOF( tMM, tCM, mMasterFIManager, mResidualDofType( 0 ), mADOF( 3 ) );  
            } 
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Dirichlet_Nitsche::Traction()
        {
            // check if the A matrices have already been assembled
            if( !mTractionEval )
            {
                return mTraction;
            }  
            
            // set the eval flag
            mTractionEval = false;          

            // get the constitutive model
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get Kij*Y,j
            Matrix< DDRMat > tKijYj;
            eval_KijYj( tCM, mMasterFIManager, mResidualDofType( 0 ), tKijYj );

            // evaluate A-matrices and store them
            mTraction = trans( tKijYj ) * mNormal;

            // return 
            return mTraction;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Dirichlet_Nitsche::dTractiondDOF()
        {
            // check if the A matrices have already been assembled
            if( !mTractionDofEval )
            {
                return mTractionDOF;
            }  
            
            // set the eval flag
            mTractionDofEval = false; 

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType( 0 )  ), 
                    "IWG_Compressible_NS_Dirichlet_Nitsche::dTractiondDOF() - check for residual DoF types failed. See error message above." );

            // check Dof dependencies
            MORIS_ASSERT( check_dof_dependencies( mSet, mResidualDofType( 0 ), mRequestedMasterGlobalDofTypes ),
                    "IWG_Compressible_NS_Dirichlet_Nitsche::dTractiondDOF() - List of Dof Dependencies not supported. See error message above." );         

            // get the constitutive model
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

            // initialize
            mTractionDOF.set_size( tNumSpaceDims + 2, ( tNumSpaceDims + 2 ) * tNumBases, 0.0 );

            // loop over dof dependencies and add contributions from property derivatives
            for (uint iDof = 0; iDof < mRequestedMasterGlobalDofTypes.size(); iDof++)
            { 
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDof );

                // get dof indices for assembly
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( iDof ), mtk::Master_Slave::MASTER );
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // get d(Kij*Y,j)_dDOF
                moris::Cell< Matrix< DDRMat > > tKijYjDOF;
                eval_KijYjDOF( tCM, mMasterFIManager, mResidualDofType( 0 ), tDofType, tKijYjDOF );

                // loop over the different colums (=cells) of tKijYjDOF and multiply & assemble into dTractiondDOF
                for ( uint iCol = 0; iCol < tNumSpaceDims + 2; iCol++ )
                {   
                    // add contribution into matrix
                    mTractionDOF( { iCol, iCol }, { tMasterDepStartIndex, tMasterDepStopIndex } ) = trans( mNormal ) * tKijYjDOF( iCol );
                }
            }

            // return 
            return mTractionDOF;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Dirichlet_Nitsche::select_matrix()
        {
            // if select matrix has been evaluated, return it imidiately
            if ( !mSelectMatrixEval )
            {
                return mSelectMat;
            }

            // else, compute it
            mSelectMatrixEval = false;

            // get the properties
            std::shared_ptr< Property > tPropPrescDof1 = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_1 ) );
            std::shared_ptr< Property > tPropPrescVel  = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_VELOCITY ) );
            std::shared_ptr< Property > tPropSelectVel = mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT_VELOCITY ) );
            std::shared_ptr< Property > tPropPrescDof3 = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_3 ) );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get number of spatial dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // initialize Selection matrix
            mSelectMat.set_size( tNumSpaceDims + 2, tNumSpaceDims + 2, 0.0 );

            // check if 
            if ( tPropPrescDof1 != nullptr )
            {
                // set selection to active
                mSelectMat( 0, 0 ) = 1.0;
            }

            if ( tPropPrescVel != nullptr )
            {
                // if selection matrix prescribed by user, use it
                if ( tPropSelectVel != nullptr )
                {
                    // check selection matrix set by user
                    MORIS_ASSERT( ( tPropSelectVel->val().n_cols() == tNumSpaceDims ) and ( tPropSelectVel->val().n_rows() == tNumSpaceDims ),
                            "IWG_Compressible_NS_Dirichlet_Nitsche::compute_residual() - size of select matrix wrong." );
                    
                    // set selection matrix as identity
                    mSelectMat( { 1, tNumSpaceDims }, { 1, tNumSpaceDims } ) = tPropSelectVel->val().matrix_data();
                }
                else // else, use identity matrix
                {
                    // get identity matrix
                    Matrix< DDRMat > tIdentity;
                    eye( tNumSpaceDims, tNumSpaceDims, tIdentity);

                    // set select matrix to identity
                    mSelectMat( { 1, tNumSpaceDims }, { 1, tNumSpaceDims } ) = tIdentity.matrix_data();
                }
            } 

            if ( tPropPrescDof3 != nullptr )
            {
                // set selection to active
                mSelectMat( tNumSpaceDims + 1, tNumSpaceDims + 1 ) = 1.0;
            }

            // return select matrix
            return mSelectMat;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Dirichlet_Nitsche::test_functions()
        {
            // if select matrix has been evaluated, return it imidiately
            if ( !mTestFunctionsEval )
            {
                return mTestFunctions;
            }

            // else, compute it
            mTestFunctionsEval = false;

            // get the FIs associated with each residual dof type
            Field_Interpolator * tFIFirstDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
            Field_Interpolator * tFIVelocity     =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFIThirdDofType =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 2 ) );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();   

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();   

            // initialize Selection matrix
            mTestFunctions.set_size( ( tNumSpaceDims + 2 ) * tNumBases, tNumSpaceDims + 2, 0.0 );
 
            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index       = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof2Index       = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 1 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index       = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 2 ), mtk::Master_Slave::MASTER );
            uint tMasterRes1StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 0 );
            uint tMasterRes1StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 1 );
            uint tMasterRes2StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof2Index )( 0, 0 );
            uint tMasterRes2StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof2Index )( 0, 1 );
            uint tMasterRes3StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // assemble test function matrix
            mTestFunctions( { tMasterRes1StartIndex, tMasterRes1StopIndex }, { 0, 0 } ) = 
                    tFIFirstDofType->N_trans().matrix_data();
            mTestFunctions( { tMasterRes2StartIndex, tMasterRes2StopIndex }, { 1, tNumSpaceDims } ) = 
                    tFIVelocity->N_trans().matrix_data();
            mTestFunctions( { tMasterRes3StartIndex, tMasterRes3StopIndex }, { tNumSpaceDims + 1, tNumSpaceDims + 1 } ) = 
                    tFIThirdDofType->N_trans().matrix_data();

            // return matrix of test functions
            return mTestFunctions;
        }

        //------------------------------------------------------------------------------
        
    } /* namespace fem */
} /* namespace moris */
