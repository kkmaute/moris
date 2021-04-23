
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

        IWG_Compressible_NS_Dirichlet_Nitsche::IWG_Compressible_NS_Dirichlet_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "PrescribedDof1" ]     = static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_1 );
            mPropertyMap[ "PrescribedVelocity" ] = static_cast< uint >( IWG_Property_Type::PRESCRIBED_VELOCITY );
            mPropertyMap[ "SelectVelocity" ]     = static_cast< uint >( IWG_Property_Type::SELECT_VELOCITY );
            mPropertyMap[ "PrescribedDof3" ]     = static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_3 );
            mPropertyMap[ "Upwind" ]             = static_cast< uint >( IWG_Property_Type::UPWIND );

            // set size for the material model pointer cell
            mMasterMM.resize( static_cast< uint >( IWG_Material_Type::MAX_ENUM ), nullptr );

            // populate the material map
            mMaterialMap[ "FluidMM" ] = static_cast< uint >( IWG_Material_Type::FLUID_MM );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "FluidCM" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_CM );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitschePenaltyParameter" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Dirichlet_Nitsche::reset_spec_eval_flags()
        {
            // reset eval flags
            mSelectMatrixEval = true;
            mTestFunctionsEval = true;
            mJumpEval = true;
            mJumpDofEval = true;
            mFluxAMatEval = true;
            mFluxAMatDofEval = true;
            mTractionEval = true;
            mTractionDofEval = true;
            mTestTractionEval = true;
            mTestTractionDofEval = true;
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Dirichlet_Nitsche::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ), 
                    "IWG_Compressible_NS_Dirichlet_Nitsche::compute_residual() - Only pressure or density primitive variables supported for now."
                    " See error message above for specifics." );

            // get indeces for residual dof types, indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterRes1StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // get the Velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the upwind property
            std::shared_ptr< Property > tPropUpwind    = mMasterProp( static_cast< uint >( IWG_Property_Type::UPWIND ) );

            // get the Nitsche stabilization parameter - is a diagonal matrix, each diagonal entry corresponding to the respective Dof Type
            std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER ) );

            // get matrix subviews for different residuals - FIXME: assuming three different Residual DoF-Types, primitive vars
            auto tRes  = mSet->get_residual()( 0 )( { tMasterRes1StartIndex, tMasterRes3StopIndex }, { 0, 0 } );

            // get number of spatial dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // Boundary terms from Ibp
            tRes -= this->test_functions() * this->select_matrix() * this->Traction();

            // Upwind Term
            if ( tPropUpwind != nullptr )
            {
                // add A-matrices
                Matrix< DDRMat > tAini = this->A_Matrix( 1 ) * mNormal( 0 ) + this->A_Matrix( 2 ) * mNormal( 1 );
                if ( tNumSpaceDims == 3 )
                {
                    tAini = tAini + this->A_Matrix( 3 ) * mNormal( 2 );
                }

                // add contribution
                tRes -= tPropUpwind->val()( 0 ) * this->test_functions() * this->select_matrix() * tAini * this->jump();
            }

            // Nitsche Penalty Term
            if ( tSPNitsche != nullptr )
            {
                // convert vector of nitsche weights to diagonal matrix
                Matrix< DDRMat > tDiagSP( tNumSpaceDims + 2, tNumSpaceDims + 2, 0.0 );
                for ( uint iEntry = 0; iEntry < tNumSpaceDims + 2; iEntry++ )
                {
                    tDiagSP( iEntry, iEntry ) = tSPNitsche->val()( iEntry );
                }

                // add contribution
                tRes += this->test_functions() * tDiagSP * this->select_matrix() * this->jump();
            }

            // FIXME: adjoint term missing

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Dirichlet_Nitsche::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Dirichlet_Nitsche::compute_jacobian( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ), 
                    "IWG_Compressible_NS_Dirichlet_Nitsche::compute_jacobian() - Only pressure or density primitive variables supported for now." );

            // get residual dof indices for assembly (FIXME: assembly only for primitive vars)
            uint tMasterDof1Index      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDof3Index      = mSet->get_dof_index_for_type( mResidualDofType( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterRes1StartIndex = mSet->get_res_dof_assembly_map()( tMasterDof1Index )( 0, 0 );
            uint tMasterRes3StopIndex  = mSet->get_res_dof_assembly_map()( tMasterDof3Index )( 0, 1 );

            // check DoF dependencies
            MORIS_ASSERT( check_dof_dependencies( mSet, mResidualDofType, mRequestedMasterGlobalDofTypes ), 
                    "IWG_Compressible_NS_Dirichlet_Nitsche::compute_jacobian - Set of DoF dependencies not suppported. See error message above." );

            // get the indices for assembly - dependent dof types
            sint tDofFirstDepIndex     = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 0 )( 0 ), mtk::Master_Slave::MASTER );
            sint tDofThirdDepIndex     = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 2 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterDep1StartIndex = mSet->get_jac_dof_assembly_map()( tMasterDof1Index )( tDofFirstDepIndex, 0 );
            uint tMasterDep3StopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDof3Index )( tDofThirdDepIndex, 1 );  

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get number of space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();   

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();    

            // get the upwind property
            std::shared_ptr< Property > tPropUpwind = mMasterProp( static_cast< uint >( IWG_Property_Type::UPWIND ) );

            // get the Nitsche stabilization parameter - is a diagonal matrix, each diagonal entry corresponding to the respective Dof Type
            std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER ) );

            // get matrix subview for jacobian - FIXME: assuming three different Residual DoF-Types, primitive vars
            auto tJac  = mSet->get_jacobian()( { tMasterRes1StartIndex, tMasterRes3StopIndex }, { tMasterDep1StartIndex, tMasterDep3StopIndex } );

            // Boundary terms from Ibp
            tJac -= this->test_functions() * this->select_matrix() * this->dTractiondDOF();

            // Upwind Term
            if ( tPropUpwind != nullptr )
            {
                // add A-matrices
                Matrix< DDRMat > tAini = this->A_Matrix( 1 ) * mNormal( 0 ) + this->A_Matrix( 2 ) * mNormal( 1 );
                if ( tNumSpaceDims == 3 )
                {
                    tAini = tAini + this->A_Matrix( 3 ) * mNormal( 2 );
                }

                // add contribution - jump derivative
                tJac -= tPropUpwind->val()( 0 ) * this->test_functions() * this->select_matrix() * tAini * this->dJumpdDOF();

                // evaluate Dof derivs of A
                this->eval_A_DOF_matrices();

                // initialize
                moris::Cell< Matrix< DDRMat > > tSumADOF( tNumSpaceDims + 2 );
                Matrix< DDRMat > tADeltaY( tNumSpaceDims + 2, ( tNumSpaceDims + 2 ) * tNumBases );

                // add A-matrix derivatives and multiply with the jump term
                for ( uint iVar = 0; iVar < tNumSpaceDims + 2; iVar++ )
                {
                    // for each row of the A's, add them by multiplying with the normal
                    tSumADOF( iVar ) = mNormal( 0 ) * mADOF( 1 )( iVar ) + mNormal( 1 ) * mADOF( 2 )( iVar );
                    if ( tNumSpaceDims == 3 )
                    {
                        tSumADOF( iVar ) = mNormal( 2 ) * mADOF( 3 )( iVar );
                    }

                    // multiply with Y jump term 
                    tADeltaY( { iVar, iVar }, { tMasterDep1StartIndex, tMasterDep3StopIndex } ) = trans( this->jump() ) * tSumADOF( iVar );

                }

                // add contribution to jacobian
                tJac -= tPropUpwind->val()( 0 ) * this->test_functions() * this->select_matrix() * tADeltaY;

            }

            // Nitsche Penalty Term 
            if ( tSPNitsche != nullptr )
            {
                // convert vector of nitsche weights to diagonal matrix
                Matrix< DDRMat > tDiagSP( tNumSpaceDims + 2, tNumSpaceDims + 2, 0.0 );
                for ( uint iEntry = 0; iEntry < tNumSpaceDims + 2; iEntry++ )
                {
                    tDiagSP( iEntry, iEntry ) = tSPNitsche->val()( iEntry );
                }

                // add contribution
                tJac += this->test_functions() * tDiagSP * this->select_matrix() * this->dJumpdDOF();

                // FIXME: assuming no dependency of the penalty paramter on the dof types
            }

            // FIXME: adjoint term missing

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Dirichlet_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Dirichlet_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Dirichlet_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Dirichlet_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Dirichlet_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
        
    } /* namespace fem */
} /* namespace moris */
