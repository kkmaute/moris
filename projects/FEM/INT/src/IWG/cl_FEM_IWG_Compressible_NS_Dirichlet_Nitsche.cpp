
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
            mPropertyMap[ "PrescribedDof1" ]      = static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_1 );
            mPropertyMap[ "PrescribedVelocity" ]  = static_cast< uint >( IWG_Property_Type::PRESCRIBED_VELOCITY );
            mPropertyMap[ "SelectVelocity" ]      = static_cast< uint >( IWG_Property_Type::SELECT_VELOCITY );
            mPropertyMap[ "PrescribedDof3" ]      = static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_3 );
            mPropertyMap[ "PressureUpwind" ]      = static_cast< uint >( IWG_Property_Type::PRESSUREUPWIND );
            mPropertyMap[ "DynamicViscosity" ]    = static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY );
            mPropertyMap[ "ThermalConductivity" ] = static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY );

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

        void IWG_Compressible_NS_Dirichlet_Nitsche::reset_child_eval_flags()
        {
            // reset eval flags
            mJumpEval = true;
            mJumpDofEval = true;
            mSelectMatrixEval = true;

            mTractionEval = true;
            mTractionDofEval = true;
            mTestTractionEval = true;
            mTestTractionDofEval = true;

            mUpwindOperatorEval = true;
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Dirichlet_Nitsche::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif
            // check residual dof types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ), 
                    "IWG_Compressible_NS_Dirichlet_Nitsche::compute_residual() - Only pressure or density primitive variables supported for now."
                    " See error message above for specifics." );

            // get number of space dimensions
            uint tNumSpaceDims = this->num_space_dims();

            // get total number of DoFs on Comp Flow Element
            uint tNumTotalBases = ( tNumSpaceDims + 2 ) * this->num_bases();

            // construct temporary Vector for residual
            Matrix< DDRMat > tTempRes( tNumTotalBases, 1, 0.0 );
            auto tRes = tTempRes( { 0, tNumTotalBases - 1 }, { 0, 0 } );

            // Boundary terms from Ibp
            tRes -= this->W_trans() * this->select_matrix() * this->Traction();

            // adjoint term
            tRes -= mBeta * this->TestTraction() * this->jump();

            // get the Nitsche stabilization parameter - is a diagonal matrix, each diagonal entry corresponding to the respective Dof Type
            std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER ) );

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
                tRes += this->W_trans() * tDiagSP * this->jump();
            }

            // get the upwind property
            std::shared_ptr< Property > tPropUpwind = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESSUREUPWIND ) );

            // Upwind Term
            if ( tPropUpwind != nullptr )
            {
                // add contribution using the upwind operator
                tRes -= tPropUpwind->val()( 0 ) * this->W_trans() * this->UpwindOperator() * this->jump();
            }

            // assemble into set residual
            this->assemble_residual( tTempRes );

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

            // check DoF dependencies
            MORIS_ASSERT( check_dof_dependencies( mSet, mResidualDofType, mRequestedMasterGlobalDofTypes ), 
                    "IWG_Compressible_NS_Dirichlet_Nitsche::compute_jacobian() - Set of DoF dependencies not suppported. See error message above." );

            // get number of space dimensions
            uint tNumSpaceDims = this->num_space_dims();

            // get total number of DoFs on Comp Flow Element
            uint tNumTotalBases = ( tNumSpaceDims + 2 ) * this->num_bases();

            // construct temporary Matrix for elemental Jacobian in standardized format
            Matrix< DDRMat > tTempJac( tNumTotalBases, tNumTotalBases, 0.0 );
            auto tJac = tTempJac( { 0, tNumTotalBases - 1 }, { 0, tNumTotalBases - 1 } );

            // Boundary terms from Ibp
            tJac -= this->W_trans() * this->select_matrix() * this->dTractiondDOF();

            // adjoint term
            tJac -= mBeta * this->TestTraction() * this->dJumpdDOF();
            tJac -= mBeta * this->dTestTractiondDOF( this->jump() );

            // get the Nitsche stabilization parameter - is a diagonal matrix, each diagonal entry corresponding to the respective Dof Type
            std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER ) );

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
                tJac += this->W_trans() * tDiagSP * this->dJumpdDOF();

                // FIXME: assuming no dependency of the penalty paramter on the dof types
            }

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the upwind property
            std::shared_ptr< Property > tPropUpwind = mMasterProp( static_cast< uint >( IWG_Property_Type::PRESSUREUPWIND ) );

            // Upwind Term
            if ( tPropUpwind != nullptr )
            {
                // add contribution
                tJac -= tPropUpwind->val()( 0 ) * this->W_trans() * ( 
                        this->UpwindOperator() * this->dJumpdDOF() + this->dUpwindOperatordY( this->jump() ) * this->W() );
            }

            // assemble jacobian into set jacobian
            this->assemble_jacobian( tTempJac );

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
