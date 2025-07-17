/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Dirichlet_Nitsche_Flux_Matrices.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Dirichlet_Nitsche.hpp"
#include "fn_FEM_IWG_Compressible_NS.hpp"

//LINALG/src
#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &IWG_Compressible_NS_Dirichlet_Nitsche::jump()
    {
        // check if the variable vectors have already been assembled
        if ( !mJumpEval )
        {
            return mJump;
        }

        // set the eval flag
        mJumpEval = false;

        // check residual DoF types
        MORIS_ASSERT( check_residual_dof_types( mResidualDofType ),
                "IWG_Compressible_NS_Dirichlet_Nitsche::eval_jump() - check for residual DoF types failed." );

        // get number of space dimensions
        uint tNumSpaceDims = this->num_space_dims();

        // get the properties
        std::shared_ptr< Property > tPropPrescDof1 = mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_1 ) );
        std::shared_ptr< Property > tPropPrescVel  = mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_VELOCITY ) );
        std::shared_ptr< Property > tPropSelectVel = mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT_VELOCITY ) );
        std::shared_ptr< Property > tPropPrescDof3 = mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_3 ) );

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
            Field_Interpolator *tFIFirstDofType = mLeaderFIManager->get_field_interpolators_for_type( this->get_primary_state_var( 0 ) );

            // compute jump term for velocity, if prescribed
            mJump( 0 ) = tFIFirstDofType->val()( 0 ) - tPropPrescDof1->val()( 0 );
        }

        if ( tPropPrescVel != nullptr )
        {
            // get field interpolator
            Field_Interpolator *tFIVelocity = mLeaderFIManager->get_field_interpolators_for_type( this->get_primary_state_var( 1 ) );

            // compute jump term for velocity, if prescribed
            mJump( { 1, tNumSpaceDims } ) = tSelectVelMat * ( tFIVelocity->val() - tPropPrescVel->val() );
        }

        if ( tPropPrescDof3 != nullptr )
        {
            // get field interpolator
            Field_Interpolator *tFIThirdDofType = mLeaderFIManager->get_field_interpolators_for_type( this->get_primary_state_var( 2 ) );

            // compute jump term for third Dof, if prescribed
            mJump( tNumSpaceDims + 1 ) = tFIThirdDofType->val()( 0 ) - tPropPrescDof3->val()( 0 );
        }

        // return jump
        return mJump;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &IWG_Compressible_NS_Dirichlet_Nitsche::dJumpdDOF()
    {
        // check if the variable vectors have already been assembled
        if ( !mJumpDofEval )
        {
            return mdJumpdDOF;
        }

        // set the eval flag
        mJumpDofEval = false;

        // FIXME: only density and pressure primitive variable sets supported in this function

        // check residual DoF types
        MORIS_ASSERT( check_residual_dof_types( mResidualDofType ),
                "IWG_Compressible_NS_Dirichlet_Nitsche::mJumpDofEval() - check for residual DoF types failed." );

        // check Dof dependencies
        MORIS_ASSERT( check_dof_dependencies( mSet, mResidualDofType, mRequestedLeaderGlobalDofTypes ),
                "IWG_Compressible_NS_Dirichlet_Nitsche::eval_dJumpdDOF() - List of Dof Dependencies not supported. See error messages above." );

        // get number of space dimensions
        uint tNumSpaceDims = this->num_space_dims();

        // get number of bases for the elements used
        uint tNumBases = this->num_bases();

        // get the properties
        std::shared_ptr< Property > tPropPrescDof1 = mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_1 ) );
        std::shared_ptr< Property > tPropPrescVel  = mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_VELOCITY ) );
        std::shared_ptr< Property > tPropSelectVel = mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT_VELOCITY ) );
        std::shared_ptr< Property > tPropPrescDof3 = mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_3 ) );

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
            Field_Interpolator *tFIFirstVar = mLeaderFIManager->get_field_interpolators_for_type( this->get_primary_state_var( 0 ) );

            // direct dependency of the state variable
            mdJumpdDOF( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tFIFirstVar->N().matrix_data();
        }

        // second (velocity) DoF
        if ( tPropPrescVel != nullptr )
        {
            // get field interpolator
            Field_Interpolator *tFIVelocity = mLeaderFIManager->get_field_interpolators_for_type( this->get_primary_state_var( 1 ) );

            // direct dependency of the state variable
            mdJumpdDOF( { 1, tNumSpaceDims }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) =
                    tSelectVelMat * tFIVelocity->N().matrix_data();
        }

        // third (temperature) DoF
        if ( tPropPrescDof3 != nullptr )
        {
            // get field interpolator
            Field_Interpolator *tFIThirdVar = mLeaderFIManager->get_field_interpolators_for_type( this->get_primary_state_var( 2 ) );

            // direct dependency of the state variable
            mdJumpdDOF( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tFIThirdVar->N().matrix_data();
        }

        // loop over dof dependencies and add contributions from property derivatives
        for ( uint iDof = 0; iDof < mRequestedLeaderGlobalDofTypes.size(); iDof++ )
        {
            // get the treated dof type
            Vector< MSI::Dof_Type > &tDofType = mRequestedLeaderGlobalDofTypes( iDof );

            // get index
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( iDof )( 0 ), mtk::Leader_Follower::LEADER );
            sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
            uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

            // first DoF type
            if ( tPropPrescDof1 != nullptr )
            {
                if ( tPropPrescDof1->check_dof_dependency( tDofType ) )
                {
                    // add contribution
                    mdJumpdDOF( { 0, 0 }, { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=
                            tPropPrescDof1->dPropdDOF( tDofType ).matrix_data();
                }
            }

            // second (velocity) DoF
            if ( tPropPrescVel != nullptr )
            {
                if ( tPropPrescVel->check_dof_dependency( tDofType ) )
                {
                    // direct dependency of the state variable
                    mdJumpdDOF( { 1, tNumSpaceDims }, { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=
                            tPropPrescVel->dPropdDOF( tDofType ).matrix_data();
                }
            }

            // third (temperature) DoF
            if ( tPropPrescDof3 != nullptr )
            {
                if ( tPropPrescDof3->check_dof_dependency( tDofType ) )
                {
                    // add contribution
                    mdJumpdDOF( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=
                            tPropPrescDof3->dPropdDOF( tDofType ).matrix_data();
                }
            }
        }    // end: loop over dof dependencies

        // return matrix
        return mdJumpdDOF;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &IWG_Compressible_NS_Dirichlet_Nitsche::select_matrix()
    {
        // if select matrix has been evaluated, return it imidiately
        if ( !mSelectMatrixEval )
        {
            return mSelectMat;
        }

        // else, compute it
        mSelectMatrixEval = false;

        // get the properties
        std::shared_ptr< Property > tPropPrescDof1 = mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_1 ) );
        std::shared_ptr< Property > tPropPrescVel  = mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_VELOCITY ) );
        std::shared_ptr< Property > tPropSelectVel = mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT_VELOCITY ) );
        std::shared_ptr< Property > tPropPrescDof3 = mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_DOF_3 ) );

        // get number of spatial dimensions
        uint tNumSpaceDims = this->num_space_dims();

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
                        "IWG_Compressible_NS_Dirichlet_Nitsche::select_matrix() - size of select matrix wrong." );

                // set selection matrix as identity
                mSelectMat( { 1, tNumSpaceDims }, { 1, tNumSpaceDims } ) = tPropSelectVel->val().matrix_data();
            }
            else    // else, use identity matrix
            {
                // get identity matrix
                Matrix< DDRMat > tIdentity;
                eye( tNumSpaceDims, tNumSpaceDims, tIdentity );

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
    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &IWG_Compressible_NS_Dirichlet_Nitsche::Traction()
    {
        // check if the A matrices have already been assembled
        if ( !mTractionEval )
        {
            return mTraction;
        }

        // set the eval flag
        mTractionEval = false;

        // number of state variables and total bases
        uint tNumStateVars = this->num_space_dims() + 2;

        // initialize traction
        mTraction.set_size( tNumStateVars, 1, 0.0 );
        auto tTraction = mTraction( { 0, tNumStateVars - 1 }, { 0, 0 } );

        // loop over K-matrices
        for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
        {
            for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
            {
                // add contribution
                tTraction += this->K( iDim, jDim ) * this->dYdx( jDim ) * mNormal( iDim );
            }
        }

        // return
        return mTraction;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &IWG_Compressible_NS_Dirichlet_Nitsche::dTractiondDOF()
    {
        // check if the A matrices have already been assembled
        if ( !mTractionDofEval )
        {
            return mTractionDOF;
        }

        // set the eval flag
        mTractionDofEval = false;

        // get the properties
        std::shared_ptr< Property > tPropMu    = mLeaderProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
        std::shared_ptr< Property > tPropKappa = mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

        // number of state variables and total bases
        uint tNumStateVars = this->num_space_dims() + 2;
        uint tTotNumBases  = tNumStateVars * this->num_bases();

        // initialize traction
        mTractionDOF.set_size( tNumStateVars, tTotNumBases, 0.0 );
        auto tTractionDOF = mTractionDOF( { 0, tNumStateVars - 1 }, { 0, tTotNumBases - 1 } );

        // initialize storage matrix for variable derivative of K
        Matrix< DDRMat > dKijdY_dYdxj;

        // loop over K-matrices
        for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
        {
            for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
            {
                // get dKij/dY * Y,j
                eval_dKdY_VR( tPropMu, tPropKappa, mLeaderFIManager, this->dYdx( jDim ), iDim, jDim, dKijdY_dYdxj );

                // add contribution
                tTractionDOF += mNormal( iDim ) * ( dKijdY_dYdxj * this->W() + this->K( iDim, jDim ) * this->dWdx( jDim ) );
            }
        }

        // return
        return mTractionDOF;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &IWG_Compressible_NS_Dirichlet_Nitsche::TestTraction()
    {
        // check if matrix has already been assembled
        if ( !mTestTractionEval )
        {
            return mTestTraction;
        }

        // set the eval flag
        mTestTractionEval = false;

        // get test traction
        mTestTraction = trans( this->dTractiondDOF() );

        // return value
        return mTestTraction;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &IWG_Compressible_NS_Dirichlet_Nitsche::dTestTractiondDOF( const Matrix< DDRMat > &aVL )
    {
        // check if matrix has already been assembled
        if ( !mTestTractionDofEval )
        {
            return mTestTractionDOF;
        }

        // set the eval flag
        mTestTractionDofEval = false;

        // get number of basis functions for complete variable set
        uint tTotNumBases = ( this->num_space_dims() + 2 ) * this->num_bases();

        // initialize test traction dof deriv
        mTestTractionDOF.set_size( tTotNumBases, tTotNumBases, 0.0 );
        auto tTestTractionDOF = mTestTractionDOF( { 0, tTotNumBases - 1 }, { 0, tTotNumBases - 1 } );

        // get the properties
        std::shared_ptr< Property > tPropMu    = mLeaderProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
        std::shared_ptr< Property > tPropKappa = mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

        // initialize storage matrix for variable derivative of K
        Matrix< DDRMat > tdKdY_jump;

        // loop over K-matrices
        for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
        {
            for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
            {
                // get dKij/dY * Y,j
                eval_VL_dKdY( tPropMu, tPropKappa, mLeaderFIManager, aVL, iDim, jDim, tdKdY_jump );

                // add contribution
                Matrix< DDRMat > tdTijdDOF = mNormal( iDim ) * trans( this->dWdx( jDim ) ) * tdKdY_jump * this->W();
                tTestTractionDOF += tdTijdDOF + trans( tdTijdDOF );
            }
        }

        // return value
        return mTestTractionDOF;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &IWG_Compressible_NS_Dirichlet_Nitsche::UpwindOperator()
    {
        // check if matrix has already been assembled
        if ( !mUpwindOperatorEval )
        {
            return mUpwindOperator;
        }

        // set the eval flag
        mTestTractionEval = false;

        // get number of space dimensions
        uint tNumSpaceDims = this->num_space_dims();

        // initialize the operator
        mUpwindOperator.set_size( tNumSpaceDims + 2, tNumSpaceDims + 2, 0.0 );

        // place normal such that it gets multiplied with the pressure difference in the velocity residual
        mUpwindOperator( { 1, tNumSpaceDims }, { 0, 0 } ) = mNormal.matrix_data();

        // for temperature residual place normal dotted against velocity such that it gets multiplied with the pressure difference
        Matrix< DDRMat > tVelVec                = this->Y()( { 1, tNumSpaceDims } );
        mUpwindOperator( tNumSpaceDims + 1, 0 ) = dot( mNormal, tVelVec );

        // return value
        return mUpwindOperator;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &IWG_Compressible_NS_Dirichlet_Nitsche::dUpwindOperatordY( const Matrix< DDRMat > &aVR )
    {
        // get number of space dimensions
        uint tNumSpaceDims = this->num_space_dims();

        // initialize the operator
        mdUpwindOperatordY.set_size( tNumSpaceDims + 2, tNumSpaceDims + 2, 0.0 );

        // place normal such that it gets multiplied with the pressure difference in the velocity residual
        mdUpwindOperatordY( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 1, tNumSpaceDims } ) = aVR( 0 ) * trans( mNormal );

        // return value
        return mdUpwindOperatordY;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
