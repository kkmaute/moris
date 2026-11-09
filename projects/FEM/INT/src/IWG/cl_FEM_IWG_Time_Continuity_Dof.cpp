/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Time_Continuity_Dof.cpp
 *
 */

#include "cl_FEM_IWG_Time_Continuity_Dof.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Enums.hpp"

// LINALG/src
#include "fn_trans.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_Time_Continuity_Dof::IWG_Time_Continuity_Dof()
    {
        // set IWG type
        mIWGType = moris::fem::IWG_Type::TIME_CONTINUITY_DOF;

        // set time continuity flag
        mTimeContinuity = true;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "WeightCurrent" ]    = static_cast< uint >( IWG_Property_Type::WEIGHT_CURRENT );
        mPropertyMap[ "WeightPrevious" ]   = static_cast< uint >( IWG_Property_Type::WEIGHT_PREVIOUS );
        mPropertyMap[ "InitialCondition" ] = static_cast< uint >( IWG_Property_Type::INITIAL_CONDITION );
        mPropertyMap[ "WeightResidual" ]   = static_cast< uint >( IWG_Property_Type::WEIGHT_RESIDUAL );
        mPropertyMap[ "Thickness" ]        = static_cast< uint >( IWG_Property_Type::THICKNESS );
        mPropertyMap[ "Lump" ]             = static_cast< uint >( IWG_Property_Type::LUMP );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Time_Continuity_Dof::compute_residual( real aWStar )
    {
        // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
        this->check_field_interpolators();
#endif
        // get leader index for residual dof type, indices for assembly
        uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get residual dof type field interpolator for current time step
        Field_Interpolator* tFICurrent =
                mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get residual dof type field interpolator for previous time step
        Field_Interpolator* tFIPrevious =
                mLeaderPreviousFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // store field manager of previous time step with field manager of current time step (previous state might be used in property)
        mLeaderFIManager->set_field_interpolator_manager_previous( mLeaderPreviousFIManager );

        // get current weight property
        const std::shared_ptr< Property >& tPropWeightCurrent =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT_CURRENT ) );

        // get previous weight property
        const std::shared_ptr< Property >& tPropWeightPrevious =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT_PREVIOUS ) );

        // get residual weight property
        const std::shared_ptr< Property >& tPropWeightResidual =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT_RESIDUAL ) );

        // get initial condition property
        const std::shared_ptr< Property >& tPropInitialCondition =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::INITIAL_CONDITION ) );

        // get thickness property
        const std::shared_ptr< Property >& tPropThickness =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // get lumping property
        const std::shared_ptr< Property >& tPropLump =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::LUMP ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // externally controlled weighting for residual (but not Jacobian) to enable
        // inexact Newton with time continuity contribution to Jacobian but not residual
        real tResWeight = 1.0;
        if ( tPropWeightResidual != nullptr )
        {
            tResWeight = tPropWeightResidual->val()( 0 );
        }

        // lumping parameter
        real tConsistent = 1.0;
        real tLump       = 0.0;

        if ( tPropLump != nullptr )
        {
            // compute weights for lumped and consistent terms
            tLump       = tPropLump->val()( 0 );
            tConsistent = 1.0 - tLump;
        }

        // FIXME set initial time
        real tInitTime = 0.0;

        // initialize jump in time
        Matrix< DDRMat > tJump = tPropWeightCurrent->val()( 0 ) * tFICurrent->val();

        // if not the first time step
        if ( mLeaderFIManager->get_IP_geometry_interpolator()->valt()( 0 ) > tInitTime )
        {
            // compute the jump
            tJump -= tPropWeightPrevious->val()( 0 ) * tFIPrevious->val();
        }
        // if first time step
        else
        {
            // compute the jump
            tJump -= tPropWeightPrevious->val()( 0 ) * tPropInitialCondition->val();
        }

        // add contribution to residual
        mSet->get_residual()( 0 )(
                { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                aWStar * tResWeight * tConsistent * ( tFICurrent->N_trans() * tJump );

        // consider lumping
        if ( tLump > 0.0 )
        {
            // generate lumped matrix
            Lumped_Matrix    tLumpCurrent;
            Lumped_Matrix    tLumpPrevious;
            Matrix< DDRMat > tInitialConditionCoeff;

            if ( mLeaderFIManager->get_IP_geometry_interpolator()->valt()( 0 ) > tInitTime )
            {
                compute_lumped_matrices(
                        tFICurrent,
                        tFIPrevious,
                        tLumpCurrent,
                        tLumpPrevious,
                        tInitialConditionCoeff,
                        Matrix< DDRMat >{},
                        false );

                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                        aWStar * tResWeight * tLump * (                                                                   //
                                tPropWeightCurrent->val()( 0 ) * tLumpCurrent * vectorize( tFICurrent->get_coeff() ) -    //
                                tPropWeightPrevious->val()( 0 ) * tLumpPrevious * vectorize( tFIPrevious->get_coeff() ) );
            }
            else
            {
                compute_lumped_matrices(
                        tFICurrent,
                        tFIPrevious,
                        tLumpCurrent,
                        tLumpPrevious,
                        tInitialConditionCoeff,
                        tPropInitialCondition->val(),
                        true );

                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                        aWStar * tResWeight * tLump * (                                                                   //
                                tPropWeightCurrent->val()( 0 ) * tLumpCurrent * vectorize( tFICurrent->get_coeff() ) -    //
                                tPropWeightPrevious->val()( 0 ) * tLumpPrevious * vectorize( tInitialConditionCoeff ) );
            }
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                "IWG_Time_Continuity_Dof::compute_residual - Residual contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Time_Continuity_Dof::compute_jacobian( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
        // check leader field interpolators
        this->check_field_interpolators();
#endif
        // get leader index for residual dof type, indices for assembly
        uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get residual dof type field interpolator for current time step
        Field_Interpolator* tFICurrent =
                mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get residual dof type field interpolator for previous time step
        Field_Interpolator* tFIPrevious =
                mLeaderPreviousFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // store field manager of previous time step with field manager of current time step (previous state might be used in property)
        mLeaderFIManager->set_field_interpolator_manager_previous( mLeaderPreviousFIManager );

        // get current weight property
        const std::shared_ptr< Property >& tPropWeightCurrent =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT_CURRENT ) );

        // get previous weight property
        const std::shared_ptr< Property >& tPropWeightPrevious =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT_PREVIOUS ) );

        // get initial condition property
        const std::shared_ptr< Property >& tPropInitialCondition =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::INITIAL_CONDITION ) );

        // get thickness property
        const std::shared_ptr< Property >& tPropThickness =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // get lumping property
        const std::shared_ptr< Property >& tPropLump =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::LUMP ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // lumping parameter
        real tConsistent = 1.0;
        real tLump       = 0.0;

        Lumped_Matrix tLumpCurrent;
        Lumped_Matrix tLumpPrevious;

        if ( tPropLump != nullptr )
        {
            // compute weights for lumped and consistent terms
            tLump       = tPropLump->val()( 0 );
            tConsistent = 1.0 - tLump;

            // generate lumped mass matrix
            if ( tLump > 0.0 )
            {
                Matrix< DDRMat > tInitialConditionCoeff;

                compute_lumped_matrices(
                        tFICurrent,
                        tFIPrevious,
                        tLumpCurrent,
                        tLumpPrevious,
                        tInitialConditionCoeff,
                        Matrix< DDRMat >{},
                        false );
            }
        }

        // get the number of leader dof type dependencies
        uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // FIXME set initial time
        real tInitTime = 0.0;

        // loop over leader dof type dependencies
        for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
        {
            // get the treated dof type
            const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

            // get the index for dof type, indices for assembly
            sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
            uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

            // if residual dof type
            if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                // add contribution to Jacobian
                mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                        aWStar * tConsistent * tPropWeightCurrent->val()( 0 ) * ( tFICurrent->N_trans() * tFICurrent->N() );

                // consider lumping
                if ( tLump > 0.0 )
                {
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                            aWStar * tLump * tPropWeightCurrent->val()( 0 ) * tLumpCurrent;
                }
            }

            // if current weight property has dependency on the dof type
            if ( tPropWeightCurrent->check_dof_dependency( tDofType ) )
            {
                // add contribution to Jacobian
                mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                        aWStar * tConsistent *    //
                        ( tFICurrent->N_trans() * tFICurrent->val() * tPropWeightCurrent->dPropdDOF( tDofType ) );

                // consider lumping
                if ( tLump > 0.0 )
                {
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                            aWStar * tLump *    //
                            ( tLumpCurrent * vectorize( tFICurrent->get_coeff() ) * tPropWeightCurrent->dPropdDOF( tDofType ) );
                }
            }

            // if previous weight property has dependency on the dof type
            if ( tPropWeightPrevious->check_dof_dependency( tDofType ) )
            {
                // if not the first time step
                if ( mLeaderFIManager->get_IP_geometry_interpolator()->valt()( 0 ) > tInitTime )
                {
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=
                            aWStar * tConsistent *    //
                            ( tFICurrent->N_trans() * tFIPrevious->val() * tPropWeightPrevious->dPropdDOF( tDofType ) );

                    // consider lumping
                    if ( tLump > 0.0 )
                    {
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=
                                aWStar * tLump *    //
                                ( tLumpPrevious * vectorize( tFIPrevious->get_coeff() ) * tPropWeightPrevious->dPropdDOF( tDofType ) );
                    }
                }
                else
                {
                    // Note: lumped and non lumped are combined
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=
                            aWStar * ( tFICurrent->N_trans() * tPropInitialCondition->val() * tPropWeightPrevious->dPropdDOF( tDofType ) );
                }
            }

            // if initial condition property has dependency on the dof type
            if ( tPropInitialCondition->check_dof_dependency( tDofType ) )
            {
                MORIS_ERROR( false,
                        " IWG_Time_Continuity_Dof::compute_jacobian - Dof dependence of initial conditions not implemented." );
            }
        }
        // FIXME add derivative for initial conditions?

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                "IWG_Time_Continuity_Dof::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Time_Continuity_Dof::compute_jacobian_and_residual( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Time_Continuity_Dof::compute_jacobian_and_residual - Not implemented." );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Time_Continuity_Dof::compute_dRdp( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Time_Continuity_Dof::compute_dRdp - Not implemented." );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Time_Continuity_Dof::compute_jacobian_previous( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
        // check leader field interpolators
        this->check_field_interpolators();
#endif

        // get leader index for residual dof type, indices for assembly
        uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get residual dof type field interpolator for current time step
        Field_Interpolator* tFICurrent =
                mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get residual dof type field interpolator for previous time step
        Field_Interpolator* tFIPrevious =
                mLeaderPreviousFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get previous weight property
        const std::shared_ptr< Property >& tPropWeightPrevious =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT_PREVIOUS ) );

        // get thickness property
        const std::shared_ptr< Property >& tPropThickness =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // get lumping property
        const std::shared_ptr< Property >& tPropLump =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::LUMP ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // lumping parameter
        real tConsistent = 1.0;
        real tLump       = 0.0;

        Lumped_Matrix tLumpPrevious;

        if ( tPropLump != nullptr )
        {
            // compute weights for lumped and consistent terms
            tLump       = tPropLump->val()( 0 );
            tConsistent = 1.0 - tLump;

            // generate lumped mass matrix
            if ( tLump > 0.0 )
            {
                compute_lumped_matrix( tFICurrent,
                        tFICurrent->get_number_of_time_bases() == 1 ? tFICurrent : tFIPrevious,
                        tLumpPrevious );
            }
        }

        // get the number of leader dof type dependencies
        uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // FIXME if not first time step
        //            MORIS_ASSERT( mLeaderFIManager->get_IP_geometry_interpolator()->valt()( 0 ) > 0.0,
        //                    "IWG_Time_Continuity_Dof::compute_jacobian_previous - should only be called after first time step" );

        // loop over leader dof type dependencies
        for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
        {
            // get the treated dof type
            const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

            // FIXME needs to be assemble on previous time step solution?
            // get the index for dof type, indices for assembly
            sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
            uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

            // if residual dof type
            if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                // add contribution to Jacobian
                mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=
                        aWStar * tConsistent * tPropWeightPrevious->val()( 0 ) * ( tFICurrent->N_trans() * tFIPrevious->N() );

                // consider lumping
                if ( tLump > 0.0 )
                {
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=
                            aWStar * tLump * tPropWeightPrevious->val()( 0 ) * tLumpPrevious;
                }
            }

            // note: there is no contribution from weight dof dependency because derivative of property
            // only with respect to current time step solution assumed
        }
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
