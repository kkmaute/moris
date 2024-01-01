/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Incompressible_NS_Pressure_Neumann.cpp
 *
 */

#include "cl_FEM_IWG_Incompressible_NS_Pressure_Neumann.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Incompressible_NS_Pressure_Neumann::IWG_Incompressible_NS_Pressure_Neumann()
        {
            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Pressure" ]           = static_cast< uint >( IWG_Property_Type::PRESSURE );
            mPropertyMap[ "TotalPressure" ]      = static_cast< uint >( IWG_Property_Type::TOTAL_PRESSURE );
            mPropertyMap[ "Density" ]            = static_cast< uint >( IWG_Property_Type::DENSITY );
            mPropertyMap[ "BackFlowPrevention" ] = static_cast< uint >( IWG_Property_Type::BACKFLOW_PREVENTION );
            mPropertyMap[ "Select" ]             = static_cast< uint >( IWG_Property_Type::SELECT );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Pressure_Neumann::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here pressure), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the velocity FI
            Field_Interpolator* tVelocityFI =    //
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the imposed pressure property
            const std::shared_ptr< Property >& tPropPressure =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) );

            // get the imposed total pressure property
            const std::shared_ptr< Property >& tPropTotalPressure =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::TOTAL_PRESSURE ) );

            // get density property
            const std::shared_ptr< Property >& tPropDensity =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::DENSITY ) );

            // get back flow prevention prevention property
            const std::shared_ptr< Property >& tPropBackflowPrevention =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BACKFLOW_PREVENTION ) );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            real tM = 1.0;
            if ( tPropSelect != nullptr )
            {
                tM = tPropSelect->val()( 0 );

                // skip computing residual if projection matrix is zero
                if ( std::abs( tM ) < MORIS_REAL_EPS ) return;
            }

            // check that either pressure or total pressure property is set
            // FIXME: this check should moved in to checking function and not executed when computing residual
            MORIS_ERROR( tPropPressure != nullptr || tPropTotalPressure != nullptr || tPropBackflowPrevention != nullptr,
                    "IWG_Incompressible_NS_Pressure_Neumann::compute_residual - Either Pressure or TotalPressure or BackflowPrevention needs to be set." );

            // check that not both pressure and total pressure are set
            // Fixme: this check should moved in to checking function and not executed when computing residual
            MORIS_ERROR( tPropPressure == nullptr || tPropTotalPressure == nullptr,
                    "IWG_Incompressible_NS_Pressure_Neumann::set_property - Can only set either Pressure or TotalPressure." );

            // compute the residual weak form

            // when pressure is imposed
            if ( tPropPressure != nullptr )
            {
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { 0, 0 } ) +=
                        aWStar * ( tVelocityFI->N_trans() * tM * tPropPressure->val()( 0 ) * mNormal );
            }

            // when total pressure is imposed
            if ( tPropTotalPressure != nullptr )
            {
                // check that density is set
                // Fixme: this check should moved in to checking function and not executed when computing residual
                MORIS_ERROR( tPropDensity != nullptr,
                        "IWG_Incompressible_NS_Pressure_Neumann::compute_residual - Density needs to be set if TotalPressure BC is used." );

                // compute imposed pressure from prescribed total pressure and velocities
                const real tVelocNorm2      = std::pow( norm( tVelocityFI->val() ), 2.0 );
                const real tImposedPressure = tPropTotalPressure->val()( 0 ) - 0.5 * tPropDensity->val()( 0 ) * tVelocNorm2;

                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { 0, 0 } ) +=
                        aWStar * ( tM * tImposedPressure * tVelocityFI->N_trans() * mNormal );
            }

            if ( tPropBackflowPrevention != nullptr )
            {
                // compute normal velocity
                const real tNormalVelocity = dot( tVelocityFI->val(), mNormal );

                // check for back flow
                if ( tNormalVelocity < 0.0 )
                {
                    // check that density is set
                    // Fixme: this check should moved in to checking function and not executed when computing residual
                    MORIS_ERROR( tPropDensity != nullptr,
                            "IWG_Incompressible_NS_Pressure_Neumann::compute_residual - Density needs to be set if back flow prevention is used." );

                    // get the density value
                    const real tDensity = tPropDensity->val()( 0 );

                    mSet->get_residual()( 0 )(
                            { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) -=                       //
                            aWStar * (                                                                         //
                                    tM * tPropBackflowPrevention->val()( 0 ) * tDensity * tNormalVelocity *    //
                                    tVelocityFI->N_trans() * tVelocityFI->val() );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Pressure_Neumann::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Pressure_Neumann::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here pressure), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the velocity FI
            Field_Interpolator* tVelocityFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the imposed pressure property
            const std::shared_ptr< Property >& tPropPressure =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) );

            // get the imposed total pressure property
            const std::shared_ptr< Property >& tPropTotalPressure =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::TOTAL_PRESSURE ) );

            // get back flow prevention prevention property
            const std::shared_ptr< Property >& tPropBackflowPrevention =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BACKFLOW_PREVENTION ) );

            // get density property
            const std::shared_ptr< Property >& tPropDensity =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::DENSITY ) );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            real tM = 1.0;
            if ( tPropSelect != nullptr )
            {
                tM = tPropSelect->val()( 0 );

                // skip computing residual if projection matrix is zero
                if ( std::abs( tM ) < MORIS_REAL_EPS ) return;
            }

            // compute the Jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // extract sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // when pressure is imposed
                if ( tPropPressure != nullptr )
                {
                    // if imposed pressure property depends on the dof type
                    if ( tPropPressure->check_dof_dependency( tDofType ) )
                    {
                        // compute the Jacobian
                        tJac += aWStar * ( tM * tVelocityFI->N_trans() * mNormal * tPropPressure->dPropdDOF( tDofType ) );
                    }
                }

                if ( tPropTotalPressure != nullptr )
                {
                    // if dof type is residual dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // compute derivative of imposed pressure with respect to velocities
                        const Matrix< DDRMat > tdpredvel =
                                -1.0 * tPropDensity->val()( 0 ) * tVelocityFI->val_trans() * tVelocityFI->N();

                        // compute the jacobian/
                        tJac += aWStar * ( tM * tVelocityFI->N_trans() * mNormal * tdpredvel );
                    }

                    // if imposed total pressure property depends on the dof type
                    if ( tPropTotalPressure->check_dof_dependency( tDofType ) )
                    {
                        // compute the Jacobian
                        tJac += aWStar * ( tM * tVelocityFI->N_trans() * mNormal * tPropTotalPressure->dPropdDOF( tDofType ) );
                    }

                    // if density property depends on the dof type
                    if ( tPropDensity->check_dof_dependency( tDofType ) )
                    {
                        // compute derivative of imposed pressure with respect to density
                        const real tdpreddens = -0.5 * std::pow( norm( tVelocityFI->val() ), 2.0 );

                        // compute the Jacobian
                        tJac += aWStar * ( tM * tdpreddens * tVelocityFI->N_trans() * mNormal * tPropDensity->dPropdDOF( tDofType ) );
                    }
                }

                if ( tPropBackflowPrevention != nullptr )
                {
                    // compute normal velocity
                    const real tNormalVelocity = dot( tVelocityFI->val(), mNormal );

                    // check for back flow
                    if ( tNormalVelocity < 0.0 )
                    {
                        // get the density value
                        const real tDensity = tPropDensity->val()( 0 );

                        // if dof type is residual dof type
                        if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                        {
                            tJac -= aWStar * (                                                  //
                                            tPropBackflowPrevention->val()( 0 ) * tDensity *    //
                                            tVelocityFI->N_trans() * tM * (                     //
                                                    tNormalVelocity * tVelocityFI->N()          //
                                                    + tVelocityFI->val() * trans( mNormal ) * tVelocityFI->N() ) );
                        }

                        // if back flow prevention property depends on the dof type
                        if ( tPropBackflowPrevention->check_dof_dependency( tDofType ) )
                        {
                            tJac -= aWStar * (                                   //
                                            tDensity * tNormalVelocity * tM *    //
                                            tVelocityFI->N_trans() * tVelocityFI->val() * tPropBackflowPrevention->dPropdDOF( tDofType ) );
                        }

                        // if density property depends on the dof type
                        if ( tPropDensity->check_dof_dependency( tDofType ) )
                        {
                            tJac -= aWStar * (                                                              //
                                            tPropBackflowPrevention->val()( 0 ) * tNormalVelocity * tM *    //
                                            tVelocityFI->N_trans() * tVelocityFI->val() * tPropDensity->dPropdDOF( tDofType ) );
                        }
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Incompressible_NS_Pressure_Neumann::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Pressure_Neumann::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Neumann::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Pressure_Neumann::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Neumann::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
