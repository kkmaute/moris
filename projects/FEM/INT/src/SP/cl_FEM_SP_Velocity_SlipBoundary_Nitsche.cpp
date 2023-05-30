/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Velocity_SlipBoundary_Nitsche.cpp
 *
 */

#include "cl_FEM_SP_Velocity_SlipBoundary_Nitsche.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_Velocity_SlipBoundary_Nitsche::SP_Velocity_SlipBoundary_Nitsche()
                {
            // set the property pointer cell size
            mLeaderProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Viscosity" ]  = static_cast< uint >( Property_Type::VISCOSITY );
            mPropertyMap[ "Density" ]    = static_cast< uint >( Property_Type::DENSITY );
            mPropertyMap[ "SlipLength" ] = static_cast< uint >( Property_Type::SLIPLENGTH );
                }

        //------------------------------------------------------------------------------

        void SP_Velocity_SlipBoundary_Nitsche::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Leader_Follower                             aIsLeader )
        {
            // switch on leader follower
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER :
                {
                    // set dof type list
                    mLeaderDofTypes = aDofTypes;

                    // loop on dof type
                    for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if velocity
                        if( tDofString == "Velocity" )
                        {
                            mLeaderDofVelocity = tDofType;
                        }
                        else
                        {
                            // create error message
                            std::string tErrMsg =
                                    std::string( "SP_Velocity_SlipBoundary_Nitsche::set_dof_type_list - Unknown aDofString : ") +
                                    tDofString;
                            MORIS_ERROR( false , tErrMsg.c_str() );
                        }
                    }
                    break;
                }

                case mtk::Leader_Follower::FOLLOWER :
                {
                    // set dof type list
                    mFollowerDofTypes = aDofTypes;
                    break;
                }

                default:
                    MORIS_ERROR( false, "SP_Velocity_SlipBoundary_Nitsche::set_dof_type_list - unknown leader follower type." );
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Leader_Follower > > SP_Velocity_SlipBoundary_Nitsche::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void SP_Velocity_SlipBoundary_Nitsche::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofVelocity );

            // get the viscosity, density, and slip length property
            const std::shared_ptr< Property > & tPropViscosity =
                    mLeaderProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            const std::shared_ptr< Property > & tPropDensity   =
                    mLeaderProp( static_cast< uint >( Property_Type::DENSITY ) );

            const std::shared_ptr< Property > & tPropSlipLength   =
                    mLeaderProp( static_cast< uint >( Property_Type::SLIPLENGTH ) );

            // check that properties are set
            MORIS_ASSERT( tPropViscosity and tPropDensity and tPropSlipLength,
                    "SP_Velocity_Dirichlet_Nitsche::eval_SP - Viscosity, density, and slip length need to be defined.\n");

            // compute infinity norm of u
            real tInfinityNorm = std::abs( tFIVelocity->val()( 0 ) );
            for( uint iDim = 0; iDim < mSpaceDim; iDim++ )
            {
                real tAbsVelocity = std::abs( tFIVelocity->val()( iDim ) );
                if ( tInfinityNorm < tAbsVelocity )
                {
                    tInfinityNorm = tAbsVelocity;
                }
            }

            // compute deltaT
            real tDeltaT = mLeaderFIManager->get_IP_geometry_interpolator()->get_time_step();

            // set size of vector of stabilization values
            mPPVal.set_size(3,1);

            // compute stabilization parameter for normal direction
            mPPVal(0) = mParameters( 0 )( 0 ) * (
                    tPropViscosity->val()( 0 ) / tElementSize +
                    tPropDensity->val()( 0 ) * tInfinityNorm / 6.0 +
                    tPropDensity->val()( 0 ) * tElementSize / ( 12.0 * mParameters( 1 )( 0 ) * tDeltaT ) );

            // compute stabilization parameters for tangential direction
            // note: stabilization parameter mParameters( 2 )( 0 ) is 1/gamma^t in Winter et al 2018
            mPPVal(1) = mParameters( 2 )( 0 ) / ( mParameters( 2 )( 0 ) * tPropSlipLength->val()( 0 ) + tElementSize );
            mPPVal(2) = tElementSize          / ( mParameters( 2 )( 0 ) * tPropSlipLength->val()( 0 ) + tElementSize );
        }

        //------------------------------------------------------------------------------

        void SP_Velocity_SlipBoundary_Nitsche::eval_dSPdLeaderDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the dof type index
            const uint tDofIndex = mLeaderGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFIDerivative =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of dofs
            const uint tNumDofs = tFIDerivative->get_number_of_space_time_coefficients();

            // set size for dSPdLeaderDof
            mdPPdLeaderDof( tDofIndex ).set_size( 3, tNumDofs );

            // get the velocity field interpolator
            Field_Interpolator * tVelocityFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofVelocity );

            // get the viscosity, density, and slip length property
            const std::shared_ptr< Property > & tPropViscosity =
                    mLeaderProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            const std::shared_ptr< Property > & tPropDensity =
                    mLeaderProp( static_cast< uint >( Property_Type::DENSITY ) );

            const std::shared_ptr< Property > & tPropSlipLength   =
                    mLeaderProp( static_cast< uint >( Property_Type::SLIPLENGTH ) );

            // compute infinity norm and its derivative
            uint tInfinityNormIndex = 0;

            real tInfinityNorm = std::abs( tVelocityFI->val()( 0 ) );

            for( uint iDim = 0; iDim < tVelocityFI->val().numel(); iDim++ )
            {
                real tAbsVelocity = std::abs( tVelocityFI->val()( iDim ) );
                if ( tInfinityNorm < tAbsVelocity )
                {
                    tInfinityNormIndex = iDim;
                    tInfinityNorm = tAbsVelocity;
                }
            }

            // if dof type == velocity
            if( aDofTypes( 0 ) == mLeaderDofVelocity )
            {
                // compute derivative of the infinity norm
                Matrix< DDRMat > tdInfinityNormdu = tVelocityFI->N().get_row( tInfinityNormIndex );
                if( tVelocityFI->val()( tInfinityNormIndex ) < 0.0 )
                {
                    tdInfinityNormdu = -1.0 * tdInfinityNormdu;
                }

                // compute contribution from velocity
                mdPPdLeaderDof( tDofIndex ).get_row(0)=
                        mParameters( 0 )( 0 ) * tPropDensity->val()( 0 ) * tdInfinityNormdu / 6.0;

                mdPPdLeaderDof( tDofIndex ).get_row(1).fill( 0.0 );
                mdPPdLeaderDof( tDofIndex ).get_row(2).fill( 0.0 );
            }
            else
            {
                mdPPdLeaderDof( tDofIndex ).fill( 0.0 );
            }

            // if viscosity depends on dof type
            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from viscosity
                mdPPdLeaderDof( tDofIndex ).get_row(0) +=
                        mParameters( 0 )( 0 ) * tPropViscosity->dPropdDOF( aDofTypes ) / tElementSize;
            }

            // if density depends on dof type
            if( tPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute deltaT
                real tDeltaT = mLeaderFIManager->get_IP_geometry_interpolator()->get_time_step();

                // compute contribution from density
                mdPPdLeaderDof( tDofIndex ).get_row(0) +=
                        mParameters( 0 )( 0 ) * tPropDensity->dPropdDOF( aDofTypes ) *
                        ( tInfinityNorm / 6.0 + tElementSize / ( 12.0 * mParameters( 1 )( 0 ) * tDeltaT ) );
            }

            // if density depends on dof type
            if( tPropSlipLength->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false,
                        "SP_Velocity_SlipBoundary_Nitsche::eval_dSPdLeaderDOF - dof dependency of slip length not implemented.\n");
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

