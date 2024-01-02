/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Velocity_Dirichlet_Nitsche.cpp
 *
 */

#include "cl_FEM_SP_Velocity_Dirichlet_Nitsche.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        SP_Velocity_Dirichlet_Nitsche::SP_Velocity_Dirichlet_Nitsche()
        {
            // set the property pointer cell size
            mLeaderProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Viscosity" ]       = static_cast< uint >( Property_Type::VISCOSITY );
            mPropertyMap[ "Density" ]         = static_cast< uint >( Property_Type::DENSITY );
            mPropertyMap[ "InvPermeability" ] = static_cast< uint >( Property_Type::INV_PERMEABILITY );
        }

        //------------------------------------------------------------------------------

        void
        SP_Velocity_Dirichlet_Nitsche::set_parameters( Vector< Matrix< DDRMat > > aParameters )
        {
            // FIXME not necessary
            // set mParameters
            mParameters = aParameters;

            // get number of parameters
            uint tParamSize = aParameters.size();

            // check for proper size of constant function parameters
            MORIS_ERROR( tParamSize >= 1 && tParamSize <= 3,
                    "SP_Velocity_Dirichlet_Nitsche::set_parameters - either 1, 2, or 3 constant parameters need to be set." );

            // set alphaN
            mAlphaN = aParameters( 0 )( 0 );

            // set alpha N flag to true
            mSetAlphaN = true;

            // if time term
            if ( tParamSize > 1 )
            {
                // set alpha_time
                mAlphaTime = aParameters( 1 )( 0 );

                // set alpha_time flag to true only if mAlphaTime > 0.0
                mSetAlphaTime = mAlphaTime > MORIS_REAL_EPS ? true : false;
            }

            // if geometry measure type is defined
            if ( tParamSize == 3 )
            {
                mGeometryFormulation = mParameters( 2 )( 0 );
            }

            // set pointer function for geometry measure
            switch ( mGeometryFormulation )
            {
                case 0:
                {
                    mGeometryMeasureFunc = &SP_Velocity_Dirichlet_Nitsche::eval_geometry_measure_edge_length;
                    break;
                }
                case 1:
                {
                    mGeometryMeasureFunc = &SP_Velocity_Dirichlet_Nitsche::eval_geometry_measure_relative_edge_length;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "SP_Velocity_Dirichlet_Nitsche::set_parameters - wrong required formulation of cluster measure." );
                    break;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        SP_Velocity_Dirichlet_Nitsche::set_dof_type_list(
                Vector< Vector< MSI::Dof_Type > >& aDofTypes,
                Vector< std::string >&                  aDofStrings,
                mtk::Leader_Follower                         aIsLeader )
        {
            // switch on leader follower
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    // set dof type list
                    mLeaderDofTypes = aDofTypes;

                    // loop on dof type
                    for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if velocity
                        if ( tDofString == "Velocity" )
                        {
                            mLeaderDofVelocity = tDofType;
                        }
                        else
                        {
                            // create error message
                            MORIS_ERROR( false,
                                    " SP_Velocity_Dirichlet_Nitsche::set_dof_type_list - Unknown aDofString : %s",
                                    tDofString.c_str() );
                        }
                    }
                    break;
                }

                case mtk::Leader_Follower::FOLLOWER:
                {
                    // set dof type list
                    mFollowerDofTypes = aDofTypes;
                    break;
                }

                default:
                    MORIS_ERROR( false,
                            "SP_Velocity_Dirichlet_Nitsche::set_dof_type_list - unknown leader follower type." );
            }
        }

        //------------------------------------------------------------------------------

        Vector< std::tuple<
                fem::Measure_Type,
                mtk::Primary_Void,
                mtk::Leader_Follower > >
        SP_Velocity_Dirichlet_Nitsche::get_cluster_measure_tuple_list()
        {
            if ( mGeometryFormulation == 1 )
            {
                return { mLeaderVolumeTuple, mInterfaceSurfaceTuple };
            }
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void
        SP_Velocity_Dirichlet_Nitsche::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = ( this->*mGeometryMeasureFunc )();

            // get the velocity FI
            Field_Interpolator* tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofVelocity );

            // get the viscosity and density property
            const std::shared_ptr< Property >& tPropViscosity =
                    mLeaderProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            const std::shared_ptr< Property >& tPropDensity =
                    mLeaderProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the inverse permeability property (Brinkman coefficient)
            const std::shared_ptr< Property >& tInvPermeabProp =
                    mLeaderProp( static_cast< uint >( Property_Type::INV_PERMEABILITY ) );

            // check that properties are set
            MORIS_ASSERT( tPropViscosity and tPropDensity,
                    "SP_Velocity_Dirichlet_Nitsche::eval_SP - Viscosity and density need to be defined.\n" );

            // compute infinity norm of u
            real tInfinityNorm = std::abs( tFIVelocity->val()( 0 ) );

            for ( uint iDim = 0; iDim < mSpaceDim; iDim++ )
            {
                real tAbsVelocity = std::abs( tFIVelocity->val()( iDim ) );
                if ( tInfinityNorm < tAbsVelocity )
                {
                    tInfinityNorm = tAbsVelocity;
                }
            }

            // compute stabilization parameter value
            mPPVal = mAlphaN * ( tPropViscosity->val()( 0 ) / tElementSize + tPropDensity->val()( 0 ) * tInfinityNorm / 6.0 );

            // add contribution from inverse permeability
            if ( tInvPermeabProp )
            {
                mPPVal += mAlphaN * tInvPermeabProp->val()( 0 ) * tElementSize / 12.0;
            }

            // if time step contribution
            if ( mSetAlphaTime )
            {
                // compute deltaT
                real tDeltaT = mLeaderFIManager->get_IP_geometry_interpolator()->get_time_step();

                // add time step contribution
                mPPVal += mAlphaN * tPropDensity->val()( 0 ) * tElementSize / ( 12.0 * mAlphaTime * tDeltaT );
            }
        }

        //------------------------------------------------------------------------------

        void
        SP_Velocity_Dirichlet_Nitsche::eval_dSPdLeaderDOF(
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get element size cluster measure value
            real tElementSize = ( this->*mGeometryMeasureFunc )();

            // get the dof type index
            uint tDofIndex = mLeaderGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator* tFIDerivative =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdLeaderDof
            mdPPdLeaderDof( tDofIndex ).set_size( 1, tFIDerivative->get_number_of_space_time_coefficients() );

            // get the velocity field interpolator
            Field_Interpolator* tVelocityFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofVelocity );

            // get the viscosity property
            const std::shared_ptr< Property >& tPropViscosity =
                    mLeaderProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            // get the density property
            const std::shared_ptr< Property >& tPropDensity =
                    mLeaderProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the inverse permeability property (Brinkman coefficient)
            const std::shared_ptr< Property >& tInvPermeabProp =
                    mLeaderProp( static_cast< uint >( Property_Type::INV_PERMEABILITY ) );

            // compute infinity norm
            uint tInfinityNormIndex = 0;
            real tInfinityNorm      = std::abs( tVelocityFI->val()( 0 ) );
            for ( uint iDim = 0; iDim < tVelocityFI->val().numel(); iDim++ )
            {
                real tAbsVelocity = std::abs( tVelocityFI->val()( iDim ) );
                if ( tInfinityNorm < tAbsVelocity )
                {
                    tInfinityNormIndex = iDim;
                    tInfinityNorm      = tAbsVelocity;
                }
            }

            // if dof type == velocity
            if ( aDofTypes( 0 ) == mLeaderDofVelocity )
            {
                // compute derivative of the infinity norm
                Matrix< DDRMat > tdInfinityNormdu = tVelocityFI->N().get_row( tInfinityNormIndex );
                if ( tVelocityFI->val()( tInfinityNormIndex ) < 0.0 )
                {
                    tdInfinityNormdu = -1.0 * tdInfinityNormdu;
                }

                // compute contribution from velocity
                mdPPdLeaderDof( tDofIndex ) =
                        mAlphaN * tPropDensity->val()( 0 ) * tdInfinityNormdu / 6.0;
            }
            else
            {
                mdPPdLeaderDof( tDofIndex ).fill( 0.0 );
            }

            // if viscosity depends on dof type
            if ( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from viscosity
                mdPPdLeaderDof( tDofIndex ) +=
                        mAlphaN * tPropViscosity->dPropdDOF( aDofTypes ) / tElementSize;
            }

            // if inverse permeability depends on dof type
            if ( tInvPermeabProp )
            {
                if ( tInvPermeabProp->check_dof_dependency( aDofTypes ) )
                {
                    MORIS_ERROR( false, "dof dependence of inverse permeability not implemented." );
                }
            }

            // if density depends on dof type
            if ( tPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from density
                mdPPdLeaderDof( tDofIndex ) +=
                        mAlphaN * tPropDensity->dPropdDOF( aDofTypes ) * ( tInfinityNorm / 6.0 );

                // if time step contribution
                if ( mSetAlphaTime )
                {
                    // compute deltaT
                    real tDeltaT = mLeaderFIManager->get_IP_geometry_interpolator()->get_time_step();

                    // add time step contribution
                    mdPPdLeaderDof( tDofIndex ) +=
                            mAlphaN * tPropDensity->dPropdDOF( aDofTypes ) *    //
                            ( tElementSize / ( 12.0 * mAlphaTime * tDeltaT ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        real
        SP_Velocity_Dirichlet_Nitsche::eval_geometry_measure_edge_length()
        {
            return mCluster->get_cluster_measure(
                                   std::get< 0 >( mElementSizeTuple ),
                                   std::get< 1 >( mElementSizeTuple ),
                                   std::get< 2 >( mElementSizeTuple ) )
                    ->val()( 0 );
        }

        //------------------------------------------------------------------------------

        real
        SP_Velocity_Dirichlet_Nitsche::eval_geometry_measure_relative_edge_length()
        {
            // get leader volume cluster measure value
            const real tLeaderVolume =    //
                    mCluster->get_cluster_measure(
                                    std::get< 0 >( mLeaderVolumeTuple ),
                                    std::get< 1 >( mLeaderVolumeTuple ),
                                    std::get< 2 >( mLeaderVolumeTuple ) )
                            ->val()( 0 );

            // get interface surface cluster measure value
            const real tInterfaceSurface =    //
                    mCluster->get_cluster_measure(
                                    std::get< 0 >( mInterfaceSurfaceTuple ),
                                    std::get< 1 >( mInterfaceSurfaceTuple ),
                                    std::get< 2 >( mInterfaceSurfaceTuple ) )
                            ->val()( 0 );

            return tLeaderVolume / tInterfaceSurface;
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
