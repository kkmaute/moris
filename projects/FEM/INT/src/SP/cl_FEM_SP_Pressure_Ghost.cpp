/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Pressure_Ghost.cpp
 *
 */

#include "cl_FEM_SP_Pressure_Ghost.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_Pressure_Ghost::SP_Pressure_Ghost()
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
        SP_Pressure_Ghost::set_parameters( moris::Vector< Matrix< DDRMat > > aParameters )
        {
            // FIXME not necessary
            // set mParameters
            mParameters = aParameters;

            // get number of parameters
            uint tParamSize = aParameters.size();

            // check for proper size of constant function parameters
            MORIS_ERROR( tParamSize >= 1 && tParamSize < 3,
                    "SP_Pressure_Ghost::set_parameters - either 1 or 2 constant parameters need to be set." );

            // set alphaP
            mAlphaP = aParameters( 0 )( 0 );

            // set alpha P flag to true
            mSetAlphaP = true;

            // if time term
            if ( tParamSize > 1 )
            {
                // set theta
                mTheta = aParameters( 1 )( 0 );

                // set theta flag to true
                mSetTheta = true;
            }
        }

        //------------------------------------------------------------------------------

        void
        SP_Pressure_Ghost::set_dof_type_list(
                moris::Vector< moris::Vector< MSI::Dof_Type > >& aDofTypes,
                moris::Vector< std::string >&                  aDofStrings,
                mtk::Leader_Follower                            aIsLeader )
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
                            // error unknown dof string
                            MORIS_ERROR( false,
                                    "SP_Pressure_Ghost::set_dof_type_list - Unknown aDofString : %s \n",
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
                    MORIS_ERROR( false, "SP_Pressure_Ghost::set_dof_type_list - unknown leader follower type." );
            }
        }

        //------------------------------------------------------------------------------

        moris::Vector< std::tuple<
                fem::Measure_Type,
                mtk::Primary_Void,
                mtk::Leader_Follower > >
        SP_Pressure_Ghost::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void
        SP_Pressure_Ghost::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                                                std::get< 0 >( mElementSizeTuple ),
                                                std::get< 1 >( mElementSizeTuple ),
                                                std::get< 2 >( mElementSizeTuple ) )
                                        ->val()( 0 );

            // get the viscosity and density property
            const std::shared_ptr< Property >& tViscosityProp =
                    mLeaderProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            const std::shared_ptr< Property >& tDensityProp =
                    mLeaderProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the inverse permeability property (Brinkman coefficient)
            const std::shared_ptr< Property >& tInvPermeabProp =
                    mLeaderProp( static_cast< uint >( Property_Type::INV_PERMEABILITY ) );

            // get the velocity FI
            Field_Interpolator* tVelocityFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofVelocity );

            // compute infinity norm of u
            real tInfinityNorm = std::abs( tVelocityFI->val()( 0 ) );
            for ( uint iDim = 0; iDim < tVelocityFI->val().numel(); iDim++ )
            {
                real tAbsVelocity = std::abs( tVelocityFI->val()( iDim ) );
                if ( tInfinityNorm < tAbsVelocity )
                {
                    tInfinityNorm = tAbsVelocity;
                }
            }

            // compute deltaP
            real tDeltaP = tViscosityProp->val()( 0 ) / tElementSize    //
                         + tDensityProp->val()( 0 ) * tInfinityNorm / 6.0;

            // add contribution from inverse permeability
            if ( tInvPermeabProp )
            {
                tDeltaP += tInvPermeabProp->val()( 0 ) * tElementSize / 12.0;
            }

            // add time step contribution
            if ( mSetTheta )
            {
                // compute deltaT
                real tDeltaT = mLeaderFIManager->get_IP_geometry_interpolator()->get_time_step();

                // add time step contribution
                tDeltaP += tDensityProp->val()( 0 ) * tElementSize / ( 12.0 * mTheta * tDeltaT );
            }

            // compute stabilization parameter value
            mPPVal = mAlphaP * std::pow( tElementSize, 2 * mOrder ) / tDeltaP;
        }

        //------------------------------------------------------------------------------

        void
        SP_Pressure_Ghost::eval_dSPdLeaderDOF(
                const moris::Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                                                std::get< 0 >( mElementSizeTuple ),
                                                std::get< 1 >( mElementSizeTuple ),
                                                std::get< 2 >( mElementSizeTuple ) )
                                        ->val()( 0 );

            // get the dof type index
            uint tDofIndex = mLeaderGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator* tFI = mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdLeaderDof
            mdPPdLeaderDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // get velocity field interpolator
            Field_Interpolator* tVelocityFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mLeaderDofVelocity );

            // get the viscosity property
            const std::shared_ptr< Property >& tViscosityProp =
                    mLeaderProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            // get the density property
            const std::shared_ptr< Property >& tDensityProp =
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

            // compute deltaP
            real tDeltaP = tViscosityProp->val()( 0 ) / tElementSize    //
                         + tDensityProp->val()( 0 ) * tInfinityNorm / 6.0;

            // add contribution from inverse permeability
            if ( tInvPermeabProp )
            {
                tDeltaP += tInvPermeabProp->val()( 0 ) * tElementSize / 12.0;
            }

            // if time step contribution
            if ( mSetTheta )
            {
                // compute deltaT
                real tDeltaT = mLeaderFIManager->get_IP_geometry_interpolator()->get_time_step();

                // add time step contribution
                tDeltaP += tDensityProp->val()( 0 ) * tElementSize / ( 12.0 * mTheta * tDeltaT );
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
                        -mAlphaP * std::pow( tElementSize, 2 * mOrder ) * tDensityProp->val()( 0 )    //
                        * tdInfinityNormdu / ( 6.0 * std::pow( tDeltaP, 2.0 ) );
            }
            else
            {
                mdPPdLeaderDof( tDofIndex ).fill( 0.0 );
            }

            // if viscosity depends on dof type
            if ( tViscosityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from viscosity
                mdPPdLeaderDof( tDofIndex ) -=
                        mAlphaP * std::pow( tElementSize, 2 * mOrder ) * tViscosityProp->dPropdDOF( aDofTypes )    //
                        / ( tElementSize * std::pow( tDeltaP, 2.0 ) );
            }

            // if density depends on dof type
            if ( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from density
                mdPPdLeaderDof( tDofIndex ) -=
                        mAlphaP * std::pow( tElementSize, 2 * mOrder ) * tDensityProp->dPropdDOF( aDofTypes )    //
                        * ( tInfinityNorm / 6.0 ) / std::pow( tDeltaP, 2.0 );

                // if time step contribution
                if ( mSetTheta )
                {
                    // compute deltaT
                    real tDeltaT = mLeaderFIManager->get_IP_geometry_interpolator()->get_time_step();

                    // add time step contribution
                    mdPPdLeaderDof( tDofIndex ) -=
                            mAlphaP * std::pow( tElementSize, 2 * mOrder ) * tDensityProp->dPropdDOF( aDofTypes )    //
                            * ( tElementSize / ( 12.0 * mTheta * tDeltaT ) ) / std::pow( tDeltaP, 2.0 );
                }
            }

            // if inverse permeability depends on dof type
            if ( tInvPermeabProp )
            {
                if ( tInvPermeabProp->check_dof_dependency( aDofTypes ) )
                {
                    MORIS_ERROR( false, "dof dependence of inverse permeability not implemented." );
                }
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
