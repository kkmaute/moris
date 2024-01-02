/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Compressible_Velocity_Dirichlet_Nitsche.cpp
 *
 */

#include "cl_FEM_SP_Compressible_Velocity_Dirichlet_Nitsche.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_Compressible_Velocity_Dirichlet_Nitsche::SP_Compressible_Velocity_Dirichlet_Nitsche()
        {
            // set the property pointer cell size
            mLeaderProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "DynamicViscosity" ] = static_cast< uint >( Property_Type::VISCOSITY );
        }

        //------------------------------------------------------------------------------

        void
        SP_Compressible_Velocity_Dirichlet_Nitsche::set_dof_type_list(
                Vector< Vector< MSI::Dof_Type > > &aDofTypes,
                Vector< std::string >                  &aDofStrings,
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
                                    "SP_Velocity_Dirichlet_Nitsche::set_dof_type_list - Unknown aDofString : %s",
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
                            "SP_Compressible_Velocity_Dirichlet_Nitsche::set_dof_type_list - unknown leader follower type." );
            }
        }

        //------------------------------------------------------------------------------

        Vector< std::tuple<
                fem::Measure_Type,
                mtk::Primary_Void,
                mtk::Leader_Follower > >
        SP_Compressible_Velocity_Dirichlet_Nitsche::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void
        SP_Compressible_Velocity_Dirichlet_Nitsche::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                                                std::get< 0 >( mElementSizeTuple ),
                                                std::get< 1 >( mElementSizeTuple ),
                                                std::get< 2 >( mElementSizeTuple ) )
                                        ->val()( 0 );

            // get the viscosity and density property
            std::shared_ptr< Property > &tPropViscosity =
                    mLeaderProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * tPropViscosity->val()( 0 ) / tElementSize;
        }

        //------------------------------------------------------------------------------

        void
        SP_Compressible_Velocity_Dirichlet_Nitsche::eval_dSPdLeaderDOF(
                const Vector< MSI::Dof_Type > &aDofTypes )
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
            Field_Interpolator *tFIDerivative =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdLeaderDof
            mdPPdLeaderDof( tDofIndex ).set_size( 1, tFIDerivative->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity property
            std::shared_ptr< Property > &tPropViscosity =
                    mLeaderProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            // if viscosity depends on dof type
            if ( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from viscosity
                mdPPdLeaderDof( tDofIndex ) +=
                        mParameters( 0 )( 0 ) * tPropViscosity->dPropdDOF( aDofTypes ) / tElementSize;
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
