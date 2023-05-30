/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.cpp
 *
 */

#include "cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        SP_Spalart_Allmaras_Nitsche_Interface::SP_Spalart_Allmaras_Nitsche_Interface()
        {
            // set has follower flag to true
            mHasFollower = true;

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mFollowerCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "SpalartAllmarasTurbulence" ] = static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Leader_Follower                             aIsLeader )
        {
            Stabilization_Parameter::set_dof_type_list( aDofTypes, aIsLeader );
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Leader_Follower > > SP_Spalart_Allmaras_Nitsche_Interface::get_cluster_measure_tuple_list()
        {
            return { mLeaderVolumeTuple, mFollowerVolumeTuple, mInterfaceSurfaceTuple };
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::eval_SP()
        {
            // get leader volume cluster measure value
            real tLeaderVolume = mCluster->get_cluster_measure(
                    std::get<0>( mLeaderVolumeTuple ),
                    std::get<1>( mLeaderVolumeTuple ),
                    std::get<2>( mLeaderVolumeTuple ) )->val()( 0 );

            // get follower volume cluster measure value
            real tFollowerVolume = mCluster->get_cluster_measure(
                    std::get<0>( mFollowerVolumeTuple ),
                    std::get<1>( mFollowerVolumeTuple ),
                    std::get<2>( mFollowerVolumeTuple ) )->val()( 0 );

            // get interface surface cluster measure value
            real tInterfaceSurface = mCluster->get_cluster_measure(
                    std::get<0>( mInterfaceSurfaceTuple ),
                    std::get<1>( mInterfaceSurfaceTuple ),
                    std::get<2>( mInterfaceSurfaceTuple ) )->val()( 0 );

            // init SP values
            mPPVal.set_size( 3, 1, 0.0 );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMLeaderSATurbulence =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );
            const std::shared_ptr< Constitutive_Model > & tCMFollowerSATurbulence =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // compute weighted property
            real tDeno =
                    tLeaderVolume / tCMLeaderSATurbulence->diffusion_coefficient()( 0 ) +
                    tFollowerVolume / tCMFollowerSATurbulence->diffusion_coefficient()( 0 );

            // compute Nitsche parameter value
            mPPVal( 0 ) = 2.0 * mParameters( 0 )( 0 ) * tInterfaceSurface / tDeno;

            // compute leader weight value
            mPPVal( 1 ) = ( tLeaderVolume / tCMLeaderSATurbulence->diffusion_coefficient()( 0 ) ) / tDeno;

            // compute follower weight value
            mPPVal( 2 ) = ( tFollowerVolume / tCMFollowerSATurbulence->diffusion_coefficient()( 0 ) ) / tDeno;
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::eval_dSPdLeaderDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get leader volume cluster measure value
            real tLeaderVolume = mCluster->get_cluster_measure(
                    std::get<0>( mLeaderVolumeTuple ),
                    std::get<1>( mLeaderVolumeTuple ),
                    std::get<2>( mLeaderVolumeTuple ) )->val()( 0 );

            // get follower volume cluster measure value
            real tFollowerVolume = mCluster->get_cluster_measure(
                    std::get<0>( mFollowerVolumeTuple ),
                    std::get<1>( mFollowerVolumeTuple ),
                    std::get<2>( mFollowerVolumeTuple ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mLeaderGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdLeaderDof( tDofIndex ).set_size(
                    3,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMLeaderSATurbulence =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );
            const std::shared_ptr< Constitutive_Model > & tCMFollowerSATurbulence =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // if leader turbulence CM depends on dof type
            if ( tCMLeaderSATurbulence->check_dof_dependency( aDofTypes ) )
            {
                // compute weighted property
                real tDeno =
                        tLeaderVolume / tCMLeaderSATurbulence->diffusion_coefficient()( 0 ) +
                        tFollowerVolume / tCMFollowerSATurbulence->diffusion_coefficient()( 0 );

                // compute derivative of Nitsche SP FIXME ????
                mdPPdLeaderDof( tDofIndex ).get_row( 0 ) =
                        this->val()( 0 ) * tLeaderVolume * tCMLeaderSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( std::pow( tCMLeaderSATurbulence->diffusion_coefficient()( 0 ), 2 ) * tDeno );

                // compute derivative of leader weight SP
                mdPPdLeaderDof( tDofIndex ).get_row( 1 ) =
                        - this->val()( 1 ) * tFollowerVolume * tCMLeaderSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( tCMLeaderSATurbulence->diffusion_coefficient()( 0 ) * tCMFollowerSATurbulence->diffusion_coefficient()( 0 ) * tDeno );

                // compute derivative of the follower weight SP
                mdPPdLeaderDof( tDofIndex ).get_row( 2 ) =
                        this->val()( 2 ) * tLeaderVolume * tCMLeaderSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( std::pow( tCMLeaderSATurbulence->diffusion_coefficient()( 0 ), 2 ) * tDeno );
            }
            else
            {
                mdPPdLeaderDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::eval_dSPdFollowerDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get leader volume cluster measure value
            real tLeaderVolume = mCluster->get_cluster_measure(
                    std::get<0>( mLeaderVolumeTuple ),
                    std::get<1>( mLeaderVolumeTuple ),
                    std::get<2>( mLeaderVolumeTuple ) )->val()( 0 );

            // get follower volume cluster measure value
            real tFollowerVolume = mCluster->get_cluster_measure(
                    std::get<0>( mFollowerVolumeTuple ),
                    std::get<1>( mFollowerVolumeTuple ),
                    std::get<2>( mFollowerVolumeTuple ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mFollowerGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mFollowerFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdFollowerDof( tDofIndex ).set_size( 3, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMLeaderSATurbulence =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );
            const std::shared_ptr< Constitutive_Model > & tCMFollowerSATurbulence =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // if follower turbulence CM depends on dof type
            if ( tCMFollowerSATurbulence->check_dof_dependency( aDofTypes ) )
            {
                // compute weighted property
                real tDeno =
                        tLeaderVolume / tCMLeaderSATurbulence->diffusion_coefficient()( 0 ) +
                        tFollowerVolume / tCMFollowerSATurbulence->diffusion_coefficient()( 0 );

                // compute derivative of Nitsche SP FIXME ????
                mdPPdFollowerDof( tDofIndex ).get_row( 0 ) =
                        this->val()( 0 ) * tFollowerVolume * tCMFollowerSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( std::pow( tCMFollowerSATurbulence->diffusion_coefficient()( 0 ), 2 ) * tDeno );

                // compute derivative of leader weight SP
                mdPPdFollowerDof( tDofIndex ).get_row( 1 ) =
                        this->val()( 1 ) * tFollowerVolume * tCMFollowerSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( std::pow( tCMFollowerSATurbulence->diffusion_coefficient()( 0 ), 2 ) * tDeno );

                // compute derivative of the follower weight SP
                mdPPdFollowerDof( tDofIndex ).get_row( 2 ) =
                        - this->val()( 2 ) * tLeaderVolume * tCMFollowerSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                        ( tCMLeaderSATurbulence->diffusion_coefficient()( 0 ) * tCMFollowerSATurbulence->diffusion_coefficient()( 0 ) * tDeno );
            }
            else
            {
                mdPPdFollowerDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

