/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Nitsche_Interface.cpp
 *
 */

#include "cl_FEM_SP_Nitsche_Interface.hpp" //FEM/INT/src
#include "cl_FEM_Cluster.hpp"              //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"              //FEM/INT/src

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

        SP_Nitsche_Interface::SP_Nitsche_Interface()
        {
            // set has follower flag to true
            mHasFollower = true;

            // set size for the property pointer cells
            mLeaderProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
            mFollowerProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = static_cast< uint >( SP_Property_Type::MATERIAL );
        }

        //------------------------------------------------------------------------------

        Vector< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Leader_Follower > > SP_Nitsche_Interface::get_cluster_measure_tuple_list()
        {
            return { mLeaderVolumeTuple, mFollowerVolumeTuple, mInterfaceSurfaceTuple };
        }

        //------------------------------------------------------------------------------

        void SP_Nitsche_Interface::eval_SP()
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
            mPPVal.set_size( 3, 1 );

            // get the leader/follower material property
            const std::shared_ptr< Property > & tLeaderPropMaterial =
                    mLeaderProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );
            const std::shared_ptr< Property > & tFollowerPropMaterial =
                    mFollowerProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // get the leader/follower material property values
            real tLeaderMaterial = tLeaderPropMaterial->val()( 0 );
            real tFollowerMaterial  = tFollowerPropMaterial->val()( 0 );

            // compute the mean property value
            real tDeno = ( tLeaderVolume / tLeaderMaterial ) + ( tFollowerVolume / tFollowerMaterial );

            // compute Nitsche SP value
            mPPVal( 0 ) = 2.0 * mParameters( 0 )( 0 ) * tInterfaceSurface / tDeno;

            // compute leader weight SP value
            mPPVal( 1 ) = ( tLeaderVolume / tLeaderMaterial ) / tDeno;

            // compute follower weight SP value
            mPPVal( 2 ) = ( tFollowerVolume / tFollowerMaterial ) / tDeno;
        }

        //------------------------------------------------------------------------------

        void SP_Nitsche_Interface::eval_dSPdLeaderDOF(
                const Vector< MSI::Dof_Type > & aDofTypes )
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

            // get the leader/follower material property
            const std::shared_ptr< Property > & tLeaderPropMaterial =
                    mLeaderProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );
            const std::shared_ptr< Property > & tFollowerPropMaterial =
                    mFollowerProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // get the leader/follower material property values
            real tLeaderMaterial = tLeaderPropMaterial->val()( 0 );
            real tFollowerMaterial  = tFollowerPropMaterial->val()( 0 );

            // compute the mean property value
            real tDeno = ( tLeaderVolume / tLeaderMaterial ) + ( tFollowerVolume / tFollowerMaterial );

            // if leader property depends on dof type
            if ( tLeaderPropMaterial->check_dof_dependency( aDofTypes ) )
            {
                // compute Nitsche SP derivative
                mdPPdLeaderDof( tDofIndex ).get_row( 0 ) =
                        this->val()( 0 ) * tLeaderVolume * tLeaderPropMaterial->dPropdDOF( aDofTypes ) /
                        ( std::pow( tLeaderMaterial, 2 ) * tDeno );

                // compute leader weight derivative
                mdPPdLeaderDof( tDofIndex ).get_row( 1 ) =
                        - this->val()( 1 ) * tFollowerVolume * tLeaderPropMaterial->dPropdDOF( aDofTypes ) /
                        ( tLeaderMaterial * tFollowerMaterial * tDeno );

                // compute follower weight derivative
                mdPPdLeaderDof( tDofIndex ).get_row( 2 ) =
                        this->val()( 2 ) * tLeaderVolume * tLeaderPropMaterial->dPropdDOF( aDofTypes ) /
                        ( std::pow( tLeaderMaterial, 2 ) * tDeno );
            }
            else
            {
                mdPPdLeaderDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------

        void SP_Nitsche_Interface::eval_dSPdFollowerDOF(
                const Vector< MSI::Dof_Type > & aDofTypes )
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
            mdPPdFollowerDof( tDofIndex ).set_size(
                    3,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the leader/follower material property
            const std::shared_ptr< Property > & tLeaderPropMaterial =
                    mLeaderProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );
            const std::shared_ptr< Property > & tFollowerPropMaterial =
                    mFollowerProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // get the leader/follower material property values
            real tLeaderMaterial = tLeaderPropMaterial->val()( 0 );
            real tFollowerMaterial  = tFollowerPropMaterial->val()( 0 );

            // compute the mean property value
            real tDeno = ( tLeaderVolume / tLeaderMaterial ) + ( tFollowerVolume / tFollowerMaterial );

            // if indirect dependency on the dof type
            if ( tFollowerPropMaterial->check_dof_dependency( aDofTypes ) )
            {
                // compute Nitsche SP derivative
                mdPPdFollowerDof( tDofIndex ).get_row( 0 ) =
                        this->val()( 0 ) * tFollowerVolume * tFollowerPropMaterial->dPropdDOF( aDofTypes ) /
                        ( std::pow( tFollowerMaterial, 2 ) * tDeno );

                // compute leader weight SP derivative
                mdPPdFollowerDof( tDofIndex ).get_row( 1 ) =
                        this->val()( 1 ) * tFollowerVolume * tFollowerPropMaterial->dPropdDOF( aDofTypes ) /
                        ( std::pow( tFollowerMaterial, 2 ) * tDeno );

                // compute follower weight SP derivative
                mdPPdFollowerDof( tDofIndex ).get_row( 2 ) =
                        - this->val()( 2 ) * tLeaderVolume * tFollowerPropMaterial->dPropdDOF( aDofTypes ) /
                        ( tLeaderMaterial * tFollowerMaterial * tDeno );
            }
            else
            {
                mdPPdFollowerDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

