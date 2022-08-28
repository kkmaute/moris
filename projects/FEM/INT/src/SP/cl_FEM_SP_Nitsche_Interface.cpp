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
            // set has slave flag to true
            mHasSlave = true;

            // set size for the property pointer cells
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
            mSlaveProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = static_cast< uint >( SP_Property_Type::MATERIAL );
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Master_Slave > > SP_Nitsche_Interface::get_cluster_measure_tuple_list()
        {
            return { mMasterVolumeTuple, mSlaveVolumeTuple, mInterfaceSurfaceTuple };
        }

        //------------------------------------------------------------------------------

        void SP_Nitsche_Interface::eval_SP()
        {
            // get master volume cluster measure value
            real tMasterVolume = mCluster->get_cluster_measure(
                    std::get<0>( mMasterVolumeTuple ),
                    std::get<1>( mMasterVolumeTuple ),
                    std::get<2>( mMasterVolumeTuple ) )->val()( 0 );

            // get slave volume cluster measure value
            real tSlaveVolume = mCluster->get_cluster_measure(
                    std::get<0>( mSlaveVolumeTuple ),
                    std::get<1>( mSlaveVolumeTuple ),
                    std::get<2>( mSlaveVolumeTuple ) )->val()( 0 );

            // get interface surface cluster measure value
            real tInterfaceSurface = mCluster->get_cluster_measure(
                    std::get<0>( mInterfaceSurfaceTuple ),
                    std::get<1>( mInterfaceSurfaceTuple ),
                    std::get<2>( mInterfaceSurfaceTuple ) )->val()( 0 );

            // init SP values
            mPPVal.set_size( 3, 1 );

            // get the master/slave material property
            const std::shared_ptr< Property > & tMasterPropMaterial =
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );
            const std::shared_ptr< Property > & tSlavePropMaterial =
                    mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // get the master/slave material property values
            real tMasterMaterial = tMasterPropMaterial->val()( 0 );
            real tSlaveMaterial  = tSlavePropMaterial->val()( 0 );

            // compute the mean property value
            real tDeno = ( tMasterVolume / tMasterMaterial ) + ( tSlaveVolume / tSlaveMaterial );

            // compute Nitsche SP value
            mPPVal( 0 ) = 2.0 * mParameters( 0 )( 0 ) * tInterfaceSurface / tDeno;

            // compute master weight SP value
            mPPVal( 1 ) = ( tMasterVolume / tMasterMaterial ) / tDeno;

            // compute slave weight SP value
            mPPVal( 2 ) = ( tSlaveVolume / tSlaveMaterial ) / tDeno;
        }

        //------------------------------------------------------------------------------

        void SP_Nitsche_Interface::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get master volume cluster measure value
            real tMasterVolume = mCluster->get_cluster_measure(
                    std::get<0>( mMasterVolumeTuple ),
                    std::get<1>( mMasterVolumeTuple ),
                    std::get<2>( mMasterVolumeTuple ) )->val()( 0 );

            // get slave volume cluster measure value
            real tSlaveVolume = mCluster->get_cluster_measure(
                    std::get<0>( mSlaveVolumeTuple ),
                    std::get<1>( mSlaveVolumeTuple ),
                    std::get<2>( mSlaveVolumeTuple ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdMasterDof( tDofIndex ).set_size(
                    3,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the master/slave material property
            const std::shared_ptr< Property > & tMasterPropMaterial =
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );
            const std::shared_ptr< Property > & tSlavePropMaterial =
                    mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // get the master/slave material property values
            real tMasterMaterial = tMasterPropMaterial->val()( 0 );
            real tSlaveMaterial  = tSlavePropMaterial->val()( 0 );

            // compute the mean property value
            real tDeno = ( tMasterVolume / tMasterMaterial ) + ( tSlaveVolume / tSlaveMaterial );

            // if master property depends on dof type
            if ( tMasterPropMaterial->check_dof_dependency( aDofTypes ) )
            {
                // compute Nitsche SP derivative
                mdPPdMasterDof( tDofIndex ).get_row( 0 ) =
                        this->val()( 0 ) * tMasterVolume * tMasterPropMaterial->dPropdDOF( aDofTypes ) /
                        ( std::pow( tMasterMaterial, 2 ) * tDeno );

                // compute master weight derivative
                mdPPdMasterDof( tDofIndex ).get_row( 1 ) =
                        - this->val()( 1 ) * tSlaveVolume * tMasterPropMaterial->dPropdDOF( aDofTypes ) /
                        ( tMasterMaterial * tSlaveMaterial * tDeno );

                // compute slave weight derivative
                mdPPdMasterDof( tDofIndex ).get_row( 2 ) =
                        this->val()( 2 ) * tMasterVolume * tMasterPropMaterial->dPropdDOF( aDofTypes ) /
                        ( std::pow( tMasterMaterial, 2 ) * tDeno );
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------

        void SP_Nitsche_Interface::eval_dSPdSlaveDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get master volume cluster measure value
            real tMasterVolume = mCluster->get_cluster_measure(
                    std::get<0>( mMasterVolumeTuple ),
                    std::get<1>( mMasterVolumeTuple ),
                    std::get<2>( mMasterVolumeTuple ) )->val()( 0 );

            // get slave volume cluster measure value
            real tSlaveVolume = mCluster->get_cluster_measure(
                    std::get<0>( mSlaveVolumeTuple ),
                    std::get<1>( mSlaveVolumeTuple ),
                    std::get<2>( mSlaveVolumeTuple ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mSlaveGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mSlaveFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdSlaveDof( tDofIndex ).set_size(
                    3,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the master/slave material property
            const std::shared_ptr< Property > & tMasterPropMaterial =
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );
            const std::shared_ptr< Property > & tSlavePropMaterial =
                    mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // get the master/slave material property values
            real tMasterMaterial = tMasterPropMaterial->val()( 0 );
            real tSlaveMaterial  = tSlavePropMaterial->val()( 0 );

            // compute the mean property value
            real tDeno = ( tMasterVolume / tMasterMaterial ) + ( tSlaveVolume / tSlaveMaterial );

            // if indirect dependency on the dof type
            if ( tSlavePropMaterial->check_dof_dependency( aDofTypes ) )
            {
                // compute Nitsche SP derivative
                mdPPdSlaveDof( tDofIndex ).get_row( 0 ) =
                        this->val()( 0 ) * tSlaveVolume * tSlavePropMaterial->dPropdDOF( aDofTypes ) /
                        ( std::pow( tSlaveMaterial, 2 ) * tDeno );

                // compute master weight SP derivative
                mdPPdSlaveDof( tDofIndex ).get_row( 1 ) =
                        this->val()( 1 ) * tSlaveVolume * tSlavePropMaterial->dPropdDOF( aDofTypes ) /
                        ( std::pow( tSlaveMaterial, 2 ) * tDeno );

                // compute slave weight SP derivative
                mdPPdSlaveDof( tDofIndex ).get_row( 2 ) =
                        - this->val()( 2 ) * tMasterVolume * tSlavePropMaterial->dPropdDOF( aDofTypes ) /
                        ( tMasterMaterial * tSlaveMaterial * tDeno );
            }
            else
            {
                mdPPdSlaveDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

