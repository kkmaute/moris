/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Convective_Ghost.cpp
 *
 */

#include "cl_FEM_SP_Convective_Ghost.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_Convective_Ghost::SP_Convective_Ghost()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Density" ] = static_cast< uint >( Property_Type::DENSITY );
        }

        //------------------------------------------------------------------------------

        void SP_Convective_Ghost::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Master_Slave                             aIsMaster )
        {
            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    // set dof type list
                    mMasterDofTypes = aDofTypes;

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
                            mMasterDofVelocity = tDofType;
                        }
                        else
                        {
                            // error unknown dof string
                            MORIS_ERROR( false ,
                                    "SP_Convective_Ghost::set_dof_type_list - Unknown aDofString : %s \n",
                                    tDofString.c_str() );
                        }
                    }
                    break;
                }

                case mtk::Master_Slave::SLAVE :
                {
                    // set dof type list
                    mSlaveDofTypes = aDofTypes;
                    break;
                }

                default:
                    MORIS_ERROR( false, "SP_Convective_Ghost::set_dof_type_list - unknown master slave type." );
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Master_Slave > > SP_Convective_Ghost::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void SP_Convective_Ghost::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the velocity FI
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the density property
            const std::shared_ptr< Property > & tDensityProp =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get absolute value of u.n
            real tNormalDispl = dot( tVelocityFI->val(), mNormal );
            real tAbsReal = std::sqrt( tNormalDispl * tNormalDispl + mEpsilon * mEpsilon ) - mEpsilon;

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * tDensityProp->val()( 0 ) * tAbsReal * std::pow( tElementSize, 2.0 );
        }

        //------------------------------------------------------------------------------

        void SP_Convective_Ghost::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof
            mdPPdMasterDof( tDofIndex ).set_size(
                    1,
                    tFI->get_number_of_space_time_coefficients() );

            // get the velocity FI
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get absolute value of u.n
            real tNormalDispl = dot( tVelocityFI->val(), mNormal );
            real tAbsReal = std::sqrt( tNormalDispl * tNormalDispl + mEpsilon * mEpsilon ) - mEpsilon;

            // get the density property
            const std::shared_ptr< Property > & tDensityProp =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // if velocity dof
            if( aDofTypes( 0 ) == mMasterDofVelocity )
            {
                // compute contribution from velocity
                mdPPdMasterDof( tDofIndex ) =
                        mParameters( 0 ) * std::pow( tElementSize, 2.0 ) * tDensityProp->val()( 0 ) *
                        tNormalDispl * trans( mNormal ) * tVelocityFI->N() /
                        ( tAbsReal + mEpsilon );
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }

            // if density depends on dof
            if( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from density
                mdPPdMasterDof( tDofIndex ) +=
                        mParameters( 0 ) * std::pow( tElementSize, 2.0 ) * tAbsReal *
                        tDensityProp->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

