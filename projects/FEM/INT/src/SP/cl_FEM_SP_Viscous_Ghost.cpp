/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Viscous_Ghost.cpp
 *
 */

#include "cl_FEM_SP_Viscous_Ghost.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_Viscous_Ghost::SP_Viscous_Ghost()
        {
            // set the property pointer cell size
            mLeaderProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Viscosity" ]       = static_cast< uint >( Property_Type::VISCOSITY );
            mPropertyMap[ "InvPermeability" ] = static_cast< uint >( Property_Type::INV_PERMEABILITY );
        }

        //------------------------------------------------------------------------------

        Vector< std::tuple<
                fem::Measure_Type,
                mtk::Primary_Void,
                mtk::Leader_Follower > >
        SP_Viscous_Ghost::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void
        SP_Viscous_Ghost::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                                                std::get< 0 >( mElementSizeTuple ),
                                                std::get< 1 >( mElementSizeTuple ),
                                                std::get< 2 >( mElementSizeTuple ) )
                                        ->val()( 0 );

            // get the viscosity property
            const std::shared_ptr< Property >& tViscosityProp =
                    mLeaderProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            // get the inverse permeability property (Brinkman coefficient)
            const std::shared_ptr< Property >& tInvPermeabProp =
                    mLeaderProp( static_cast< uint >( Property_Type::INV_PERMEABILITY ) );

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * tViscosityProp->val()( 0 )    //
                   * std::pow( tElementSize, 2.0 * ( mOrder - 1.0 ) + 1.0 );

            if ( tInvPermeabProp )
            {
                mPPVal += mParameters( 0 ) * tInvPermeabProp->val()( 0 )    //
                        * std::pow( tElementSize, 2.0 * ( mOrder - 1.0 ) + 3.0 );
            }
        }

        //------------------------------------------------------------------------------

        void
        SP_Viscous_Ghost::eval_dSPdLeaderDOF(
                const Vector< MSI::Dof_Type >& aDofTypes )
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

            // get the viscosity property
            const std::shared_ptr< Property >& tViscosityProp =
                    mLeaderProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            // get the inverse permeability property (Brinkman coefficient)
            const std::shared_ptr< Property >& tInvPermeabProp =
                    mLeaderProp( static_cast< uint >( Property_Type::INV_PERMEABILITY ) );

            // if viscosity depends on dof type
            if ( tViscosityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from viscosity
                mdPPdLeaderDof( tDofIndex ) =
                        mParameters( 0 ) * std::pow( tElementSize, 2.0 * ( mOrder - 1.0 ) + 1.0 )    //
                        * tViscosityProp->dPropdDOF( aDofTypes );
            }
            else
            {
                mdPPdLeaderDof( tDofIndex ).fill( 0.0 );
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
