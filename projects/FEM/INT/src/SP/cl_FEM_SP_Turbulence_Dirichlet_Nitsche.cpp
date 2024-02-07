/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Turbulence_Dirichlet_Nitsche.cpp
 *
 */

#include "cl_FEM_SP_Turbulence_Dirichlet_Nitsche.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        SP_Turbulence_Dirichlet_Nitsche::SP_Turbulence_Dirichlet_Nitsche()
        {
            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "SpalartAllmarasTurbulence" ] = static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        }

        //------------------------------------------------------------------------------

        void SP_Turbulence_Dirichlet_Nitsche::set_dof_type_list(
                Vector< Vector< MSI::Dof_Type > > & aDofTypes,
                Vector< std::string >                  & aDofStrings,
                mtk::Leader_Follower                             aIsLeader )
        {
            Stabilization_Parameter::set_dof_type_list( aDofTypes, aIsLeader );
        }

        //------------------------------------------------------------------------------

        Vector< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Leader_Follower > > SP_Turbulence_Dirichlet_Nitsche::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void SP_Turbulence_Dirichlet_Nitsche::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMSATurbulence =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * tCMSATurbulence->diffusion_coefficient() / tElementSize;
        }

        //------------------------------------------------------------------------------

        void SP_Turbulence_Dirichlet_Nitsche::eval_dSPdLeaderDOF(
                const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mLeaderGlobalDofTypeMap( tDofType );

            // get FI for derivative dof type
            Field_Interpolator * tFIDer =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdLeaderDof( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMSATurbulence =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // if turbulence CM depends on dof
            if( tCMSATurbulence->check_dof_dependency( aDofTypes ) )
            {
                // add contribution from diffusion coefficient
                mdPPdLeaderDof( tDofIndex ) =
                        mParameters( 0 ) * tCMSATurbulence->ddiffusioncoeffdu( aDofTypes ) / tElementSize;
            }
            else
            {
                mdPPdLeaderDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

