/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_FS_Struc_Interface.cpp
 *
 */

#include "cl_FEM_IWG_FS_Struc_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_FS_Struc_Interface::IWG_FS_Struc_Interface(){}

        //------------------------------------------------------------------------------

        void IWG_FS_Struc_Interface::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get leader field interpolator for the residual dof type
            // FIXME protect dof type
            Field_Interpolator * tFIFluidPressure =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::P );

            // get leader field interpolator for the residual dof type
            Field_Interpolator * tFISolidDispl =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // compute follower residual
            mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex },
                    { 0, 0 } ) -= aWStar * (
                            trans( tFISolidDispl->N() ) * tFIFluidPressure->val() );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_FS_Struc_Interface::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_FS_Struc_Interface::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get leader field interpolator for the residual dof type
            // FIXME protect dof type
            Field_Interpolator * tFIFluidPressure =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::P );

            // get leader field interpolator for the residual dof type
            Field_Interpolator * tFISolidDispl =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get number of leader dof dependencies
            uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through leader constitutive models
            for( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == MSI::Dof_Type::P )
                {
                    mSet->get_jacobian()(
                            { tFollowerResStartIndex,  tFollowerResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * (
                                    trans( tFISolidDispl->N() ) * tFIFluidPressure->N() );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_FS_Struc_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_FS_Struc_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_FS_Struc_Interface::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_FS_Struc_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_FS_Struc_Interface::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

