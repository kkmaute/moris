/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Helmholtz_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Helmholtz_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        IWG_Helmholtz_Bulk::IWG_Helmholtz_Bulk()
        {
            //FIXME set the Helmholtz filter parameter
            mFilterParam = 1.0;

            //            // set the residual dof type
            //            mResidualDofType = { MSI::Dof_Type::VX };
            //
            //            // set the active dof type
            //            mLeaderDofTypes = {{ MSI::Dof_Type::VX }};
        }

        //------------------------------------------------------------------------------
        void IWG_Helmholtz_Bulk::compute_residual( real tWStar )
        {
            //FIXME set unfiltered velocity values at nodes
            Matrix< DDRMat > tVHat  = mNodalWeakBCs;

            // set field interpolator
            Field_Interpolator* vN = mLeaderFI( 0 );

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );

            // compute the residual
            mSet->get_residual()( 0 )( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
                             += ( mFilterParam * trans( vN->dnNdxn( 1 ) ) * vN->gradx( 1 )
                                     + trans( vN->N() ) * ( vN->val() - vN->N() * tVHat ) ) * tWStar;
        }

        //------------------------------------------------------------------------------
        void IWG_Helmholtz_Bulk::compute_jacobian( real tWStar )
        {
            // set field interpolator
            Field_Interpolator* vN = mLeaderFI( 0 );

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );

            // compute the jacobian
            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                    { mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) } )
                    += ( mFilterParam * trans( vN->dnNdxn( 1 ) ) * vN->dnNdxn( 1 )
                            + trans( vN->N() ) * vN->N() ) * tWStar;
        }

        //------------------------------------------------------------------------------
        void IWG_Helmholtz_Bulk::compute_jacobian_and_residual(
                Vector< Vector< Matrix< DDRMat > > > & aJacobians,
                Vector< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, " IWG_Helmholtz_Bulk::compute_jacobian_and_residual - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

