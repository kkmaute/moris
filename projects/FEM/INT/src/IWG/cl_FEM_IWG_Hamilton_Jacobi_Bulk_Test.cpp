/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.cpp
 *
 */

#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Hamilton_Jacobi_Bulk_Test::IWG_Hamilton_Jacobi_Bulk_Test(){}

        //------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk_Test::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get index for residual dof type and indices for residual assembly
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator * tFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // FIXME should be a property
            // velocity field value
            Matrix< DDRMat > tVN( 1, tFI->get_number_of_fields(), 1.0 );

            //compute the residual
            mSet->get_residual()( 0 )( { tStartRow, tEndRow }, { 0, 0 } )
                    += trans( tFI->N() ) * ( tFI->gradt( 1 ) + tVN * tFI->gradx( 1 ) ) * aWStar;
        }

        //------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk_Test::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get index for residual dof type and indices for residual assembly
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator * tFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // FIXME should be a property
            // velocity field value
            Matrix< DDRMat > tVN( 1, tFI->get_number_of_fields(), 1.0 );

            // get indices for jacobian assembly
            uint tStartCol = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 );
            uint tEndCol   = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 );

            // compute the jacobian
            mSet->get_jacobian()( { tStartRow, tEndRow }, { tStartCol, tEndCol } )
                    += trans( tFI->N() ) * ( tFI->dnNdtn( 1 ) + tVN * tFI->dnNdxn( 1 ) ) * aWStar;
        }

        //------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk_Test::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Hamilton_Jacobi_Bulk_Test::compute_jacobian_and_residual - Not implemented.");
        }

        //------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk_Test::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Hamilton_Jacobi_Bulk_Test::compute_dRdp - Not implemented." );
        }
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

