/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Olsson_CLS_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Olsson_CLS_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        IWG_Olsson_CLS_Bulk::IWG_Olsson_CLS_Bulk()
        {
            //FIXME set field upper and lower bound
            mPhiUB = 1.0;
            mPhiLB = 0.0;

            //FIXME set Olsson CLS epsilon parameter
            mEpsilon = 1.0;

            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::LS2 };

            // set the active dof type
            mLeaderDofTypes = {{ MSI::Dof_Type::LS2 },
                               { MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY } };
        }

//------------------------------------------------------------------------------
        void IWG_Olsson_CLS_Bulk::compute_residual( real tWStar )
        {
            // check leader field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // set field interpolators
            Field_Interpolator* phi  = mLeaderFI( 0 );
            Field_Interpolator* nPhi = mLeaderFI( 1 );

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );

            // compute residual
            mSet->get_residual()( 0 )( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
                    +=      ( trans( phi->N() ) * phi->gradt( 1 )
                           - trans( phi->dnNdxn( 1 ) ) * ( ( phi->val()( 0 ) - mPhiLB ) * ( mPhiUB - phi->val()( 0 ) )
                           - mEpsilon  * dot( phi->gradx( 1 ), nPhi->val() ) ) * trans( nPhi->val() ) ) * tWStar;
        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Bulk::compute_jacobian( real tWStar )
        {
            // check leader field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // set field interpolators
            Field_Interpolator* phi  = mLeaderFI( 0 );
            Field_Interpolator* nPhi = mLeaderFI( 1 );

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );

            //compute the jacobians
            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                  { mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) } )
                    += ( trans( phi->N() )  * phi->dnNdtn( 1 )
                                 - trans( phi->dnNdxn( 1 ) ) * ( ( mPhiUB + mPhiLB - 2 * phi->val()( 0 ) ) * trans( nPhi->val() ) * phi->N()
                                                          - mEpsilon * ( trans( nPhi->val() ) * nPhi->val() * phi->dnNdxn( 1 ) ) ) ) * tWStar;

//           // build the global shape functions matrix for vectorial field nPhi //FIXME
//           uint tNBasesNPhi  = nPhi->get_number_of_space_time_bases();
//           uint tNFieldsNPhi = nPhi->get_number_of_fields();
//           Matrix< DDRMat > tNNPhi( tNFieldsNPhi, tNFieldsNPhi * tNBasesNPhi, 0.0 );
//           for( uint i = 0; i < tNFieldsNPhi; i++ )
//           {
//               tNNPhi({i,i},{i * tNBasesNPhi, (i+1) * tNBasesNPhi - 1}) = nPhi->N().get_row( 0 );
//           }
//
//
//            aJacobians( 0 )( 1 ) = ( - trans( phi->dnNdxn( 1 ) ) *(
//                                    ( phi->val()( 0 ) - mPhiLB ) * ( mPhiUB - phi->val()( 0 ) ) * tNNPhi
//                                   - mEpsilon * trans( nPhi->val() ) * trans( phi->gradx( 1 ) ) * tNNPhi
//                                   - mEpsilon * dot( phi->gradx( 1 ), nPhi->val() ) * tNNPhi ) ) * tWStar;

        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Bulk::compute_jacobian_and_residual( Vector< Vector< Matrix< DDRMat > > > & aJacobians,
                                                                 Vector< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Olsson_CLS_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

