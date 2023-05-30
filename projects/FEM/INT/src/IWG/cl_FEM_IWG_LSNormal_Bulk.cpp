/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_LSNormal_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_LSNormal_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        IWG_LSNormal_Bulk::IWG_LSNormal_Bulk()
        {
//            // set the residual dof type
//            mResidualDofType = { MSI::Dof_Type::NLSX,
//                                 MSI::Dof_Type::NLSY,
//                                 MSI::Dof_Type::NLSZ };
//
//            // set the active dof type
//            mLeaderDofTypes = {{ MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY, MSI::Dof_Type::NLSZ },
//                               { MSI::Dof_Type::LS1 }};
        }

//------------------------------------------------------------------------------
        void IWG_LSNormal_Bulk::compute_residual( real tWStar )
        {
            // check leader field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // set field interpolators
            Field_Interpolator* nPhi = mLeaderFI( 0 );
            Field_Interpolator* phi  = mLeaderFI( 1 );

            // build the global shape functions matrix for vectorial field nPhi
            uint tNBasesNPhi  = nPhi->get_number_of_space_time_bases();
            uint tNFieldsNPhi = nPhi->get_number_of_fields();
            Matrix< DDRMat > tNNPhi( tNFieldsNPhi, tNFieldsNPhi * tNBasesNPhi, 0.0 );
            for( uint i = 0; i < tNFieldsNPhi; i++ )
            {
                tNNPhi({i,i},{i * tNBasesNPhi, (i+1) * tNBasesNPhi - 1}) = nPhi->N().get_row( 0 );
            }

            // compute norm( phi )
            real tNormPhi = norm( phi->gradx( 1 ) );
            if( tNormPhi < 1.0e-12 )
            {
                tNormPhi = 1.0e-12;
            }

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );

            // compute residual
            mSet->get_residual()( 0 )( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
                    += trans( tNNPhi ) * ( trans( nPhi->val() ) - phi->gradx( 1 ) / tNormPhi ) * tWStar;
        }

//------------------------------------------------------------------------------
        void IWG_LSNormal_Bulk::compute_jacobian( real tWStar )
        {
            // check leader field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // set field interpolators
            Field_Interpolator* nPhi = mLeaderFI( 0 );
            Field_Interpolator* phi  = mLeaderFI( 1 );

            // build the global shape functions matrix for vectorial field nPhi
            uint tNBasesNPhi  = nPhi->get_number_of_space_time_bases();
            uint tNFieldsNPhi = nPhi->get_number_of_fields();
            Matrix< DDRMat > tNNPhi( tNFieldsNPhi, tNFieldsNPhi * tNBasesNPhi, 0.0 );
            for( uint i = 0; i < tNFieldsNPhi; i++ )
            {
                tNNPhi({i,i},{i * tNBasesNPhi, (i+1) * tNBasesNPhi - 1}) = nPhi->N().get_row( 0 );
            }

            // compute norm( phi ) and derivative wrt phiHat
            real tNormPhi                     = norm( phi->gradx( 1 ) );
            Matrix< DDRMat > tDNormPhiDPhiHat = trans( phi->dnNdxn( 1 ) ) * phi->gradx( 1 ) / tNormPhi;

            // If all values of level set in this element are the same,
            // then gradient is zero, protect from going to NAN/inf
            if( tNormPhi < 1.0e-12 )
            {
                tNormPhi = 1.0e-12;
                uint tNPhiBases = phi->get_number_of_space_time_bases();
                tDNormPhiDPhiHat.set_size( tNPhiBases, 1, 0.0 );
            }

            // set the jacobian size
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );

            // compute j_nPhi_nPhi
            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                  { mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) } )
                    += trans( tNNPhi ) * tNNPhi * tWStar;

            // compute j_nPhi_phi
//            aJacobians( 0 )( 1 ) = - trans( tNNPhi ) * ( phi->dnNdxn( 1 ) * tNormPhi - phi->gradx( 1 ) * trans( tDNormPhiDPhiHat ) ) / std::pow( tNormPhi, 2 ) * tWStar;
        }

//------------------------------------------------------------------------------
        void IWG_LSNormal_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                               moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
//            // check leader field interpolators
//            this->check_dof_field_interpolators();
//            this->check_dv_field_interpolators();
//
//            // set field interpolators
//            Field_Interpolator* nPhi = mLeaderFI( 0 );
//            Field_Interpolator* phi  = mLeaderFI( 1 );
//
//            // build the global shape functions matrix for vectorial field nPhi
//            uint tNBasesNPhi  = nPhi->get_number_of_space_time_bases();
//            uint tNFieldsNPhi = nPhi->get_number_of_fields();
//            Matrix< DDRMat > tNNPhi( tNFieldsNPhi, tNFieldsNPhi * tNBasesNPhi, 0.0 );
//            for( uint i = 0; i < tNFieldsNPhi; i++ )
//            {
//                tNNPhi({i,i},{i * tNBasesNPhi, (i+1) * tNBasesNPhi - 1}) = nPhi->N().get_row( 0 );
//            }
//
//            // compute norm( phi ) and derivative wrt phiHat
//            real tNormPhi                     = norm( phi->gradx( 1 ) );
//            Matrix< DDRMat > tDNormPhiDPhiHat = trans( phi->gradx( 1 ) ) * phi->dnNdxn( 1 ) / tNormPhi;
//
//            // If all values of level set in this element are the same,
//            // then gradient is zero, protect from going to NAN/inf
//            if( tNormPhi < 1.0e-12 )
//            {
//                tNormPhi = 1.0e-12;
//                uint tNPhiBases = phi->get_number_of_space_time_bases();
//                tDNormPhiDPhiHat.set_size( tNPhiBases, 1, 0.0 );
//            }
//
//            //set resiaul size
//            this->set_residual( aResidual );
//
//            // compute residual
//            aResidual( 0 ) = trans( tNNPhi ) * ( trans( nPhi->val() ) - phi->gradx( 1 ) / tNormPhi );
//
//            // set jacobain size
//            this->set_jacobian( aJacobians );
//
//            // compute j_nPhi_nPhi
//            aJacobians( 0 )( 0 ) = trans( tNNPhi ) * tNNPhi;
//
//            // compute j_nPhi_phi
//            aJacobians( 0 )( 1 ) = - trans( tNNPhi ) * ( phi->dnNdxn( 1 ) * tNormPhi - phi->gradx( 1 ) * trans( tDNormPhiDPhiHat ) ) / std::pow( tNormPhi, 2 ) ;

        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

