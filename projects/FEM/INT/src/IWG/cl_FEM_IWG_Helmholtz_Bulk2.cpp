/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Helmholtz_Bulk2.cpp
 *
 */

#include "cl_FEM_IWG_Helmholtz_Bulk2.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Helmholtz_Bulk2::IWG_Helmholtz_Bulk2()
        {
            //FIXME set the Helmholtz filter parameter
            mFilterParam = 1.0;

            //FIXME set the sharpness parameter for smoothed Heaviside function
            mSharpParam = 1.0;

            //FIXME set the element size
            mHe = 1.0;

            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::VX };

            // set the active dof type
            mLeaderDofTypes = { { MSI::Dof_Type::VX } ,
                                { MSI::Dof_Type::LS1 } };
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk2::compute_residual( real tWStar )
        {
            // set field interpolator
            Field_Interpolator* vN  = mLeaderFI( 0 );
            Field_Interpolator* phi = mLeaderFI( 1 );

            //FIXME set unfiltered velocity value
            uint tVNBases = vN->get_number_of_space_time_bases();
            uint tVNFields = vN->get_number_of_fields();
            Matrix< DDRMat > aVHat( tVNBases, tVNFields, 1.0 );

            // compute norm( phi )
            real tNormPhi = norm( phi->gradx( 1 ) );
            if( tNormPhi < 1.0e-12 )
            {
                tNormPhi = 1.0e-12;
            }

            // compute Dirac filter
            real tDiracFilter = 0.0;
            if ( phi->val()( 0 ) < 3.0 * mHe )
            {
                real tTanh   = std::tanh( mSharpParam * phi->val()( 0 ) / tNormPhi );
                tDiracFilter = 0.5 * mSharpParam * ( 1.0 - std::pow( tTanh, 2.0 ) );
            }

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );

            // compute the residual
            mSet->get_residual()( 0 )( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
                         += ( mFilterParam * trans( vN->dnNdxn( 1 ) ) * vN->gradx( 1 )
                           + trans( vN->N() ) * ( vN->val() - vN->N() * aVHat ) * tDiracFilter  ) * tWStar;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk2::compute_jacobian( real tWStar )
        {
            // set field interpolator
            Field_Interpolator* vN  = mLeaderFI( 0 );
            Field_Interpolator* phi = mLeaderFI( 1 );

            //FIXME set unfiltered velocity value
            uint tVNBases = vN->get_number_of_space_time_bases();
            uint tVNFields = vN->get_number_of_fields();
            Matrix< DDRMat > aVHat( tVNBases, tVNFields, 1.0 );

            // compute norm( phi ) and derivative wrt phiHat
            real tNormPhi                     = norm( phi->gradx( 1 ) );
            Matrix< DDRMat > tDNormPhiDPhiHat = trans( phi->dnNdxn( 1 ) ) * phi->gradx( 1 ) / tNormPhi;

            // If all values of level set in this element are the same,
            // then gradient is zero, protect from going to NAN/inf
            if( tNormPhi < 1.0e-12 )
            {
                tNormPhi = 1.0e-12;
                uint tNPhiCoeff = phi->get_number_of_space_time_coefficients();
                tDNormPhiDPhiHat.set_size( tNPhiCoeff, 1, 0.0 );
            }

            // compute Dirac filter and derivative wrt phi / norm( phi )
            real tDiracFilter = 0.0;
//            real tDDiracFilter = 0.0;
            if ( phi->val()( 0 ) < 3.0 * mHe )
            {
                real tTanh = std::tanh( mSharpParam * phi->val()( 0 ) / tNormPhi );
                tDiracFilter  = 0.5 * mSharpParam * ( 1.0 - std::pow( tTanh, 2.0 ) );
//                tDDiracFilter = - 2 * mSharpParam * tDiracFilter * tTanh;
            }

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            // compute the jacobian j_vN_vN
            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                  { mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) } )
                    += mFilterParam * trans( vN->dnNdxn( 1 ) ) * vN->dnNdxn( 1 ) + trans( vN->N() ) * vN->N() * tDiracFilter * tWStar;

            // compute the jacobian j_vN_phi //FIXME
//            aJacobians( 0 )( 1 ) = trans( vN->N() ) * ( vN->val() - vN->N() * aVHat )
//                                 * tDDiracFilter * ( phi->N() * tNormPhi  - phi->val()( 0 ) * trans( tDNormPhiDPhiHat ) ) / std::pow( tNormPhi, 2 ) ;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk2::compute_jacobian_and_residual( Vector< Vector< Matrix< DDRMat > > > & aJacobians,
                                                                 Vector< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, " IWG_Helmholtz_Bulk2::compute_jacobian_and_residual - Not implemented.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

