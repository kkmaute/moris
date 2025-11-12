/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Nonlinear_Bedding.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IWG_Nonlinear_Bedding::IWG_Nonlinear_Bedding()
    {
        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Bedding" ]   = static_cast< uint >( IWG_Property_Type::BEDDING );
        mPropertyMap[ "Bedding_Threshold" ] = static_cast< uint >( IWG_Property_Type::BEDDING_THRESHOLD );

    }

    //------------------------------------------------------------------------------

    void IWG_Nonlinear_Bedding::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for dof type
            Field_Interpolator* tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // get bedding property
            const std::shared_ptr< Property >& tPropBedding =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // get bedding threshold property
            const std::shared_ptr< Property >& tPropBeddingThreshold =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING_THRESHOLD ) );

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // get displacement vector
            Matrix< DDRMat > tDisp = tDisplacementFI->val();

            //get displacement transpose times displacement
            Matrix< DDRMat > tDispTDisp = trans( tDisp ) * tDisp;

            // Get basis function matrix
            Matrix< DDRMat > tN = tDisplacementFI->N();

            //get bedding and threshold value
            real tBeddingValue = tPropBedding->val()( 0 );
            real tBeddingThresholdValue = tPropBeddingThreshold->val()( 0 );

            // compute the residual
            tRes += aWStar * 0.5 * tBeddingValue * (( tDispTDisp( 0 ) ) * (2.0 / ( 1.0 + std::pow((tDispTDisp( 0 ) / tBeddingThresholdValue ), 2)) ) * ( trans( tDisp ) * tN ) +
                    2.0 * std::tanh( tDispTDisp( 0 ) / tBeddingThresholdValue ) * ( trans( tDisp ) * tN ) );
            

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Nonlinear_Bedding::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for given dof type
            Field_Interpolator * tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get bedding property
            const std::shared_ptr< Property > & tPropBedding =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // get bedding threshold property
            const std::shared_ptr< Property > & tPropBeddingThreshold =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING_THRESHOLD ) );

            // get the number of leader dof dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over the leader dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // get displacement vector
                Matrix< DDRMat > tDisp = tDisplacementFI->val();

                //get displacement transpose times displacement
                Matrix< DDRMat > tDispTDisp = trans( tDisp ) * tDisp;

                // Get basis function matrix
                Matrix< DDRMat > tN = tDisplacementFI->N();

                //get bedding and threshold value
                real tBeddingValue = tPropBedding->val()( 0 );
                real tBeddingThresholdValue = tPropBeddingThreshold->val()( 0 );

                // Compute the jacobian
                tJac += aWStar * 0.5 * tBeddingValue * ( 4.0 * tDispTDisp( 0 ) * ( trans( tN ) * ( tN ) ) * std::pow( ( 1.0 + std::pow( ( tDispTDisp( 0 )/ tBeddingThresholdValue ) , 2 ) ), 2 ) *  ( ( 1.0 +  std::pow( ( tDispTDisp( 0 )/ tBeddingThresholdValue ) , 2)  ) +
                                                           tDispTDisp( 0 ) * ( 2.0 * tDispTDisp( 0 )/ tBeddingThresholdValue )) + 2.0 * tDispTDisp( 0 ) * ( trans( tN ) * tN ) * ( 1.0 / ( 1.0 + std::pow( (tDispTDisp( 0 )/tBeddingThresholdValue ), 2 ))) +
                                                           2.0 * std::tanh( tDispTDisp( 0 )/ tBeddingThresholdValue ) * ( trans( tN ) * tN ) + 4.0 * tDispTDisp( 0 ) * ( trans( tN ) * tN ) * ( 1.0 / ( 1.0 + std::pow( ( tDispTDisp( 0 )/tBeddingThresholdValue ), 2 ) ) ) );

            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Nonlinear_Bedding::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Nonlinear_Bedding::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Nonlinear_Bedding::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Nonlinear_Bedding::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
