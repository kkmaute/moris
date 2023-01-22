/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_L2.cpp
 *
 */

#include "cl_FEM_IWG_L2.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "op_times.hpp"    //LINALG/src
#include "fn_norm.hpp"     //LINALG/src
#include "fn_trans.hpp"    //LINALG/src
#include "fn_dot.hpp"      //LINALG/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        IWG_L2::IWG_L2()
        {
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "L2coefficient" ] = static_cast< uint >( IWG_Property_Type::L2COEFFICIENT );
            mPropertyMap[ "H1coefficient" ] = static_cast< uint >( IWG_Property_Type::H1COEFFICIENT );
            mPropertyMap[ "Diffusion" ]     = static_cast< uint >( IWG_Property_Type::DIFFUSION );
            mPropertyMap[ "Source" ]        = static_cast< uint >( IWG_Property_Type::SOURCE );
        }
        //------------------------------------------------------------------------------

        void
        IWG_L2::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get the field interpolator for residual dof type
            Field_Interpolator* tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
       
            // get L2 coefficient
            const std::shared_ptr< Property >& tPropL2Term =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::L2COEFFICIENT ) );

            real tL2Term = 1.0;

            if ( tPropL2Term != nullptr )
            {
                    tL2Term = tPropL2Term->val()(0);
            }

            // get H1 coefficient property
            const std::shared_ptr< Property >& tPropH1Term =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::H1COEFFICIENT ) );

            // get diffusion coefficient property
            const std::shared_ptr< Property >& tPropDiffusion =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIFFUSION ) );

            // get source property
            const std::shared_ptr< Property >& tPropSource =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SOURCE ) );

            // extract residual sub-vector
            auto tRes = mSet->get_residual()( 0 )( { tResStartIndex, tResStopIndex } );

            // add L2 contribution to residual
            if ( tPropSource != nullptr )
            {
                tRes += aWStar * tL2Term * tFI->N_trans() * ( tFI->val() - tPropSource->val() );
            }
            else    // FIXME: replace mNodalWeakBCs with field interpolator
            {
                tRes += aWStar * tL2Term * tFI->N_trans() * ( tFI->val() - tFI->N() * mNodalWeakBCs );
            }

            // add H1 contribution to residual
            if ( tPropH1Term != nullptr )
            {
                if ( tPropSource != nullptr )
                {
                    tRes += aWStar * tPropH1Term->val() * trans( tFI->dnNdxn( 1 ) ) * ( tFI->gradx( 1 ) - tPropSource->dnPropdxn( 1 ) );
                }
                else
                {
                    tRes += aWStar * tPropH1Term->val() * trans( tFI->dnNdxn( 1 ) ) * ( tFI->gradx( 1 ) - tFI->dnNdxn( 1 ) * mNodalWeakBCs );
                }
            }

            // add diffusion term to residual
            if ( tPropDiffusion != nullptr )
            {
                tRes += aWStar * tPropDiffusion->val() * trans( tFI->dnNdxn( 1 ) ) * tFI->gradx( 1 );
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG_L2::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get the field interpolator for residual dof type
            Field_Interpolator* tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

             // get L2 coefficient
            const std::shared_ptr< Property >& tPropL2Term =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::L2COEFFICIENT ) );

            real tL2Term = 1.0;

            if ( tPropL2Term != nullptr )
            {
                    tL2Term = tPropL2Term->val()(0);
            }

            // get H1 coefficient property
            const std::shared_ptr< Property >& tPropH1Term =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::H1COEFFICIENT ) );

            // get diffusion coefficient property
            const std::shared_ptr< Property >& tPropDiffusion =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIFFUSION ) );

            // get source property
            const std::shared_ptr< Property >& tPropSource =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SOURCE ) );

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over master dof type dependencies
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 1 );

                // extract jacobian sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tResStartIndex, tResStopIndex },
                        { tDepStartIndex, tDepStopIndex } );

                // if residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute Jacobian
                    tJac += aWStar * tL2Term * tFI->N_trans() * tFI->N();

                    // add contribution from H1 term to Jacobian
                    if ( tPropH1Term != nullptr )
                    {
                        tJac += aWStar * tPropH1Term->val() * trans( tFI->dnNdxn( 1 ) ) * tFI->dnNdxn( 1 );
                    }

                    // add contribution from diffusion term to Jacobian
                    if ( tPropDiffusion != nullptr )
                    {
                        tJac += aWStar * tPropDiffusion->val() * trans( tFI->dnNdxn( 1 ) ) * tFI->dnNdxn( 1 );
                    }
                }

                if ( tPropSource != nullptr )
                {
                    if ( tPropSource->check_dof_dependency( tDofType ) )
                    {
                        tJac -= aWStar * tL2Term * tFI->N_trans() * tPropSource->dPropdDOF( tDofType );
                    }

                    if ( tPropH1Term != nullptr )
                    {
                        if ( tPropH1Term->check_space_dependency( 1 ) )
                        {
                            MORIS_ERROR( false, "IWG_L2::compute_jacobian - H1 contribution for spatially varying properties not implemented." );
                            // tJac -= aWStar * tPropH1Term->val() * trans( tFI->dnNdxn( 1 ) ) * !! missing functionality of property !!;
                        }
                    }
                }

                // add H1 contribution to residual
                if ( tPropH1Term != nullptr )
                {
                    if ( tPropH1Term->check_dof_dependency( tDofType ) )
                    {
                        MORIS_ERROR( false, "IWG_L2::compute_jacobian - dof dependence of H1 penalty not implemented." );
                    }
                }

                // add diffusion term to residual
                if ( tPropDiffusion != nullptr )
                {
                    if ( tPropDiffusion->check_dof_dependency( tDofType ) )
                    {
                        MORIS_ERROR( false, "IWG_L2::compute_jacobian - dof dependence of diffusion not implemented." );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG_L2::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_L2::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_L2::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_L2::compute_dRdp - not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

