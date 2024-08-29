/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_L2.cpp
 *
 */

#include "cl_FEM_IWG_User_Defined.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "op_times.hpp"    //LINALG/src
#include "fn_norm.hpp"     //LINALG/src
#include "fn_trans.hpp"    //LINALG/src
#include "fn_dot.hpp"      //LINALG/src

namespace moris::fem
{
    //------------------------------------------------------------------------------
    IWG_User_Defined::IWG_User_Defined()
    {
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "WeakForm" ]  = static_cast< uint >( IWG_Property_Type::IWG_FORMULATION );
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );
    }
    //------------------------------------------------------------------------------

    void
    IWG_User_Defined::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get weak form
            const std::shared_ptr< Property >& tPropWeakform =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::IWG_FORMULATION ) );

            MORIS_ASSERT( tPropWeakform,
                    "IWG_User_Defined::compute_residual - WeakForm property is does not exist" );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // compute residual of weak form in property
            mSet->get_residual()( 0 )( { tResStartIndex, tResStopIndex } ) += aWStar * tPropWeakform->val();
        }

        //------------------------------------------------------------------------------

        void
        IWG_User_Defined::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get weak form
            const std::shared_ptr< Property >& tPropWeakform =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::IWG_FORMULATION ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over leader dof type dependencies
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 1 );

                // compute jacobian of weak form in property
                if ( tPropWeakform->check_dof_dependency( tDofType ) )
                {
                    mSet->get_jacobian()(
                            { tResStartIndex, tResStopIndex },
                            { tDepStartIndex, tDepStopIndex } ) += aWStar * tPropWeakform->dPropdDOF( tDofType );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG_User_Defined::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_L2::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_User_Defined::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_L2::compute_dRdp - not implemented." );
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
