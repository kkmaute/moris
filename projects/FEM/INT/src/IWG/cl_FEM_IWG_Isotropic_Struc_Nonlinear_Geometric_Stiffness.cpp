/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Nonlinear_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Nonlinear_Geometric_Stiffness.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_CM_Struc_Nonlinear_Isotropic.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"
#include "cl_FEM_Set.hpp"

// #include "fn_trans.hpp"
// #include "fn_norm.hpp"
// #include "fn_eye.hpp"
// #include "fn_cond.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Nonlinear_Geometric_Stiffness::IWG_Isotropic_Struc_Nonlinear_Geometric_Stiffness(
                enum CM_Function_Type aStressType,
                enum CM_Function_Type aStrainType )
        {
            mIWGType = moris::fem::IWG_Type::STRUC_NON_LINEAR_GEOMETRIC_STIFFNESS;

            // assign stress and strain type to evaluate the IWG
            mStressType = aStressType;
            mStrainType = aStrainType;

            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Nonlinear_Geometric_Stiffness::compute_residual( real aWStar )
        {
            MORIS_ERROR( false,
                    "IWG_Isotropic_Struc_Nonlinear_Geometric_Stiffness::compute_residual - should not be called" );
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Nonlinear_Geometric_Stiffness::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get elasticity CM and cast into
            const std::shared_ptr< Constitutive_Model >& tCMElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // cast constitutive model base class pointer to either linear or nonlinear elasticity model
            CM_Struc_Nonlinear_Isotropic* tCMNonLinarElasticityPtr = nullptr;
            CM_Struc_Linear*              tCMLinarElasticityPtr    = nullptr;
            if ( dynamic_cast< CM_Struc_Nonlinear_Isotropic* >( tCMElasticity.get() ) )
            {
                tCMNonLinarElasticityPtr = dynamic_cast< CM_Struc_Nonlinear_Isotropic* >( tCMElasticity.get() );
            }
            else
            {
                tCMLinarElasticityPtr = dynamic_cast< CM_Struc_Linear* >( tCMElasticity.get() );
            }

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            MORIS_ASSERT( !( tCMElasticity->get_plane_type() == Model_Type::AXISYMMETRIC
                                  and tPropThickness == nullptr ),
                    "IWG_Isotropic_Struc_Nonlinear_Bulk::compute_jacobian - must define axis of rotation "
                    "using IWG \"Thickness\" property as {{x1,y1},{x2,y2}} if using axisymmetric formulation" );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get the number of leader dof dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over the leader dof dependencies
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if constitutive model depends on the dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // compute the contribution to Jacobian (note negative sign, used to build eigen problem Ke x = -v Kg x)
                    if ( tCMNonLinarElasticityPtr )
                    {
                        tJac -= aWStar * tCMNonLinarElasticityPtr->GeometricStiffness( tDofType, mStressType );
                    }
                    else
                    {
                        tJac -= aWStar * tCMLinarElasticityPtr->GeometricStiffness( tDofType );
                    }
                }
            }
            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Nonlinear_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Nonlinear_Geometric_Stiffness::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Bulk::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Nonlinear_Geometric_Stiffness::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Bulk::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
