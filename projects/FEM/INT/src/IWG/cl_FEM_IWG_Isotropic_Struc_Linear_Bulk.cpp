/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Linear_Bulk::IWG_Isotropic_Struc_Linear_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Load" ]      = static_cast< uint >( IWG_Property_Type::LOAD );
            mPropertyMap[ "Bedding" ]   = static_cast< uint >( IWG_Property_Type::BEDDING );
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Bulk::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here displacement), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for dof type
            Field_Interpolator* tDisplacementFI =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // get body load property
            const std::shared_ptr< Property > & tPropLoad =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::LOAD ) );

            // get bedding property
            const std::shared_ptr< Property > & tPropBedding =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // get elasticity CM
            const std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            MORIS_ASSERT( !(tCMElasticity->get_plane_type() == Model_Type::AXISYMMETRIC
                    and tPropThickness == nullptr),
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_residual - must define axis of rotation "
                    "using IWG \"Thickness\" property as {{x1,y1},{x2,y2}} if using axisymmetric formulation");

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex } );

            // compute the residual
            tRes += aWStar * ( trans( tCMElasticity->testStrain() ) * tCMElasticity->flux() );

            // if body load
            if ( tPropLoad != nullptr )
            {
                // compute body load contribution
                tRes -= aWStar * ( trans( tDisplacementFI->N() ) * tPropLoad->val() );
            }

            // if bedding
            if ( tPropBedding != nullptr )
            {
                // compute body load contribution
                tRes += aWStar * ( trans( tDisplacementFI->N() ) * tDisplacementFI->val() * tPropBedding->val() );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here displacement), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for given dof type
            Field_Interpolator * tDisplacementFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get body load property
            const std::shared_ptr< Property > & tPropLoad =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::LOAD ) );

            // get bedding property
            const std::shared_ptr< Property > & tPropBedding =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // get elasticity CM
            const std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            MORIS_ASSERT( !(tCMElasticity->get_plane_type() == Model_Type::AXISYMMETRIC
                    and tPropThickness == nullptr),
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian - must define axis of rotation "
                    "using IWG \"Thickness\" property as {{x1,y1},{x2,y2}} if using axisymmetric formulation");

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // get the number of master dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the master dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } );

                // if body load
                if ( tPropLoad != nullptr )
                {
                    // if property depends on the dof type
                    if ( tPropLoad->check_dof_dependency( tDofType ) )
                    {
                        // compute the contribution to Jacobian
                        tJac -= aWStar * ( trans( tDisplacementFI->N() ) * tPropLoad->dPropdDOF( tDofType ) );
                    }
                }

                // if bedding
                if ( tPropBedding != nullptr )
                {
                    if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // if dof type is displacement, add bedding contribution
                        tJac += aWStar * ( trans( tDisplacementFI->N() ) *  tDisplacementFI->N() * tPropBedding->val()(0) );
                    }

                    // consider contributions from dependency of bedding parameter on DOFs
                    if ( tPropBedding->check_dof_dependency( tDofType ) )
                    {
                        tJac += aWStar * (
                                trans( tDisplacementFI->N() ) *  tDisplacementFI->val() * tPropBedding->dPropdDOF( tDofType ) );
                    }
                }

                // if constitutive model depends on the dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // compute the contribution to Jacobian
                    tJac += aWStar * ( trans( tCMElasticity->testStrain() ) * tCMElasticity->dFluxdDOF( tDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Bulk::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Bulk::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

