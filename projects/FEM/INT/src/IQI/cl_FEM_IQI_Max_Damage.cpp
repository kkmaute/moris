/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Max_Damage.cpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_DAMAGE_CPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_DAMAGE_CPP_

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Max_Damage.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic_Damage.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
        IQI_Max_Damage::IQI_Max_Damage()
        {
            // set FEM IQI type
            mFEMIQIType = fem::IQI_Type::MAX_DAMAGE;

            // set the property pointer Vector size
            mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "ReferenceValue" ] = static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE );
            mPropertyMap[ "Shift" ]          = static_cast< sint >( IQI_Property_Type::SHIFT );
            mPropertyMap[ "Exponent" ]       = static_cast< uint >( IQI_Property_Type::EXPONENT );
            mPropertyMap[ "Select" ]         = static_cast< uint >( IQI_Property_Type::SELECT );

            // set size for the constitutive model pointer Vector
            mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElasticDamage" ] = static_cast< uint >( IQI_Constitutive_Type::ELASTIC_DAMAGE );
        }

        //------------------------------------------------------------------------------
        void
        IQI_Max_Damage::compute_QI( Matrix< DDRMat >& aQI )
        {
            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            // check if properties set
            MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) ) != nullptr,
                    "IQI_Max_Damage - no reference value set" );

            MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) ) != nullptr,
                    "IQI_Max_Damage - no shift set" );

            MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) ) != nullptr,
                    "IQI_Max_Damage - no exponent set" );

            // get property values
            real tRefValue = mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) )->val()( 0 );
            real tShift    = mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) )->val()( 0 );
            real tExponent = mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) )->val()( 0 );
            real tActivity = 1.0;

            MORIS_ERROR( tRefValue != 0.0,
                    "IQI_Max_Damage - Reference Value set to zero" );

            if ( mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) ) != nullptr )
                tActivity = mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) )->val()( 0 );

            // compute max damage
            aQI = { { tActivity * std::pow( std::max( tCMElasticityDamagePtr->smooth_damage()( 0 ) / tRefValue - tShift, 0.0 ), tExponent ) } };
        }

        //------------------------------------------------------------------------------

        void
        IQI_Max_Damage::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // create matrix storage for IQI
            Matrix<DDRMat> tQI( 1, 1 );

            // compute QI
            this->compute_QI( tQI );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI( 0 );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Max_Damage::compute_dQIdu( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            // check if properties set
            MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) ) != nullptr,
                    "IQI_Max_Damage - no reference value set" );

            MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) ) != nullptr,
                    "IQI_Max_Damage - no shift set" );

            MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) ) != nullptr,
                    "IQI_Max_Damage - no exponent set" );

            // get property values
            real tRefValue = mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) )->val()( 0 );
            real tShift    = mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) )->val()( 0 );
            real tExponent = mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) )->val()( 0 );
            real tActivity = 1.0;

            MORIS_ERROR( tRefValue != 0.0,
                    "IQI_Max_Damage - Reference Value set to zero");

            if ( mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) ) != nullptr )
                tActivity = mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) )->val()( 0 );

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // compute dQIdu for indirect dof dependencies
            for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDof );

                // get leader index for residual dof type, indices for assembly
                uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
                uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

                if ( tCMElasticityDamage->check_dof_dependency( tDofType ))
                {
                    // compute max damage derivative
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex }, { 0, 0 } ) +=
                            tActivity * aWStar * ( tExponent * trans( tCMElasticityDamagePtr->dSmoothDamagedu( tDofType ) ) / tRefValue ) *    //
                            std::pow( std::max( tCMElasticityDamagePtr->smooth_damage()( 0 ) / tRefValue - tShift, 0.0 ), tExponent - 1.0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI_Max_Damage::compute_dQIdu( moris::Vector< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            // check if properties set
            MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) ) != nullptr,
                    "IQI_Max_Damage - no reference value set" );

            MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) ) != nullptr,
                    "IQI_Max_Damage - no shift set" );

            MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) ) != nullptr,
                    "IQI_Max_Damage - no exponent set" );

            // get property values
            real tRefValue = mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) )->val()( 0 );
            real tShift    = mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) )->val()( 0 );
            real tExponent = mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) )->val()( 0 );
            real tActivity = 1.0;

            if ( mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) ) != nullptr )
                tActivity = mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) )->val()( 0 );

            MORIS_ERROR( tRefValue != 0.0,
                    "IQI_Max_Damage - Reference Value set to zero" );

            if ( tCMElasticityDamage->check_dof_dependency( aDofType ) )
            {
                // compute max damage derivative
                adQIdu = { { tActivity * ( tExponent * trans( tCMElasticityDamagePtr->dSmoothDamagedu( aDofType ) ) / tRefValue ) *    //
                             std::pow( std::max( tCMElasticityDamagePtr->smooth_damage()( 0 ) / tRefValue - tShift, 0.0 ), tExponent - 1.0 ) } };
            }
        }

        //------------------------------------------------------------------------------

    }    // namespace fem
}    // namespace moris
#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_DAMAGE_CPP_ */
