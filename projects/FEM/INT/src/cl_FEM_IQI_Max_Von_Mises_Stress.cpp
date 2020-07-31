/*
 * cl_FEM_IQI_Max_Von_Mises_Stress.cpp
 *
 *  Created on: Jul 30, 2020
 *      Author: wunsch
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Max_Von_Mises_Stress.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IQI_Max_Von_Mises_Stress::IQI_Max_Von_Mises_Stress()
        {
            // set IQI type
            mIQIType = vis::Output_Type::MAX_STRESS;

            // set FEM IQI type
            mFEMIQIType = fem::IQI_Type::MAX_VON_MISES_STRESS;

            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "ReferenceValue" ]    = IQI_Property_Type::REFERENCE_VALUE;
            mPropertyMap[ "Exponent" ]          = IQI_Property_Type::EXPONENT;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IQI_Constitutive_Type::ELAST_LIN_ISO;
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Von_Mises_Stress::set_property( std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster)
        {
            // can only be master
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IQI_Max_Von_Mises_Stress::set_property - can only be master." );

            // check that aPropertyString makes sense
            if ( mPropertyMap.find( aPropertyString ) == mPropertyMap.end() )
            {
                std::string tErrMsg =
                        std::string( "CM_Diffusion_Linear_Isotropic_Phase_Change::set_property - Unknown aPropertyString : ") +
                        aPropertyString;

                MORIS_ERROR( false , tErrMsg.c_str() );
            }

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Von_Mises_Stress::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster)
        {
            // can only be master
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IQI_Max_Von_Mises_Stress::set_constitutive model - can only be master." );

            // check that aPropertyString makes sense
            if ( mConstitutiveMap.find( aConstitutiveString ) == mConstitutiveMap.end() )
            {
                std::string tErrMsg =
                        std::string( "IQI_Max_Von_Mises_Stress::set_constitutive_model - Unknown aConstitutiveString : ") +
                        aConstitutiveString;

                MORIS_ERROR( false , tErrMsg.c_str() );
            }

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Von_Mises_Stress::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get indices for properties, CM and SP
            uint tElastLinIsoIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );

            // get stress
            Matrix< DDRMat > tStressVector = mMasterCM( tElastLinIsoIndex )->flux();

            // initialize stress contributions to von mises stress
            real tShearStressContribution = 0.0;
            real tNormalStressContribution = 0.0;

            // pull apart stress vector into components
            uint tNumStressComponents = tStressVector.length();
            switch  (tNumStressComponents)
            {
                // 2D plane stress
                case 3:
                    tNormalStressContribution = std::pow( tStressVector( 0 ) - tStressVector( 1 ), 2.0 );
                    tShearStressContribution = std::pow( tStressVector( 2 ), 2.0 );
                    break;

                // 2D plane strain
                case 4:
                    tNormalStressContribution =
                            std::pow( tStressVector( 0 ) - tStressVector( 1 ), 2.0 ) +
                            std::pow( tStressVector( 1 ) - tStressVector( 2 ), 2.0 ) +
                            std::pow( tStressVector( 2 ) - tStressVector( 0 ), 2.0 );
                    tShearStressContribution = std::pow( tStressVector( 3 ), 2.0 );
                    break;

                // 3D
                case 6:
                    tNormalStressContribution =
                            std::pow( tStressVector( 0 ) - tStressVector( 1 ), 2.0 ) +
                            std::pow( tStressVector( 1 ) - tStressVector( 2 ), 2.0 ) +
                            std::pow( tStressVector( 2 ) - tStressVector( 0 ), 2.0 );
                    tShearStressContribution =
                            std::pow( tStressVector( 3 ), 2.0 ) +
                            std::pow( tStressVector( 4 ), 2.0 ) +
                            std::pow( tStressVector( 5 ), 2.0 );
                    break;

                // Unknown size - error
                default:
                    MORIS_ERROR( false, "IQI_Max_Von_Mises_Stress::compute_QI - Stress vector of unknown size; 3, 4 or 6 components expected." );
            }

            // compute von mises stress value
            real tVMStress = std::sqrt( 0.5 * tNormalStressContribution + 3.0 * tShearStressContribution );

            // check if properties set
            MORIS_ERROR(mMasterProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) ) != nullptr,
                    "IQI_Max_Dof - no reference value set");

            MORIS_ERROR(mMasterProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) ) != nullptr,
                    "IQI_Max_Dof - no exponent set");

            // get property values
            real tRefValue = mMasterProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) )->val()( 0 );
            real tExponent = mMasterProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) )->val()( 0 );

            // check if dof index was set (for the case of vector field)
            if( mMasterDofTypes( 0 ).size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Max_Dof::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            aQI = {{ std::pow( 1/tRefValue * tVMStress - 1.0, tExponent ) }};
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Von_Mises_Stress::compute_QI( moris::real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get indices for properties, CM and SP
            uint tElastLinIsoIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );

            // get stress
            Matrix< DDRMat > tStressVector = mMasterCM( tElastLinIsoIndex )->flux();

            // initialize stress contributions to von mises stress
            real tShearStressContribution = 0.0;
            real tNormalStressContribution = 0.0;

            // pull apart stress vector into components
            uint tNumStressComponents = tStressVector.length();
            switch  (tNumStressComponents)
            {
                // 2D plane stress
                case 3:
                    tNormalStressContribution = std::pow( tStressVector( 0 ) - tStressVector( 1 ), 2.0 );
                    tShearStressContribution = std::pow( tStressVector( 2 ), 2.0 );
                    break;

                // 2D plane strain
                case 4:
                    tNormalStressContribution =
                            std::pow( tStressVector( 0 ) - tStressVector( 1 ), 2.0 ) +
                            std::pow( tStressVector( 1 ) - tStressVector( 2 ), 2.0 ) +
                            std::pow( tStressVector( 2 ) - tStressVector( 0 ), 2.0 );
                    tShearStressContribution = std::pow( tStressVector( 3 ), 2.0 );
                    break;

                // 3D
                case 6:
                    tNormalStressContribution =
                            std::pow( tStressVector( 0 ) - tStressVector( 1 ), 2.0 ) +
                            std::pow( tStressVector( 1 ) - tStressVector( 2 ), 2.0 ) +
                            std::pow( tStressVector( 2 ) - tStressVector( 0 ), 2.0 );
                    tShearStressContribution =
                            std::pow( tStressVector( 3 ), 2.0 ) +
                            std::pow( tStressVector( 4 ), 2.0 ) +
                            std::pow( tStressVector( 5 ), 2.0 );
                    break;

                // Unknown size - error
                default:
                    MORIS_ERROR( false, "IQI_Max_Von_Mises_Stress::compute_QI - Stress vector of unknown size; 3, 4 or 6 components expected." );
            }

            // compute von mises stress value
            real tVMStress = std::sqrt( 0.5 * tNormalStressContribution + 3.0 * tShearStressContribution );

            // check if properties are set
            MORIS_ERROR(mMasterProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) ) != nullptr,
                    "IQI_Max_Dof - no reference value set");

            MORIS_ERROR(mMasterProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) ) != nullptr,
                    "IQI_Max_Dof - no exponent set");

            // get property values
            real tRefValue = mMasterProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) )->val()( 0 );
            real tExponent = mMasterProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) )->val()( 0 );

            // check if dof index was set (for the case of vector field)
            if( mMasterDofTypes( 0 ).size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Max_Dof::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            mSet->get_QI()( tQIIndex ).matrix_data() += { aWStar *  std::pow( 1/tRefValue * tVMStress - 1.0, tExponent ) };
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Von_Mises_Stress::compute_dQIdDof( Matrix< DDRMat > & adQIdDof )
        {
            MORIS_ERROR(0,"IQI_Max_Von_Mises_Stress::compute_dQIdDof - Derivatives of stress wrt dof not implemented");
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



