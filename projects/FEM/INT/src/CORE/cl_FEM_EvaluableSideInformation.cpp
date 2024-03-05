//
// Created by frank on 3/4/24.
//

#include "cl_FEM_EvaluableSideInformation.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_Vector.hpp"
#include "moris_typedefs.hpp"
#include <map>

namespace moris::fem
{
    void EvaluableSideInformation::set_field_interpolator_manager( Field_Interpolator_Manager *aFIManager )
    {
        auto tUpdateFieldInterpolatorManager = [ aFIManager ]( auto &tItem ) { tItem.second->set_field_interpolator_manager( aFIManager ); };
        mFIManager                           = aFIManager;
        std::for_each( mConstitutiveModels.begin(), mConstitutiveModels.end(), tUpdateFieldInterpolatorManager );
        std::for_each( mMaterialModel.begin(), mMaterialModel.end(), tUpdateFieldInterpolatorManager );
        std::for_each( mProperties.begin(), mProperties.end(), tUpdateFieldInterpolatorManager );
    }

    void EvaluableSideInformation::set_fem_set( fem::Set *aSet )
    {
        mSet                 = aSet;
        auto tUpdatePointers = [ aSet ]( auto &tItem ) { tItem.second->set_set_pointer( aSet ); };
        std::for_each( mConstitutiveModels.begin(), mConstitutiveModels.end(), tUpdatePointers );
        std::for_each( mMaterialModel.begin(), mMaterialModel.end(), tUpdatePointers );
        std::for_each( mProperties.begin(), mProperties.end(), tUpdatePointers );
    }

    void EvaluableSideInformation::set_property( std::shared_ptr< Property > aProperty, std::string aPropertyName )
    {
        mProperties[ aPropertyName ] = aProperty;
    }

    void EvaluableSideInformation::set_material_model( std::shared_ptr< fem::Material_Model > const &aMaterialModel, std::string aMaterialModelString )
    {
        mMaterialModel[ aMaterialModelString ] = aMaterialModel;
    }

    void EvaluableSideInformation::set_constitutive_models( std::shared_ptr< fem::Constitutive_Model > const &aConstitutiveModel, std::string aConstitutiveModelstring )
    {
        mConstitutiveModels[ aConstitutiveModelstring ] = aConstitutiveModel;
    }

    Vector< Vector< MSI::Dof_Type > > const &EvaluableSideInformation::get_global_dof_types(
            std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter )
    {
        if ( !mGlobalDofTypes.has_value() )
        {
            mGlobalDofTypes = this->build_global_types< MSI::Dof_Type >( aStabilizationParameter, mDofTypes );
        }
        return mGlobalDofTypes.value();
    }

    Vector< Vector< gen::PDV_Type > > const &EvaluableSideInformation::get_global_dv_types(
            std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter )
    {
        if ( !mGlobalDvTypes.has_value() )
        {
            mGlobalDvTypes = this->build_global_types< gen::PDV_Type >( aStabilizationParameter, mDvTypes );
        }
        return mGlobalDvTypes.value();
    }
    Vector< Vector< mtk::Field_Type > > const &EvaluableSideInformation::get_global_field_types(
            std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter )
    {
        if ( !mGlobalFieldTypes.has_value() )
        {
            mGlobalFieldTypes = this->build_global_types< mtk::Field_Type >( aStabilizationParameter, mFieldTypes );
        }
        return mGlobalFieldTypes.value();
    }

    Vector< Vector< MSI::Dof_Type > > const &EvaluableSideInformation::get_requested_global_dof_types(
            bool                                                                            aIsStaggered,
            std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter )
    {
        if ( !mRequestedGlobalDofTypes.has_value() )
        {
            mRequestedGlobalDofTypes = this->build_requested_dof_types( aIsStaggered, aStabilizationParameter );
        }
        return mRequestedGlobalDofTypes.value();
    }

    Vector< Vector< MSI::Dof_Type > > EvaluableSideInformation::build_requested_dof_types(
            bool                                                                            aIsStaggered,
            std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter )
    {
        Vector< Vector< MSI::Dof_Type > >       tRequestedGlobalDofTypes;
        Vector< Vector< MSI::Dof_Type > > const tGlobalDofTypes = this->get_global_dof_types( aStabilizationParameter );

        Vector< MSI::Dof_Type > tRequestedDofTypes;
        if ( aIsStaggered )
        {
            tRequestedDofTypes = mSet->get_secondary_dof_types();    // if residual evaluation
        }
        else
        {
            tRequestedDofTypes = mSet->get_requested_dof_types();    // if jacobian evaluation or default (not staggered)
        }

        for ( auto tDofTypes : tRequestedDofTypes )
        {
            // loop over the IWG leader dof types groups
            for ( uint Ik = 0; Ik < tGlobalDofTypes.size(); Ik++ )
            {
                // if requested dof type matches IWG leader dof type
                if ( tGlobalDofTypes( Ik )( 0 ) == tDofTypes )
                {
                    // add the IWG leader dof type to the requested dof list
                    tRequestedGlobalDofTypes.push_back( tGlobalDofTypes( Ik ) );
                    break;
                }
            }
        }
        return tRequestedGlobalDofTypes;
    }

    void EvaluableSideInformation::get_non_unique_dof_dv_and_field_types(
            std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter,
            Vector< MSI::Dof_Type >                                                        &aDofTypes,
            Vector< gen::PDV_Type >                                                        &aDvTypes,
            Vector< mtk::Field_Type >                                                      &aFieldTypes )
    {
        // append the direct dof, dv, and field types
        std::for_each( mDofTypes.begin(), mDofTypes.end(), [ &aDofTypes ]( auto &tDofType ) { aDofTypes.append( tDofType ); } );
        std::for_each( mDvTypes.begin(), mDvTypes.end(), [ &aDvTypes ]( auto &tDvType ) { aDvTypes.append( tDvType ); } );
        std::for_each( mFieldTypes.begin(), mFieldTypes.end(), [ &aFieldTypes ]( auto &tFieldType ) { aFieldTypes.append( tFieldType ); } );

        // append the dof, dv, and field types from the specs (properties, material models, constitutive models and stabilization parameters)
        populate_non_unique_types_from_spec< Property >( mProperties, aDofTypes, aDvTypes, aFieldTypes );
        populate_non_unique_types_from_spec< Material_Model >( mMaterialModel, aDofTypes, aDvTypes, aFieldTypes );
        populate_non_unique_types_from_spec< Constitutive_Model >( mConstitutiveModels, aDofTypes, aDvTypes, aFieldTypes );
        populate_non_unique_types_from_spec< Stabilization_Parameter >( aStabilizationParameter, aDofTypes, aDvTypes, aFieldTypes );
    }

    void EvaluableSideInformation::reset_eval_flags()
    {
        auto tResetEvalFlags = []( auto &tItem ) { tItem.second->reset_eval_flags(); };
        std::for_each( mProperties.begin(), mProperties.end(), tResetEvalFlags );
        std::for_each( mConstitutiveModels.begin(), mConstitutiveModels.end(), tResetEvalFlags );
        std::for_each( mMaterialModel.begin(), mMaterialModel.end(), tResetEvalFlags );
    }


}    // namespace moris::fem