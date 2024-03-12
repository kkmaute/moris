//
// Created by frank on 3/4/24.
//

#include "cl_FEM_EvaluableSideInformation.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Vector.hpp"
#include "moris_typedefs.hpp"
#include <algorithm>
#include <map>

namespace moris::fem
{
    void EvaluableSideInformation::set_fi_manager( Field_Interpolator_Manager *aFIManager )
    {
        auto tUpdateFieldInterpolatorManager = [ aFIManager, this ]( auto &tItem ) {
            if ( tItem.second != nullptr )
            {
                tItem.second->set_field_interpolator_manager( aFIManager );
            }
        };
        mFIManager = aFIManager;
        std::for_each( mConstitutiveModels.begin(), mConstitutiveModels.end(), tUpdateFieldInterpolatorManager );
        std::for_each( mMaterialModel.begin(), mMaterialModel.end(), tUpdateFieldInterpolatorManager );
        std::for_each( mProperties.begin(), mProperties.end(), tUpdateFieldInterpolatorManager );
    }

    void EvaluableSideInformation::set_fem_set( fem::Set *aSet )
    {
        mSet                 = aSet;
        auto tUpdatePointers = [ aSet ]( auto &tItem ) {
            if ( tItem.second != nullptr )
            {
                tItem.second->set_set_pointer( aSet );
            }
        };
        std::for_each( mConstitutiveModels.begin(), mConstitutiveModels.end(), tUpdatePointers );
        std::for_each( mMaterialModel.begin(), mMaterialModel.end(), tUpdatePointers );
        std::for_each( mProperties.begin(), mProperties.end(), tUpdatePointers );
    }

    void EvaluableSideInformation::set_property( std::shared_ptr< Property > aProperty, std::string aPropertyType )
    {
        ensure_valid_option( mProperties, aPropertyType, "property" );
        mProperties[ aPropertyType ] = aProperty;
    }

    void EvaluableSideInformation::set_material_model( std::shared_ptr< fem::Material_Model > const &aMaterialModel, std::string aMaterialModelType )
    {
        ensure_valid_option( mMaterialModel, aMaterialModelType, "material_model" );
        mMaterialModel[ aMaterialModelType ] = aMaterialModel;
    }

    void EvaluableSideInformation::set_constitutive_model( std::shared_ptr< fem::Constitutive_Model > const &aConstitutiveModel, std::string aConstitutiveModelType )
    {
        ensure_valid_option( mConstitutiveModels, aConstitutiveModelType, "constitutive_model" );
        mConstitutiveModels[ aConstitutiveModelType ] = aConstitutiveModel;
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

    Vector< Vector< MSI::Dof_Type > > const &EvaluableSideInformation::get_requested_global_dof_types( std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter )
    {
        if ( !mRequestedGlobalDofTypes.has_value() )
        {
            mRequestedGlobalDofTypes = this->build_requested_dof_types( aStabilizationParameter );
        }
        return mRequestedGlobalDofTypes.value();
    }

    Vector< Vector< MSI::Dof_Type > > EvaluableSideInformation::build_requested_dof_types(
            std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter )
    {
        Vector< Vector< MSI::Dof_Type > >       tRequestedGlobalDofTypes;
        Vector< Vector< MSI::Dof_Type > > const tGlobalDofTypes = this->get_global_dof_types( aStabilizationParameter );

        Vector< MSI::Dof_Type > tRequestedDofTypes;
        if ( mIsStaggered )
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
            Vector< mtk::Field_Type >                                                      &aFieldTypes ) const
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
        auto tResetEvalFlags = []( auto &tItem ) {
            if ( tItem.second != nullptr )
            {
                tItem.second->reset_eval_flags();
            }
        };
        std::for_each( mProperties.begin(), mProperties.end(), tResetEvalFlags );
        std::for_each( mConstitutiveModels.begin(), mConstitutiveModels.end(), tResetEvalFlags );
        std::for_each( mMaterialModel.begin(), mMaterialModel.end(), tResetEvalFlags );
    }

    void EvaluableSideInformation::set_is_staggered( bool aIsStaggered )
    {
        // invalidate requested global dof types if the staggered flag changes
        if ( mIsStaggered != aIsStaggered )
        {
            mRequestedGlobalDofTypes.reset();
        }
        mIsStaggered = aIsStaggered;
    }

    template<>
    uint EvaluableSideInformation::get_number_of_unique_types< MSI::Dof_Type >() const { return mSet->get_num_unique_dof_types(); }
    template<>
    uint EvaluableSideInformation::get_number_of_unique_types< gen::PDV_Type >() const { return mSet->get_num_unique_dv_types(); }
    template<>
    uint EvaluableSideInformation::get_number_of_unique_types< mtk::Field_Type >() const { return mSet->get_num_unique_field_types(); }


    template< typename T >
    sint EvaluableSideInformation::get_type_index_from_set( T aType ) { return 0; }

    template<>
    sint EvaluableSideInformation::get_type_index_from_set< MSI::Dof_Type >( MSI::Dof_Type aType ) { return mSet->get_index_from_unique_dof_type_map( aType ); }
    template<>
    sint EvaluableSideInformation::get_type_index_from_set< gen::PDV_Type >( gen::PDV_Type aType ) { return mSet->get_index_from_unique_dv_type_map( aType ); }
    template<>
    sint EvaluableSideInformation::get_type_index_from_set< mtk::Field_Type >( mtk::Field_Type aType ) { return mSet->get_index_from_unique_field_type_map( aType ); }


    template< typename P, typename T >
    Vector< Vector< T > > const
    EvaluableSideInformation::get_types_from_spec( const std::shared_ptr< P > &aSpec ) const { return Vector< Vector< T > >{}; }
    template<>
    Vector< Vector< MSI::Dof_Type > > const
    EvaluableSideInformation::get_types_from_spec( const std::shared_ptr< Property > &aProperty ) const { return aProperty->get_dof_type_list(); }
    template<>
    Vector< Vector< gen::PDV_Type > > const
    EvaluableSideInformation::get_types_from_spec( const std::shared_ptr< Property > &aProperty ) const { return aProperty->get_dv_type_list(); }
    template<>
    Vector< Vector< mtk::Field_Type > > const
    EvaluableSideInformation::get_types_from_spec( const std::shared_ptr< Property > &aProperty ) const { return aProperty->get_field_type_list(); }
    template<>
    Vector< Vector< MSI::Dof_Type > > const
    EvaluableSideInformation::get_types_from_spec( const std::shared_ptr< Constitutive_Model > &aConstitutiveModel ) const { return aConstitutiveModel->get_global_dof_type_list(); }
    template<>
    Vector< Vector< gen::PDV_Type > > const
    EvaluableSideInformation::get_types_from_spec( const std::shared_ptr< Constitutive_Model > &aConstitutiveModel ) const { return aConstitutiveModel->get_global_dv_type_list(); }
    template<>
    Vector< Vector< mtk::Field_Type > > const
    EvaluableSideInformation::get_types_from_spec( const std::shared_ptr< Constitutive_Model > &aConstitutiveModel ) const { return aConstitutiveModel->get_global_field_type_list(); }
    template<>
    Vector< Vector< MSI::Dof_Type > > const
    EvaluableSideInformation::get_types_from_spec( const std::shared_ptr< Material_Model > &aMaterialModel ) const { return aMaterialModel->get_global_dof_type_list(); }
    template<>
    Vector< Vector< MSI::Dof_Type > > const
    EvaluableSideInformation::get_types_from_spec( const std::shared_ptr< fem::Stabilization_Parameter > &aStabilizationParameter ) const { return aStabilizationParameter->get_global_dof_type_list( mLeaderFollower ); }
    template<>
    Vector< Vector< gen::PDV_Type > > const
    EvaluableSideInformation::get_types_from_spec( const std::shared_ptr< fem::Stabilization_Parameter > &aStabilizationParameter ) const { return aStabilizationParameter->get_global_dv_type_list( mLeaderFollower ); }

    template< typename P >
    EvaluableSideInformation::NonUniqueTypes
    EvaluableSideInformation::get_non_unique_types_from_spec( const std::shared_ptr< P > &aSpec ) const { return NonUniqueTypes{ {}, {}, {} }; }

    template<>
    EvaluableSideInformation::NonUniqueTypes
    EvaluableSideInformation::get_non_unique_types_from_spec( const std::shared_ptr< Property > &aProperty ) const
    {
        Vector< MSI::Dof_Type >   tActiveDofTypes;
        Vector< gen::PDV_Type >   tActiveDvTypes;
        Vector< mtk::Field_Type > tActiveFieldTypes;
        aProperty->get_non_unique_dof_dv_and_field_types( tActiveDofTypes, tActiveDvTypes, tActiveFieldTypes );
        return NonUniqueTypes{ tActiveDofTypes, tActiveDvTypes, tActiveFieldTypes };
    }
    template<>
    EvaluableSideInformation::NonUniqueTypes
    EvaluableSideInformation::get_non_unique_types_from_spec( const std::shared_ptr< Constitutive_Model > &aConstitutiveModel ) const
    {
        Vector< MSI::Dof_Type >   tActiveDofTypes;
        Vector< gen::PDV_Type >   tActiveDvTypes;
        Vector< mtk::Field_Type > tActiveFieldTypes;
        aConstitutiveModel->get_non_unique_dof_dv_and_field_types( tActiveDofTypes, tActiveDvTypes, tActiveFieldTypes );
        return NonUniqueTypes{ tActiveDofTypes, tActiveDvTypes, tActiveFieldTypes };
    }
    template<>
    EvaluableSideInformation::NonUniqueTypes
    EvaluableSideInformation::get_non_unique_types_from_spec( const std::shared_ptr< Material_Model > &aMaterialModel ) const
    {
        Vector< MSI::Dof_Type > tActiveDofTypes;
        aMaterialModel->get_non_unique_dof_types( tActiveDofTypes );
        return NonUniqueTypes{ tActiveDofTypes, {}, {} };
    }
    template<>
    EvaluableSideInformation::NonUniqueTypes
    EvaluableSideInformation::get_non_unique_types_from_spec( const std::shared_ptr< fem::Stabilization_Parameter > &aStabilizationParameter ) const
    {
        Vector< MSI::Dof_Type > tActiveDofTypes;
        Vector< gen::PDV_Type > tActiveDvTypes;
        aStabilizationParameter->get_non_unique_dof_and_dv_types( tActiveDofTypes, tActiveDvTypes );
        return NonUniqueTypes{ tActiveDofTypes, tActiveDvTypes, {} };
    };

    template< typename T >
    EvaluableSideInformation::TypeList< T > EvaluableSideInformation::build_global_types(
            std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > aStabilizationParameter,
            TypeList< T > const                                                     &aLocalTypes )
    {
        uint const       tNumberOfUniqueTypes = get_number_of_unique_types< T >();
        Matrix< DDSMat > tCheckList( tNumberOfUniqueTypes, 1, -1 );    // (used to avoid repeating a dof or a dv type)

        TypeList< T > tGlobalTypes;
        for ( uint iType = 0; iType < aLocalTypes.size(); iType++ )
        {
            sint tTypeIndex          = get_type_index_from_set< T >( aLocalTypes( iType )( 0 ) );    // get set index for dof type
            tCheckList( tTypeIndex ) = 1;                                                            // put the dof type in the checklist
            tGlobalTypes.push_back( aLocalTypes( iType ) );                                          // put the dof type in the global type list
        }

        populate_global_types_from_spec< Property, T >( mProperties, tGlobalTypes, tCheckList );
        populate_global_types_from_spec< Constitutive_Model, T >( mConstitutiveModels, tGlobalTypes, tCheckList );
        populate_global_types_from_spec< Material_Model, T >( mMaterialModel, tGlobalTypes, tCheckList );
        populate_global_types_from_spec< fem::Stabilization_Parameter, T >( aStabilizationParameter, tGlobalTypes, tCheckList );

        tGlobalTypes.shrink_to_fit();
        return tGlobalTypes;
    }

    template< typename P, typename T >
    void EvaluableSideInformation::populate_global_types_from_spec( std::map< std::string, std::shared_ptr< P > > const &aMap, TypeList< T > &aGlobalTypes, Matrix< DDSMat > aCheckList )
    {
        for ( auto const &[ _, tSpec ] : aMap )
        {
            if ( tSpec != nullptr )
            {
                // get types (dof, dv, field) from specification (property, constitutive model, material model, stabilization parameter)
                const TypeList< T > tActiveTypes = get_types_from_spec< P, T >( tSpec );
                for ( uint iType = 0; iType < tActiveTypes.size(); iType++ )
                {
                    sint tDofTypeIndex = get_type_index_from_set< T >( tActiveTypes( iType )( 0 ) );    // get set index for type (dof, dv, field)
                    if ( aCheckList( tDofTypeIndex ) != 1 )                                             // if type enum not in the list
                    {
                        aCheckList( tDofTypeIndex ) = 1;                    // put the dof type in the checklist
                        aGlobalTypes.push_back( tActiveTypes( iType ) );    // put the dof type in the global type list
                    }
                }
            }
        }
    }
}    // namespace moris::fem