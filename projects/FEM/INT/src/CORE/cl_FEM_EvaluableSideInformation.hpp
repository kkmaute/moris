/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_EvaluableSideInformation.hpp
 *
 */
#pragma once
#include <map>
#include <memory>
#include <optional>
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Material_Model.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_Logger.hpp"
#include "cl_MTK_Enums.hpp"
#include <cl_MSI_Dof_Type_Enums.hpp>
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_Matrix_Arma_Dynamic.hpp"
#include <cl_Matrix.hpp>
#include "linalg_typedefs.hpp"

namespace moris::fem
{
    class EvaluableSideInformation
    {
      public:
        EvaluableSideInformation( mtk::Leader_Follower const aLeaderFollower )
                : mLeaderFollower( aLeaderFollower ){};

        template< typename T >
        using TypeList = Vector< Vector< T > >;

        [[nodiscard]] std::string const                 &get_phase_name() const { return mPhaseName; }
        void                                             set_phase_name( std::string const &aPhaseName ) { mPhaseName = aPhaseName; }
        [[nodiscard]] TypeList< MSI::Dof_Type > const   &get_dof_types() const { return mDofTypes; }
        void                                             set_dof_types( TypeList< MSI::Dof_Type > const &aDofTypes ) { mDofTypes = aDofTypes; }
        [[nodiscard]] TypeList< MSI::Dof_Type > const   &get_global_dof_types( std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter );
        [[nodiscard]] TypeList< MSI::Dof_Type > const   &get_requested_global_dof_types( bool aIsStaggered, std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter );
        [[nodiscard]] TypeList< gen::PDV_Type > const   &get_dv_types() const { return mDvTypes; }
        void                                             set_dv_types( TypeList< gen::PDV_Type > const &aDvTypes ) { mDvTypes = aDvTypes; }
        [[nodiscard]] TypeList< gen::PDV_Type > const   &get_global_dv_types( std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter );
        [[nodiscard]] TypeList< mtk::Field_Type > const &get_field_types() const { return mFieldTypes; }
        void                                             set_field_types( TypeList< mtk::Field_Type > const &aFieldTypes ) { mFieldTypes = aFieldTypes; }
        [[nodiscard]] TypeList< mtk::Field_Type > const &get_global_field_types( std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter );

        void get_non_unique_dof_dv_and_field_types(
                std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter,
                Vector< MSI::Dof_Type >                                                        &aDofTypes,
                Vector< gen::PDV_Type >                                                        &aDvTypes,
                Vector< mtk::Field_Type >                                                      &aFieldTypes );

        template< typename EnumType >
        void                                                                      init_property( std::string const &aPropertyName, EnumType const aPropertyType );
        void                                                                      set_property( std::shared_ptr< Property > aProperty, std::string aPropertyType );
        [[nodiscard]] std::map< std::string, std::shared_ptr< Property > > const &get_properties() const { return mProperties; }
        template< typename EnumType >
        std::shared_ptr< Property > const &get_property( EnumType const aPropertyType ) const { return mProperties.at( mPropertyTypeToString.at( static_cast< uint >( aPropertyType ) ) ); }

        template< typename EnumType >
        void                                                                                     init_constitutive_model( std::string const &aConstitutiveModelName, EnumType const aConstitutiveModelType );
        void                                                                                     set_constitutive_model( std::shared_ptr< fem::Constitutive_Model > const &aConstitutiveModel, std::string aConstitutiveModelstring );
        [[nodiscard]] std::map< std::string, std::shared_ptr< fem::Constitutive_Model > > const &get_constitutive_models() const { return mConstitutiveModels; }
        template< typename EnumType >
        std::shared_ptr< fem::Constitutive_Model > const &get_constitutive_model( EnumType const aConstitutiveModelType ) const { return mConstitutiveModels.at( mConstitutiveTypeToString.at( static_cast< uint >( aConstitutiveModelType ) ) ); }

        template< typename EnumType >
        void                                                                                 init_material_model( std::string const &aMaterialModelName, EnumType const aMaterialModelType );
        void                                                                                 set_material_model( std::shared_ptr< fem::Material_Model > const &aMaterialModel, std::string aMaterialModelString );
        [[nodiscard]] std::map< std::string, std::shared_ptr< fem::Material_Model > > const &get_material_models() const { return mMaterialModel; }
        template< typename EnumType >
        std::shared_ptr< fem::Material_Model > const &get_material_model( EnumType const aMaterialModelType ) const { return mMaterialModel.at( mMaterialTypeToString.at( static_cast< uint >( aMaterialModelType ) ) ); }

        [[nodiscard]] Field_Interpolator_Manager *get_fi_manager() const { return mFIManager; }
        void                                      set_fi_manager( Field_Interpolator_Manager *aFIManager );

        [[nodiscard]] Field_Interpolator_Manager *get_eigen_vector_fi_manager() const { return mEigenFIManager; }
        void                                      set_eigen_vector_fi_manager( Field_Interpolator_Manager *aFIManager ) { mEigenFIManager = aFIManager; }

        [[nodiscard]] Field_Interpolator_Manager *get_previous_fi_manager() const { return mPreviousFIManager; }
        void                                      set_previous_fi_manager( Field_Interpolator_Manager *aPreviousFIManager ) { mPreviousFIManager = aPreviousFIManager; }

        void set_fem_set( fem::Set *aSet );

        void reset_eval_flags();

        template< typename T >
        void ensure_valid_option( std::map< std::string, std::shared_ptr< T > > aMap, std::string const &aOption, std::string const &tContext ) const;

      private:
        Vector< Vector< MSI::Dof_Type > > build_requested_dof_types( bool aIsStaggered, std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter );

        template< typename T >
        TypeList< T > build_global_types( std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > aStabilizationParameter, TypeList< T > const &aLocalTypes );

        template< typename T >
        uint get_number_of_unique_types() const { return 0; }

        template< typename T >
        sint get_type_index_from_set( T aType );

        template< typename P, typename T >
        const TypeList< T > get_types_from_spec( const std::shared_ptr< P > &aSpec ) const { return TypeList< T >{}; }

        template< typename P, typename T >
        void populate_global_types_from_spec( std::map< std::string, std::shared_ptr< P > > const &aMap, TypeList< T > &aGlobalTypes, Matrix< DDSMat > aCheckList );

        using NonUniqueTypes = std::tuple< Vector< MSI::Dof_Type >, Vector< gen::PDV_Type >, Vector< mtk::Field_Type > >;
        template< typename P >
        NonUniqueTypes get_non_unique_types_from_spec( const std::shared_ptr< P > &aSpec ) const { return NonUniqueTypes{ {}, {}, {} }; }

        template< typename P >
        void populate_non_unique_types_from_spec(
                std::map< std::string, std::shared_ptr< P > > const &aSpecificationMap,
                Vector< MSI::Dof_Type >                              aDofTypes,
                Vector< gen::PDV_Type >                              aDvTypes,
                Vector< mtk::Field_Type >                            aFieldTypes ) const
        {
            for ( auto const &[ _, tSpec ] : aSpecificationMap )
            {
                if ( tSpec != nullptr )
                {
                    NonUniqueTypes tActiveTypes = get_non_unique_types_from_spec< P >( tSpec );
                    aDofTypes.append( std::get< 0 >( tActiveTypes ) );
                    aDvTypes.append( std::get< 1 >( tActiveTypes ) );
                    aFieldTypes.append( std::get< 2 >( tActiveTypes ) );
                }
            }
        }

      private:
        mtk::Leader_Follower const mLeaderFollower;
        Set                       *mSet = nullptr;

        std::string                 mPhaseName;
        TypeList< MSI::Dof_Type >   mDofTypes;
        TypeList< gen::PDV_Type >   mDvTypes;
        TypeList< mtk::Field_Type > mFieldTypes;

        std::optional< TypeList< MSI::Dof_Type > >   mRequestedGlobalDofTypes;
        std::optional< TypeList< MSI::Dof_Type > >   mGlobalDofTypes;
        std::optional< TypeList< gen::PDV_Type > >   mGlobalDvTypes;
        std::optional< TypeList< mtk::Field_Type > > mGlobalFieldTypes;

        std::map< std::string, std::shared_ptr< Property > > mProperties;
        std::map< uint, std::string >                        mPropertyTypeToString;

        std::map< std::string, std::shared_ptr< fem::Constitutive_Model > > mConstitutiveModels;
        std::map< uint, std::string >                                       mConstitutiveTypeToString;

        std::map< std::string, std::shared_ptr< fem::Material_Model > > mMaterialModel;
        std::map< uint, std::string >                                   mMaterialTypeToString;

        Field_Interpolator_Manager *mFIManager         = nullptr;
        Field_Interpolator_Manager *mEigenFIManager    = nullptr;
        Field_Interpolator_Manager *mPreviousFIManager = nullptr;
    };

    template< typename EnumType >
    void EvaluableSideInformation::init_property( std::string const &aPropertyName, EnumType const aPropertyType )
    {
        mPropertyTypeToString[ static_cast< uint >( aPropertyType ) ] = aPropertyName;
        mProperties[ aPropertyName ]                                  = nullptr;
    }

    template< typename EnumType >
    void EvaluableSideInformation::init_constitutive_model( std::string const &aConstitutiveModelName, EnumType const aConstitutiveModelType )
    {
        mConstitutiveTypeToString[ static_cast< uint >( aConstitutiveModelType ) ] = aConstitutiveModelName;
        mConstitutiveModels[ aConstitutiveModelName ]                              = nullptr;
    }

    template< typename EnumType >
    void EvaluableSideInformation::init_material_model( std::string const &aMaterialModelName, EnumType const aMaterialModelType )
    {
        mMaterialTypeToString[ static_cast< uint >( aMaterialModelType ) ] = aMaterialModelName;
        mMaterialModel[ aMaterialModelName ]                               = nullptr;
    }

    template< typename T >
    void EvaluableSideInformation::ensure_valid_option( std::map< std::string, std::shared_ptr< T > > aMap, std::string const &aOption, std::string const &tContext ) const
    {
        if ( aMap.find( aOption ) == aMap.end() )
        {
            std::string tValidOptions;
            for ( auto const &[ tName, _ ] : aMap )
            {
                tValidOptions += tName + ", ";
            }
            if ( !tValidOptions.empty() )
            {
                tValidOptions.pop_back();
                tValidOptions.pop_back();
            }
            MORIS_LOG_ERROR( "Invalid option '%s' for '%s'. Valid options are:\n    %s", aOption.c_str(), tContext.c_str(), tValidOptions.c_str() );
        }
    }

    template< typename T >
    sint EvaluableSideInformation::get_type_index_from_set( T aType ) { return 0; }


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
