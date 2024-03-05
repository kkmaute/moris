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
#include <functional>
#include <map>
#include <memory>
#include <cl_MSI_Dof_Type_Enums.hpp>
#include <cl_Matrix.hpp>
#include <optional>
#include "GEN_Data_Types.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Material_Model.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_Matrix_Arma_Dynamic.hpp"
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


        [[nodiscard]] std::map< std::string, std::shared_ptr< Property > > const                &get_properties() const { return mProperties; }
        void                                                                                     set_property( std::shared_ptr< Property > aProperty, std::string aPropertyString );
        [[nodiscard]] std::map< std::string, std::shared_ptr< fem::Constitutive_Model > > const &get_constitutive_models() const { return mConstitutiveModels; }
        void                                                                                     set_constitutive_models( std::shared_ptr< fem::Constitutive_Model > const &aConstitutiveModel, std::string aConstitutiveModelstring );
        [[nodiscard]] std::map< std::string, std::shared_ptr< fem::Material_Model > > const     &get_material_models() const { return mMaterialModel; }
        void                                                                                     set_material_model( std::shared_ptr< fem::Material_Model > const &aMaterialModel, std::string aMaterialModelString );

        [[nodiscard]] Field_Interpolator_Manager *get_field_interpolator_manager() const { return mFIManager; }
        void                                      set_field_interpolator_manager( Field_Interpolator_Manager *aFIManager );
        [[nodiscard]] Field_Interpolator_Manager *get_eigen_vector_field_interpolator_manager() const { return mEigenFIManager; }
        void                                      set_eigen_fi_manager( Field_Interpolator_Manager *aFIManager ) { mEigenFIManager = aFIManager; }
        [[nodiscard]] Field_Interpolator_Manager *get_previous_fi_manager() const { return mPreviousFIManager; }
        void                                      set_previous_fi_manager( Field_Interpolator_Manager *aPreviousFIManager ) { mPreviousFIManager = aPreviousFIManager; }

        void set_fem_set( fem::Set *aSet );

        void reset_eval_flags();

      private:
        Vector< Vector< MSI::Dof_Type > > build_requested_dof_types( bool aIsStaggered, std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const &aStabilizationParameter );

        template< typename T >
        TypeList< T > build_global_types( std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > aStabilizationParameter, TypeList< T > const &aLocalTypes );

        template< typename T >
        uint get_number_of_unique_types() const { return 0; }
        template<>
        uint get_number_of_unique_types< MSI::Dof_Type >() const { return mSet->get_num_unique_dof_types(); }
        template<>
        uint get_number_of_unique_types< gen::PDV_Type >() const { return mSet->get_num_unique_dv_types(); }
        template<>
        uint get_number_of_unique_types< mtk::Field_Type >() const { return mSet->get_num_unique_field_types(); }

        template< typename T >
        sint get_type_index_from_set( T aType );
        template<>
        sint get_type_index_from_set< MSI::Dof_Type >( MSI::Dof_Type aType ) { return mSet->get_index_from_unique_dof_type_map( aType ); }
        template<>
        sint get_type_index_from_set< gen::PDV_Type >( gen::PDV_Type aType ) { return mSet->get_index_from_unique_dv_type_map( aType ); }
        template<>
        sint get_type_index_from_set< mtk::Field_Type >( mtk::Field_Type aType ) { return mSet->get_index_from_unique_field_type_map( aType ); }

        template< typename P, typename T >
        const TypeList< T > &get_types_from_spec( const std::shared_ptr< P > &aSpec ) const { return TypeList< T >{ {} }; }
        template<>
        const TypeList< MSI::Dof_Type > &get_types_from_spec( const std::shared_ptr< Property > &aProperty ) const { return aProperty->get_dof_type_list(); }
        template<>
        const TypeList< gen::PDV_Type > &get_types_from_spec( const std::shared_ptr< Property > &aProperty ) const { return aProperty->get_dv_type_list(); }
        template<>
        const TypeList< mtk::Field_Type > &get_types_from_spec( const std::shared_ptr< Property > &aProperty ) const { return aProperty->get_field_type_list(); }
        template<>
        const TypeList< MSI::Dof_Type > &get_types_from_spec( const std::shared_ptr< Constitutive_Model > &aConstitutiveModel ) const { return aConstitutiveModel->get_global_dof_type_list(); }
        template<>
        const TypeList< gen::PDV_Type > &get_types_from_spec( const std::shared_ptr< Constitutive_Model > &aConstitutiveModel ) const { return aConstitutiveModel->get_global_dv_type_list(); }
        template<>
        const TypeList< mtk::Field_Type > &get_types_from_spec( const std::shared_ptr< Constitutive_Model > &aConstitutiveModel ) const { return aConstitutiveModel->get_global_field_type_list(); }
        template<>
        const TypeList< MSI::Dof_Type > &get_types_from_spec( const std::shared_ptr< Material_Model > &aMaterialModel ) const { return aMaterialModel->get_global_dof_type_list(); }
        template<>
        const TypeList< MSI::Dof_Type > &get_types_from_spec( const std::shared_ptr< fem::Stabilization_Parameter > &aStabilizationParameter ) const { return aStabilizationParameter->get_global_dof_type_list( mLeaderFollower ); }
        template<>
        const TypeList< gen::PDV_Type > &get_types_from_spec( const std::shared_ptr< fem::Stabilization_Parameter > &aStabilizationParameter ) const { return aStabilizationParameter->get_global_dv_type_list( mLeaderFollower ); }

        template< typename P, typename T >
        void populate_global_types_from_spec( std::map< std::string, std::shared_ptr< P > > const &aMap, TypeList< T > &aGlobalTypes, Matrix< DDSMat > aCheckList );

        using NonUniqueTypes = std::tuple< Vector< MSI::Dof_Type >, Vector< gen::PDV_Type >, Vector< mtk::Field_Type > >;
        template< typename P >
        NonUniqueTypes get_non_unique_types_from_spec( const std::shared_ptr< P > &aSpec ) const { return NonUniqueTypes{ {}, {}, {} }; }
        template<>
        NonUniqueTypes get_non_unique_types_from_spec( const std::shared_ptr< Property > &aProperty ) const
        {
            Vector< MSI::Dof_Type >   tActiveDofTypes;
            Vector< gen::PDV_Type >   tActiveDvTypes;
            Vector< mtk::Field_Type > tActiveFieldTypes;
            aProperty->get_non_unique_dof_dv_and_field_types( tActiveDofTypes, tActiveDvTypes, tActiveFieldTypes );
            return NonUniqueTypes{ tActiveDofTypes, tActiveDvTypes, tActiveFieldTypes };
        }
        template<>
        NonUniqueTypes get_non_unique_types_from_spec( const std::shared_ptr< Constitutive_Model > &aConstitutiveModel ) const
        {
            Vector< MSI::Dof_Type >   tActiveDofTypes;
            Vector< gen::PDV_Type >   tActiveDvTypes;
            Vector< mtk::Field_Type > tActiveFieldTypes;
            aConstitutiveModel->get_non_unique_dof_dv_and_field_types( tActiveDofTypes, tActiveDvTypes, tActiveFieldTypes );
            return NonUniqueTypes{ tActiveDofTypes, tActiveDvTypes, tActiveFieldTypes };
        }
        template<>
        NonUniqueTypes get_non_unique_types_from_spec( const std::shared_ptr< Material_Model > &aMaterialModel ) const
        {
            Vector< MSI::Dof_Type > tActiveDofTypes;
            aMaterialModel->get_non_unique_dof_types( tActiveDofTypes );
            return NonUniqueTypes{ tActiveDofTypes, {}, {} };
        }
        template<>
        NonUniqueTypes get_non_unique_types_from_spec( const std::shared_ptr< fem::Stabilization_Parameter > &aStabilizationParameter ) const
        {
            Vector< MSI::Dof_Type > tActiveDofTypes;
            Vector< gen::PDV_Type > tActiveDvTypes;
            aStabilizationParameter->get_non_unique_dof_and_dv_types( tActiveDofTypes, tActiveDvTypes );
            return NonUniqueTypes{ tActiveDofTypes, tActiveDvTypes, {} };
        };

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

        std::map< std::string, std::shared_ptr< Property > >                mProperties;
        std::map< std::string, std::shared_ptr< fem::Constitutive_Model > > mConstitutiveModels;
        std::map< std::string, std::shared_ptr< fem::Material_Model > >     mMaterialModel;

        Field_Interpolator_Manager *mFIManager         = nullptr;
        Field_Interpolator_Manager *mEigenFIManager    = nullptr;
        Field_Interpolator_Manager *mPreviousFIManager = nullptr;
    };


    template< typename T >
    sint EvaluableSideInformation::get_type_index_from_set( T aType ) { return 0; }


    template< typename T >
    EvaluableSideInformation::TypeList< T > EvaluableSideInformation::build_global_types(
            std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > aStabilizationParameter,
            TypeList< T > const                                                     &aLocalTypes )
    {
        uint const       tNumberOfUniqueTypes = get_number_of_unique_types< T >();
        TypeList< T >    tGlobalTypes( tNumberOfUniqueTypes );
        Matrix< DDSMat > tCheckList( tNumberOfUniqueTypes, 1, -1 );    // (used to avoid repeating a dof or a dv type)

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
                const TypeList< T > &tActiveTypes = get_types_from_spec< P, T >( tSpec );
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
