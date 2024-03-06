//
// Created by frank on 12/14/23.
//

#ifndef MORIS_CL_FEM_MODEL_INITIALIZER_HPP
#define MORIS_CL_FEM_MODEL_INITIALIZER_HPP

#include "cl_Vector.hpp"
#include "cl_Param_List.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_Material_Model.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_IQI.hpp"
#include "cl_Library_IO.hpp"
#include "cl_FEM_Set_User_Info.hpp"
#include <string>

namespace moris::fem
{
    constexpr uint FEM_PARAMETER_PROP_INDEX        = 0;
    constexpr uint FEM_PARAMETER_CM_INDEX          = 1;
    constexpr uint FEM_PARAMETER_SP_INDEX          = 2;
    constexpr uint FEM_PARAMETER_IWG_INDEX         = 3;
    constexpr uint FEM_PARAMETER_IQI_INDEX         = 4;
    constexpr uint FEM_PARAMETER_COMPUTATION_INDEX = 5;
    constexpr uint FEM_PARAMETER_FIELD_INDEX       = 6;
    constexpr uint FEM_PARAMETER_PHASE_INDEX       = 7;
    constexpr uint FEM_PARAMETER_MM_INDEX          = 8;

    class Model_Initializer
    {

      public:
        Model_Initializer(
                Vector< Vector< ParameterList > >                aParameterList,
                mtk::Mesh_Pair const                            *aMeshPair,
                std::shared_ptr< Library_IO >                    aLibrary,
                uint                                             aSpatialDimension,
                bool                                             aUseNewGhostSets,
                std::unordered_map< MSI::Dof_Type, moris_index > aDofTypeToBsplineMeshIndex )
                : mParameterList( aParameterList )
                , mMeshPair( aMeshPair )
                , mLibrary( aLibrary )
                , mSpatialDimension( aSpatialDimension )
                , mUseNewGhostSets( aUseNewGhostSets )
                , mDofTypeToBsplineMeshIndex( aDofTypeToBsplineMeshIndex ){};

        virtual void initialize();

        virtual ~Model_Initializer() = default;

        using SetSpecification = std::tuple< std::string, bool, bool >;    // MeshSetName, TimeContinuity, TimeBoundary
        template< typename T >
        using NameToObjectMap = std::map< std::string, std::shared_ptr< T > >;

        const std::map< SetSpecification, fem::Set_User_Info > &get_set_info() const { return mSetInfo; }
        NameToObjectMap< moris::fem::Field > const             &get_fields() const { return mFields; }
        std::map< uint, std::string > const                    &get_field_type_to_name_map() const { return mFieldTypeToName; }
        NameToObjectMap< moris::fem::IQI > const               &get_iqis() const { return mIQIs; }

      protected:
        void create_properties();

        void create_fields();

        virtual void create_material_models() = 0;

        virtual void create_constitutive_models() = 0;

        virtual void create_stabilization_parameters() = 0;

        virtual void create_iwgs() = 0;

        virtual void create_iqis() = 0;

        virtual void create_set_info() = 0;

        virtual void print_physics_model();

        // utility functions

        template< typename T >
        Vector< Vector< T > >
        parameter_to_vec_of_vec( ParameterList const &aParameterList, std::string const &aPropertyName, map< std::string, T > const &aTypeMap ) const;

        std::pair< Vector< Vector< moris::MSI::Dof_Type > >, Vector< std::string > >
        read_dof_dependencies( ParameterList const &aParameterList, mtk::Leader_Follower const aLeaderFollower ) const
        {
            return read_mapped_dependencies( aParameterList, "dof_dependencies", mMSIDofTypeMap, aLeaderFollower );
        };

        std::pair< Vector< Vector< gen::PDV_Type > >, Vector< std::string > >
        read_dv_dependencies( ParameterList const &aParameterList, mtk::Leader_Follower const aLeaderFollower ) const
        {
            return read_mapped_dependencies( aParameterList, "dv_dependencies", mMSIDvTypeMap, aLeaderFollower );
        };

        Vector< std::pair< std::shared_ptr< Property >, std::string > >
        read_properties( ParameterList const &aParameterList, mtk::Leader_Follower const aLeaderFollower ) const;

        std::string
        get_leader_follower_key( std::string aKey, mtk::Leader_Follower const aLeaderFollower ) const;

        Vector< fem::PropertyFunc >
        load_library_property_functions( ParameterList const &aParameterList, std::string const &aPropertyName );

        /**
         * @brief check the ghost set names and conrrect them to automatically include
         * the correct B-spline mesh index when the new ghost sets are used
         * @param aMeshSetName
         * @param aDofType
         */
        void
        check_and_set_ghost_set_names( std::string &aMeshSetName, enum MSI::Dof_Type aDofType );


        // typedefs
        typedef void ( *FEM_Function )(
                moris::Matrix< moris::DDRMat >           &aPropMatrix,
                Vector< moris::Matrix< moris::DDRMat > > &aParameters,
                moris::fem::Field_Interpolator_Manager   *aFIManager );

        // data
        Vector< Vector< ParameterList > >                mParameterList;
        mtk::Mesh_Pair const                            *mMeshPair;
        std::shared_ptr< Library_IO >                    mLibrary;
        uint                                             mSpatialDimension;
        bool                                             mUseNewGhostSets;
        std::unordered_map< MSI::Dof_Type, moris_index > mDofTypeToBsplineMeshIndex;

        std::map< SetSpecification, fem::Set_User_Info > mSetInfo;
        moris::map< std::string, MSI::Dof_Type >         mMSIDofTypeMap   = moris::MSI::get_msi_dof_type_map();
        moris::map< std::string, gen::PDV_Type >         mMSIDvTypeMap    = gen::get_pdv_type_map();
        moris::map< std::string, mtk::Field_Type >       mMTKFieldTypeMap = mtk::get_field_type_map();
        std::map< uint, std::string >                    mFieldTypeToName;

        NameToObjectMap< fem::Property >                mProperties;
        NameToObjectMap< fem::Field >                   mFields;
        NameToObjectMap< fem::Material_Model >          mMaterialModels;
        NameToObjectMap< fem::Constitutive_Model >      mConstitutiveModels;
        NameToObjectMap< fem::Stabilization_Parameter > mStabilizationParameters;
        NameToObjectMap< fem::IWG >                     mIWGs;
        NameToObjectMap< fem::IQI >                     mIQIs;

      private:
        template< typename T >
        std::pair< Vector< Vector< T > >, Vector< std::string > >
        read_mapped_dependencies( ParameterList const &aParameterList, std::string aKey, map< std::string, T > const &aTypeMap, mtk::Leader_Follower aLeaderFollower ) const;
    };

    template< typename T >
    Vector< Vector< T > > Model_Initializer::parameter_to_vec_of_vec( ParameterList const &aParameterList, std::string const &aPropertyName, map< std::string, T > const &aTypeMap ) const
    {
        return string_to_cell_of_cell< T >( aParameterList.get< std::string >( aPropertyName ), aTypeMap );
    }

    template< typename T >
    std::pair< Vector< Vector< T > >, Vector< std::string > >
    Model_Initializer::read_mapped_dependencies(
            ParameterList const         &aParameterList,
            std::string                  aKey,
            map< std::string, T > const &aTypeMap,
            mtk::Leader_Follower const   aLeaderFollower ) const
    {
        std::string const tLeaderFollowerKey     = get_leader_follower_key( aKey, aLeaderFollower );
        auto const &[ tTypeString, tNameString ] = aParameterList.get< std::pair< std::string, std::string > >( tLeaderFollowerKey );
        Vector< Vector< T > > const tTypes       = string_to_cell_of_cell( tTypeString, aTypeMap );
        Vector< std::string > const tTypeNames   = string_to_cell< std::string >( tNameString );
        return { tTypes, tTypeNames };
    }

    template< typename T >
    void ensure_existing_parameter( std::map< std::string, std::shared_ptr< T > > const &aMap, std::string const &aName, std::string const &aContext )
    {
        bool const tExists = aMap.find( aName ) != aMap.end();
        if ( !tExists )
        {
            std::string tPossibleValues;
            for ( auto const &tPair : aMap )
            {
                tPossibleValues += tPair.first + ", ";
            }
            if ( tPossibleValues.size() > 0 )
            {
                tPossibleValues = tPossibleValues.substr( 0, tPossibleValues.size() - 2 );
            }
            MORIS_ERROR( tExists,
                    "Could not create FEM Model with given input: The value '%s' does not exist in the map of %s.\n"
                    "Possible values are: %s\n",
                    aName.c_str(),
                    aContext.c_str(),
                    tPossibleValues.c_str() );
        }
    }
}    // namespace moris::fem


#endif    // MORIS_CL_FEM_MODEL_INITIALIZER_HPP
