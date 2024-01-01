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

        Vector< fem::Set_User_Info > const                   &get_set_info() const { return mSetInfo; }
        Vector< std::shared_ptr< moris::fem::Field > > const &get_fields() const { return mFields; }
        Vector< moris::sint > const                          &get_field_types() const { return mFieldTypes; }
        Vector< std::shared_ptr< moris::fem::IQI > > const   &get_iqis() const { return mIQIs; }

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
        property_to_vec_of_vec( ParameterList const &aParameterList, std::string const &aPropertyName, map< std::string, T > const &aTypeMap ) const;

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
        Vector< fem::Set_User_Info >                     mSetInfo;

        moris::map< std::string, MSI::Dof_Type >   mMSIDofTypeMap = moris::MSI::get_msi_dof_type_map();
        moris::map< std::string, PDV_Type >        mMSIDvTypeMap  = get_pdv_type_map();
        moris::map< std::string, mtk::Field_Type > mFieldTypeMap  = mtk::get_field_type_map();
        Vector< moris::sint >                      mFieldTypes;

        template< typename T >
        using PointerCell = Vector< std::shared_ptr< T > >;
        PointerCell< fem::Property >                mProperties;
        PointerCell< fem::Field >                   mFields;
        PointerCell< fem::Material_Model >          mMaterialModels;
        PointerCell< fem::Constitutive_Model >      mConstitutiveModels;
        PointerCell< fem::Stabilization_Parameter > mStabilizationParameters;
        PointerCell< fem::IWG >                     mIWGs;
        PointerCell< fem::IQI >                     mIQIs;

        using StringIndexMap = std::map< std::string, uint >;
        StringIndexMap mPropertyMap;
        StringIndexMap mFieldMap;
        StringIndexMap mMaterialModelMap;
        StringIndexMap mConstitutiveModelMap;
        StringIndexMap mStabilizationParameterMap;
        StringIndexMap mIWGMap;
        StringIndexMap mIQIMap;
    };
}    // namespace moris::fem


#endif    // MORIS_CL_FEM_MODEL_INITIALIZER_HPP
