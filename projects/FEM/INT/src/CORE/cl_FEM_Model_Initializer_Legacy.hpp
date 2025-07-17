/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_FEM_Model_Initializer_Legacy.hpp
 *
 */

#ifndef MORIS_CL_FEM_MODEL_INITIALIZER_LEGACY_HPP
#define MORIS_CL_FEM_MODEL_INITIALIZER_LEGACY_HPP

#include <utility>

#include "cl_Library_IO.hpp"
#include "cl_FEM_Model_Initializer.hpp"
namespace moris::fem
{
    class Model_Initializer_Legacy : public Model_Initializer
    {
      public:
        Model_Initializer_Legacy(
                const Module_Parameter_Lists        &aParameterList,
                std::shared_ptr< Library_IO >                    aLibrary,
                mtk::Mesh_Pair const                            *aMeshPair,
                uint                                             aSpatialDimension,
                bool                                             aUseNewGhostSets,
                std::unordered_map< MSI::Dof_Type, moris_index > aDofTypeToBsplineMeshIndex )
                : Model_Initializer(
                          aParameterList,
                          aMeshPair,
                          std::move( aLibrary ),
                          aSpatialDimension,
                          aUseNewGhostSets,
                          std::move( aDofTypeToBsplineMeshIndex ) ) {};

        ~Model_Initializer_Legacy() override = default;

      private:
        void create_fem_set_info_from_iwgs( std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex );
        void create_fem_set_info_from_iqis( std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex );

        void set_iwg_dof_dependencies( Parameter_List const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower aLeaderFollowerType ) const;
        void set_iwg_dv_dependencies( Parameter_List const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType ) const;
        void set_iwg_field_types( Parameter_List const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType ) const;
        void set_iwg_properties( Parameter_List const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType );
        void set_iwg_material_models( Parameter_List const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType );
        void set_iwg_constitutive_models( Parameter_List const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType );
        void set_iwg_stabilization_parameters( Parameter_List const &aIWGParameter, std::shared_ptr< IWG > &aIWG );
        void set_iwg_residual_dof_type( Parameter_List const &aIWGParameter, std::shared_ptr< IWG > &aIWG ) const;

      protected:
        void create_material_models() override;
        void create_constitutive_models() override;
        void create_stabilization_parameters() override;
        void create_iwgs() override;
        void create_iqis() override;
        void create_set_info() override;
    };
}    // namespace moris::fem

#endif    // MORIS_CL_FEM_MODEL_INITIALIZER_LEGACY_HPP
