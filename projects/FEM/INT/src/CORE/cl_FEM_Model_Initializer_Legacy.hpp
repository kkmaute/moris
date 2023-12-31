//
// Created by frank on 12/14/23.
//

#ifndef MORIS_CL_FEM_MODEL_INITIALIZER_LEGACY_HPP
#define MORIS_CL_FEM_MODEL_INITIALIZER_LEGACY_HPP

#include "cl_Cell.hpp"
#include "cl_Param_List.hpp"
#include "cl_Library_IO.hpp"
#include "cl_FEM_Model_Initializer.hpp"
namespace moris::fem
{
    class Model_Initializer_Legacy : public Model_Initializer
    {
      public:
        Model_Initializer_Legacy(
                moris::Cell< moris::Cell< ParameterList > >      aParameterList,
                std::shared_ptr< Library_IO >                    aLibrary,
                mtk::Mesh_Pair const*                            aMeshPair,
                uint                                             aSpatialDimension,
                bool                                             aUseNewGhostSets,
                std::unordered_map< MSI::Dof_Type, moris_index > aDofTypeToBsplineMeshIndex )
                : Model_Initializer( aParameterList, aMeshPair, aLibrary, aSpatialDimension, aUseNewGhostSets, aDofTypeToBsplineMeshIndex ){};

        virtual ~Model_Initializer_Legacy() = default;

      private:
        void create_fem_set_info_from_iwgs( bool aIsAnalyticalSA, FDScheme_Type const &aFDSchemeForSA, real const aFDPerturbation, std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex );
        void create_fem_set_info_from_iqis( bool const aIsAnalyticalSA, FDScheme_Type const &aFDSchemeForSA, real const aFDPerturbation, std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex );

        void set_iwg_dof_dependencies( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower aLeaderFollowerType ) const;
        void set_iwg_dv_dependencies( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType ) const;
        void set_iwg_field_types( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType ) const;
        void set_iwg_properties( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType );
        void set_iwg_material_models( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType );
        void set_iwg_constitutive_models( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType );
        void set_iwg_stabilization_parameters( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG );
        void set_iwg_residual_dof_type( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG ) const;

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
