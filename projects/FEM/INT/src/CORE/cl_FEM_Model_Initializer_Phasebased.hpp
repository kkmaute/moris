//
// Created by frank on 12/14/23.
//

#ifndef MORIS_CL_FEM_MODEL_INITIALIZER_PHASEBASED_HPP
#define MORIS_CL_FEM_MODEL_INITIALIZER_PHASEBASED_HPP

#include "cl_FEM_Model_Initializer.hpp"
#include "cl_FEM_Phase_User_Info.hpp"
#include <map>
#include "cl_FEM_Set_User_Info.hpp"
namespace moris::fem
{
    class Model_Initializer_Phasebased : public Model_Initializer
    {
      public:
        Model_Initializer_Phasebased(
                Vector< Vector< ParameterList > >                aParameterList,
                std::shared_ptr< Library_IO >                    aLibrary,
                mtk::Mesh_Pair const                            *aMeshPair,
                uint                                             aSpatialDimension,
                bool                                             aUseNewGhostSets,
                std::unordered_map< MSI::Dof_Type, moris_index > aDofTypeToBsplineMeshIndex )
                : Model_Initializer(
                          aParameterList,
                          aMeshPair,
                          aLibrary,
                          aSpatialDimension,
                          aUseNewGhostSets,
                          aDofTypeToBsplineMeshIndex ){};

      protected:

      public:
        void initialize() override;

      protected:
        void create_material_models() override;
        void create_constitutive_models() override;
        void create_stabilization_parameters() override;
        void create_iwgs() override;
        void create_iqis() override;
        void create_set_info() override;
        void print_physics_model() override;

      private:
        std::string
        read_phase_name( ParameterList const &aParameterList, mtk::Leader_Follower const aLeaderFollower ) const;

        Vector< std::pair< std::shared_ptr< Constitutive_Model >, std::string > >
        read_constitutive_models( ParameterList const &aParameterList, mtk::Leader_Follower const &tLeaderFollower );

        Vector< std::pair< std::shared_ptr< Stabilization_Parameter >, std::string > >
        read_stabilization_parameters( ParameterList const &aParameterList );

        void create_phases();

        void get_mesh_set_names(
                fem::Element_Type      aBulkType,
                const std::string     &aLeaderPhaseName,
                const std::string     &aFollowerPhaseName,
                const std::string     &aNeighborPhasesString,
                const std::string     &aSideOrdinalsString,
                bool                   aIsGhost,
                Vector< std::string > &aMeshSetNames );

        std::map< std::string, std::shared_ptr< fem::Phase_User_Info > > mPhases;
    };
}    // namespace moris::fem


#endif    // MORIS_CL_FEM_MODEL_INITIALIZER_PHASEBASED_HPP
