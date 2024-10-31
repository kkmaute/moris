/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_FEM_Model_Initializer_Phasebased.hpp
 *
 */

#ifndef MORIS_CL_FEM_MODEL_INITIALIZER_PHASEBASED_HPP
#define MORIS_CL_FEM_MODEL_INITIALIZER_PHASEBASED_HPP

#include "cl_FEM_Model_Initializer.hpp"
#include "cl_FEM_Phase_User_Info.hpp"
#include <map>
#include <utility>
#include "cl_FEM_Set_User_Info.hpp"
namespace moris::fem
{
    class Model_Initializer_Phasebased : public Model_Initializer
    {
      public:
        Model_Initializer_Phasebased(
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
        void create_phases();

        void get_mesh_set_names(
                fem::Element_Type      aBulkType,
                const std::string     &aLeaderPhaseName,
                const std::string     &aFollowerPhaseName,
                const std::string     &aNeighborPhaseString,
                const std::string     &aSideOrdinalsString,
                bool                   aIsGhost,
                Vector< std::string > &aMeshSetNames );

        Vector< fem::Phase_User_Info > mPhaseInfo;
        std::map< std::string, uint >  mPhaseMap;
    };
}    // namespace moris::fem

#endif    // MORIS_CL_FEM_MODEL_INITIALIZER_PHASEBASED_HPP
