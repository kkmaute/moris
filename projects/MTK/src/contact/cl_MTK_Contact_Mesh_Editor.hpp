/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Contact_Mesh_Editor.hpp
 *
 */

#ifndef MORIS_CL_MTK_CONTACT_MESH_EDITOR_HPP
#define MORIS_CL_MTK_CONTACT_MESH_EDITOR_HPP


#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_ContactMapper_Ray.hpp"
#include "cl_MTK_Integrator.hpp"
#include "fn_assert.hpp"

namespace moris::mtk
{
    class Contact_Mesh_Editor
    {
      public:    // member functions
        Contact_Mesh_Editor(
                Integration_Mesh_DataBase_IG*   aIGMesh,
                Interpolation_Mesh_DataBase_IP* aIPMesh,
                const Integrator&               aIntegrationRule,
                moris::Cell< Side_Set* >&       aCandidateSideSet )
                : mIGMesh( aIGMesh )
                , mIPMesh( aIPMesh )
                , mIntegrationRule( aIntegrationRule )
                , mSideSets( aCandidateSideSet )
                , mPointMapper( ContactMapper_Ray( aIGMesh, aCandidateSideSet, moris::Cell< std::pair< moris_index, moris_index > >() ) ){};

        void update_nonconformal_side_sets();

        void update_displacements( Matrix< DDRMat >& aDisplacements );

      private:    // member functions

      private:    // data
        Integration_Mesh_DataBase_IG*   mIGMesh;
        Interpolation_Mesh_DataBase_IP* mIPMesh;
        Integrator                      mIntegrationRule;
        moris::Cell< Side_Set* >&       mSideSets;
        ContactMapper_Ray               mPointMapper;
    };
}    // namespace moris::mtk
#endif    // MORIS_CL_MTK_CONTACT_MESH_EDITOR_HPP
