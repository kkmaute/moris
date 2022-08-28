/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Periodic_Boundary_Condition_Helper.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_PERIODIC_BOUNDARY_CONDITION_HELPER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_PERIODIC_BOUNDARY_CONDITION_HELPER_HPP_

#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_Param_List.hpp"
namespace moris
{
    namespace mtk
    {
        class Periodic_Boundary_Condition_Helper
        {
            public:
            Periodic_Boundary_Condition_Helper( std::shared_ptr<Mesh_Manager> aMeshManager,
                                                moris_index                   aMeshIndex,
                                                moris::ParameterList &        aParameterList);

            void
            setup_periodic_boundary_conditions();

            private:
            std::shared_ptr<Mesh_Manager>             mMeshManager;
            moris_index                               mMeshIndex;
            moris::Cell< moris::Cell< std::string > > mMeshSideSetPairs;
        };
    }
}

#endif
