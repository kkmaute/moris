/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Mesh_Cleanup.hpp
 *
 */

#ifndef SRC_XTK_CL_XTK_MESH_CLEANUP_HPP_
#define SRC_XTK_CL_XTK_MESH_CLEANUP_HPP_

#include <unordered_map>
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"

#include "cl_Param_List.hpp"

namespace xtk
{
    class Model;
}

using namespace moris;

namespace xtk
{
    struct Mesh_Cleanup_Parameters
    {
        public:
        bool        mDeactivateOneBPChildMeshes;
    };

    class Mesh_Cleanup
    {
        public:
        Mesh_Cleanup(Model*                aModel,
                     moris::ParameterList* aParamList);

        void
        perform();

        private:
        Model*                  mModel;
        Mesh_Cleanup_Parameters mMeshCleanupParameters;

        void
        cleanup_cut_mesh();

        void
        select_candidate_child_meshes_for_cleanup(std::unordered_map<moris_index,moris_index>  & aRemoveChildMeshes);

        void
        get_vector_of_child_meshes_for_removal(
            std::unordered_map<moris_index,moris_index> const & aChildMeshesToDeleteMap,
            moris::Cell<moris_index> & aChildMeshesToDelete);

        void
        finalize_cleanup();

    };
}

#endif
