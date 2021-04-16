/*
 * cl_XTK_Mesh_Cleanup.hpp
 *
 *  Created on: Mar 25, 2021
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_MESH_CLEANUP_HPP_
#define SRC_XTK_CL_XTK_MESH_CLEANUP_HPP_

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
        std::string mMethod;
        moris::real mSnapTolerance;
    };

    struct Mesh_Quality_Data
    {
        moris::Matrix<moris::DDRMat>              mVolumes;
        moris::Cell<moris::Matrix<moris::DDRMat>> mEdgeLengths; // cell here because of various number of edges
    };

    class Mesh_Cleanup
    {
        public:
        Mesh_Cleanup(Model*                aModel,
                     moris::ParameterList* aParamList);

        void
        perform();


        // metric computations
        void
        compute_mesh_quality(Mesh_Quality_Data & aMeshQuality);

        void
        compute_mesh_volumes(moris::Matrix<moris::DDRMat> & aVolumes);

        void 
        compute_mesh_side_lengths(moris::Cell<moris::Matrix<moris::DDRMat>> & aEdgeLengths);
        

        
        private:
        Model*                  mModel;
        Mesh_Cleanup_Parameters mMeshCleanupParameters;

        void
        cleanup_cut_mesh();

        void
        finalize_cleanup();

    };
}


#endif