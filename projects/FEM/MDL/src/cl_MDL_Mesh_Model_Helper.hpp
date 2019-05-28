/*
 * cl_MDL_Mesh_Model_Helper.hpp
 *
 *  Created on: Mai 22, 2018
 *      Author: Schmidt
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_MDL_MESH_MODEL_HELPER_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_MDL_MESH_MODEL_HELPER_HPP_

#include "typedefs.hpp"                       //MRS/COR/src
#include "cl_Cell.hpp"                        //MRS/CON/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris
{
//------------------------------------------------------------------------------
    namespace mtk
    {
    class Mesh_Manager;
    }

    namespace fem
    {
        class Node_Base;
    }

    namespace mdl
    {
//------------------------------------------------------------------------------

        class Mesh_Model_Helper
        {
        private:

            mtk::Mesh_Manager * mMesh = nullptr;

            //color and vertexID

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Mesh_Model_Helper( mtk::Mesh_Manager * aMesh );

//------------------------------------------------------------------------------

            ~Mesh_Model_Helper();

//------------------------------------------------------------------------------

            void compute_unique_vertex_lists();

//------------------------------------------------------------------------------

            void create_nodes();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_MDL_MESH_MODEL_HELPER_HPP_ */
