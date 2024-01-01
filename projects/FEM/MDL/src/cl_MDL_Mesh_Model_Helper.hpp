/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MDL_Mesh_Model_Helper.hpp
 *
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_MDL_MESH_MODEL_HELPER_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_MDL_MESH_MODEL_HELPER_HPP_

#include "moris_typedefs.hpp"                       //MRS/COR/src
#include "cl_Vector.hpp"                        //MRS/CNT/src

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

    namespace MSI
    {
        class Model_Solver_Interface;
        class MSI_Solver_Interface;
        class Equation_Set;
        class Equation_Object;
    }

    namespace mdl
    {
//------------------------------------------------------------------------------

        class Mesh_Model_Helper
        {
        private:

            mtk::Mesh_Manager * mMesh = nullptr;

            Vector< moris::Matrix< DDUMat > > mColorListBlock;
            Vector< moris::Matrix< DDUMat > > mColorListSideSet;

            mtk::Interpolation_Mesh* mInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   mIntegrationMesh = nullptr;

            Vector< fem::Node_Base* >           mNodes;

            Vector< Matrix< DDSMat > > VertexIndOnColor;

            Matrix< DDSMat >  mVerticesOnBlock;
            Matrix< DDSMat >  mVerticesOnSideSet;

            Matrix< DDSMat >  mNodeToVertexIndMap;
            Vector< Matrix< DDSMat > >  mVertexColorToNodeIndMap;

            Vector< MSI::Equation_Object* >     mElements;
            Vector< MSI::Equation_Set * >      mElementBlocks;

            Vector< Cell < fem::Node_Base* > > mNodeSets;

            // color to physics map

            //color and vertexID

            void compute_max_num_vertices_on_color( Matrix< DDSMat > & aNumMaxVertPerColor );

            void compute_unique_vertex_list_on_color( const Matrix< DDSMat > & aNumMaxVertPerColor );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Mesh_Model_Helper(       mtk::Mesh_Manager * aMesh,
                               const moris_index         aMeshPairIndex);

//------------------------------------------------------------------------------

            ~Mesh_Model_Helper();

            void set_block_by_color( const uint             tColor,
                                     const Matrix< DDUMat > tBlockIndex );

//------------------------------------------------------------------------------

            void set_side_set_by_color( const uint             tColor,
                                        const Matrix< DDUMat > tSideSetIndex );

//------------------------------------------------------------------------------

            void compute_unique_node_lists();

//------------------------------------------------------------------------------

            void assign_node_ids();

//------------------------------------------------------------------------------

            void create_node_set( const moris::uint               aColor,
                                  const Matrix< moris::IndexMat > aNodeSetIndex);

//------------------------------------------------------------------------------

           void create_elements();

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//            void create_nodes();

            Vector< fem::Node_Base* > & get_nodes()
            {
                return mNodes;
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */

#endif /* PROJECTS_FEM_MDL_SRC_CL_MDL_MESH_MODEL_HELPER_HPP_ */

