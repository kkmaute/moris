/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Mesh_Editor.hpp
 *
 */

#pragma once

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "cl_TOL_Memory_Map.hpp"

namespace moris::mtk
{
    class Vertex;
    class Vertex_Interpolation;

    struct Interpolation_Mesh_Info
    {
        // mMTK mesh vertices which is being copied here for convience
        moris::Cell< mtk::Vertex const* > mVertices;

        // they are stored consecutively for each vertex
        moris::Cell< mtk::Vertex_Interpolation const* > mVertexInterpolations;

        // cell to vertex connectivity for lagrange cells
        moris::Cell< moris_index > mCellToVertexIndicies;

        uint mSpatialDim;

        uint mNumCells;

        uint mNumInterpolations;

        uint mNumLocalInterpolations;
    };
}// namespace moris::mtk

namespace moris::mtk
{
    class Interpolation_Mesh;
    class Interpolation_Mesh_DataBase_IP;

    class Interpolation_Mesh_Editor
    {
      private:
        // an mtk mesh that will be worked on
        moris::mtk::Interpolation_Mesh& mInputMesh;

        // temporay data that can be deleted later
        Interpolation_Mesh_Info* mIPMeshInfo;

        Interpolation_Mesh_DataBase_IP* mOutputMesh;

      public:
        // ----------------------------------------------------------------------------

        Interpolation_Mesh_Editor( moris::mtk::Interpolation_Mesh& aMTKMesh );

        // ----------------------------------------------------------------------------

        ~Interpolation_Mesh_Editor();

        // ----------------------------------------------------------------------------

        Interpolation_Mesh_DataBase_IP*
        perform();

        // ----------------------------------------------------------------------------

        /**
         * @brief  generate vertex data including coordinates for lagrange nodes// T-matrices for bspline nodes
         *
         */
        void
        generate_vertex_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief  generate connectivity information
         *
         */
        void
        generate_cell_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief initialize the size of t-matrix data and vertex coordinates size
         *
         */
        void
        initialize_vertex_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief  initialize the cell to vertex connectivity informaiton
         *
         */

        void
        initialize_cell_data();

        // --------------------------------------------------------------------------

        /**
         * @brief Get the num cells
         *
         * @return uint
         */
        uint
        get_num_cells();

        // --------------------------------------------------------------------------

        /**
         * @brief Get the memory usage
         *
         * @return moris::Memory_Map
         */

        moris::Memory_Map
        get_memory_usage();

        // --------------------------------------------------------------------------

        /**
         * @brief delete the unused member data
         *
         */

        void
        free_memory();

        // --------------------------------------------------------------------------

        /**
         * @brief Create a enriched mesh indices o
         *
         */

        void
        create_enriched_mesh_indices();

        // ----------------------------------------------------------------------------
        /**
         * @brief Create a ip vertices (bspline vertices)
         *
         */

        void create_ip_vertices();

        // ----------------------------------------------------------------------------

        /**
         * @brief Create the vertices lagrange
         *
         */

        void create_vertices();

        // ----------------------------------------------------------------------------

        /**
         * @brief create lagrange ip cells
         *
         */

        void create_cells();

        // ----------------------------------------------------------------------------

        /**
         * @brief Create a communication table from the old mesh
         *
         */

        void
        create_communication_table();

        // ----------------------------------------------------------------------------

        /**
         * @brief Create the vertex glb id to loc vertex ind map object
         *
         */

        void
        create_vertex_glb_id_to_loc_vertex_ind_map();

        // ----------------------------------------------------------------------------
        /**
         * @brief Create adof maps from the old mesh for different bspline meshes
         *
         */
        void
        create_adof_map();

          ////////----------------------------------------------------------------------------
            // Checking and debugging functions
            ////////----------------------------------------------------------------------------

            // ----------------------------------------------------------------------------

            /**
             * @brief check the maps are equal
             *
             */

            bool
            check_maps();

            // ----------------------------------------------------------------------------

            /**
             * @brief check the vertices are equal
             *
             */

            bool
            check_vertices();

            // ----------------------------------------------------------------------------

            /**
             * @brief check the t-matrices are equal
             *
             */

            bool
            check_t_matrices();

            // ----------------------------------------------------------------------------

            /**
             * @brief check the cell are equal
             *
             */

            bool
            check_cells();

            // ----------------------------------------------------------------------------

            /**
             * @brief checks the vertices, cells and t-matrices
             *
             */
            void
            check_input_output_mesh();
    };
}// namespace moris::mtk
