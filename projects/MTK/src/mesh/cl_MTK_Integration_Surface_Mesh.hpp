/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Surface_Mesh.hpp
 *
 */

#pragma once

#include "fn_MTK_Integration_Surface_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_XTK_Enums.hpp"
#include "moris_typedefs.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_Json_Object.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include <ostream>

namespace moris::mtk
{
    /**
     * @brief This class is used to extract (or better: to provide a view of) a surface mesh from a given integration mesh.
     * This class does not store any vertex or cell data but only provides the necessary information to access the data in the integration mesh.
     * The surface mesh will use local indices that only refer to the vertices/facets in the surface mesh.
     * A corresponding mapping between the global and local indices is provided.
     */
    class Integration_Surface_Mesh : public Surface_Mesh
    {

      public:    // constructors
        Integration_Surface_Mesh(
                Integration_Surface_Mesh_Data const &aData );

        // methods

        [[nodiscard]] Matrix< DDRMat > initialize_vertex_coordinates( Integration_Mesh const *aIGMesh );

        [[nodiscard]] Matrix< DDRMat > initialize_vertex_coordinates_from_side_sets(
                Integration_Mesh const           *aIGMesh,
                const Vector< Side_Set const * > &aSideSets );

        [[nodiscard]] Vector< Vector< moris_index > > get_cell_to_vertex_indices(
                Integration_Mesh const           *aIGMesh,
                const Vector< Side_Set const * > &aSideSets ) const;

        Vector< Side_Set const * >
        obtain_sidesets_from_names( Integration_Mesh_DataBase_IG const *aIGMesh, const Vector< std::string > &aSideSetNames );

        void set_all_displacements( const Matrix< DDRMat > &aDisplacements ) override;

        /**
         * @brief Returns the indices of all neighboring vertices for each vertex in the surface mesh.
         * @return A list of lists. The outer list contains the neighbor-lists for each vertex. The inner list contains the indices of the neighbors.
         */
        [[nodiscard]] Vector< Vector< moris_index > > get_vertex_neighbors() const;

        /**
         * @brief Returns the indices of all neighboring vertices for the vertex with the given local index.
         * @param aLocalVertexIndex The local index of the vertex in the surface mesh.
         * @return A list of indices of the neighbors of the vertex with the given index.
         */
        [[nodiscard]] Vector< moris_index > get_vertex_neighbors( moris_index aLocalVertexIndex ) const;


        /**
         * @brief Returns the facet measure (length/area) for each facet in the surface mesh.
         * @return A (n x 1) matrix where n is the number of facets in the surface mesh.
         */
        [[nodiscard]] const Matrix< DDRMat > &get_facet_measure() const;

        /**
         * @brief Returns the averaged vertex normals for each vertex in the surface mesh.
         * @details The vertex normals are averaged over all facets that are connected to the vertex, weighted by the respective facet measure.
         * @return A (d x n) matrix where d is the dimension of the mesh (holding the normal components) and n is the number of vertices in the surface mesh.
         */
        [[nodiscard]] Matrix< DDRMat > get_vertex_normals() const;

        /**
         * @brief Returns the global index of a vertex with the given local index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aLocalVertexIndex The local index of the vertex in the surface mesh.
         * @return The global index of the vertex in the whole mesh.
         */
        [[nodiscard]] moris_index get_global_vertex_index( moris_index aLocalVertexIndex ) const;

        /**
         * @brief Returns the global index of a cell with the given local index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aLocalCellIndex The local index of the cell in the surface mesh.
         * @return The global index of the cell in the whole mesh.
         */
        [[nodiscard]] moris_index get_global_cell_index( moris_index aLocalCellIndex ) const;

        /**
         * @brief Returns the local index of a vertex with the given global index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aGlobalVertexIndex The global index of the vertex in the whole mesh.
         * @return The local index of the vertex in the surface mesh.
         */
        [[nodiscard]] moris_index get_local_vertex_index( moris_index aGlobalVertexIndex ) const;

        [[nodiscard]] moris_index get_cluster_of_cell( moris_index aLocalCellIndex ) const;

        /**
         * @brief Returns all local vertices that are part of the cell with the given local index.
         * @param aLocalCellIndex The local index of the cell in the surface mesh.
         * @return A list of local vertex indices that are part of the cell with the given index.
         */
        [[nodiscard]] Vector< moris_index > get_vertices_of_cell( moris_index aLocalCellIndex ) const;

        /**
         * @brief Returns all local cells that are neighbors of the vertex with the given local index.
         * @param aLocalVertexIndex The local index of the vertex in the surface mesh.
         * @return A list of local cell indices that are neighbors of the vertex with the given index.
         */
        [[nodiscard]] Vector< moris_index > get_cells_of_vertex( moris_index aLocalVertexIndex ) const;

        /**
         * @brief Returns the (averaged) vertex normals of all vertices that are part of the cell with the given local index.
         * @param aLocalCellIndex The local index of the cell in the surface mesh.
         * @return A (d x n) matrix where d is the dimension of the mesh and n is the number of vertices in the cell.
         */
        [[nodiscard]] Matrix< DDRMat > get_vertex_normals_of_cell( moris_index aLocalCellIndex ) const;

        /**
         * @brief Returns the local index of a cell with the given global index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aGlobalCellIndex The global index of the cell in the whole mesh.
         * @return The local index of the cell in the surface mesh.
         */
        [[nodiscard]] moris_index get_local_cell_index( moris_index aGlobalCellIndex ) const;

        [[nodiscard]] uint get_spatial_dimension() const override;

        Json to_json() const;

        //--------------------------------------------------------------------------------
        // XQI Related functions
        //--------------------------------------------------------------------------------

        real compute_XQI( xtk::XQI_Type aType ) const;

        //--------------------------------------------------------------------------------

        /**
         * Computes the sensitivities of the requested XQI type wrt to PDVs.
         *
         * @param aType The type of XQI for which the sensitivities are requested.
         * @param aVertexPDVIDs The PDV IDs associated with each vertex in the surface mesh. -1 if no PDV is associated with the vertex.
         * @param aSensitivities The distributed vector where the sensitivities will be stored. It is assumed that this vector is already initialized and has the correct map.
         * @param aRequestIndex The vector index in the Dist_Vector to store the sensitivities.
         */
        void compute_XQI_sensitivities(
                const xtk::XQI_Type                    aType,
                const Vector< Vector< moris_index > > &aVertexPDVIDs,
                sol::Dist_Vector                      *aSensitivities,
                const uint                             aRequestIndex ) const;

      private:    // methods
        void initialize_facet_measure();

        void initialize_vertex_normals();

        /**
         * @brief Contains information about the surface mesh including mapping to the original IG mesh and other useful maps
         */
        Integration_Surface_Mesh_Data mData;

        /**
         * @brief Stores the averaged vertex normals for each vertex in the surface mesh. The indices are the indices of the vertices in the surface mesh, not the global indices!
         */
        Matrix< DDRMat > mVertexNormals = Matrix< DDRMat >( 0, 0 );

        /**
         * @brief Is used to store the measure of each facet.
         */
        Matrix< DDRMat > mFacetMeasure = Matrix< DDRMat >( 0, 0 );
    };
}    // namespace moris::mtk
