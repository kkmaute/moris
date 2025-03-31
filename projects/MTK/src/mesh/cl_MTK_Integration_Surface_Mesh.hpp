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

#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "moris_typedefs.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_Json_Object.hpp"
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
                Integration_Mesh const      *aIGMesh,
                const Vector< std::string > &aSideSetNames );

        Integration_Surface_Mesh(
                Integration_Mesh const          *aIGMesh,
                Vector< mtk::Side_Set const * > &aSideSets );

        Integration_Surface_Mesh(
                Integration_Mesh_DataBase_IG const *aIGMesh,
                const Vector< std::string >        &aSideSetNames );

        // methods
        
        [[nodiscard]] Matrix< DDRMat > initialize_vertex_coordinates( Integration_Mesh const *aIGMesh );

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

      private:    // methods
        /**
         * @brief This function initializes the surface mesh from an IG Mesh and side set names. It is called by the constructor.
         * @details The surface mesh is initialized only once. During the initialization, the mapping between the global
         * to the local indices is created. The local indices are the indices of the vertices/facets in the surface mesh (starting at 0).
         */
        void initialize_from_side_sets( const Vector< mtk::Side_Set const * > &aSideSets );

        /**
         * @brief Assigns the neighbors to each vertex in the surface mesh after the first pass over all vertices.
         * This has to be done after all vertices have already been added to the list of vertices as the local indices of neighbors might otherwise not be known.
         * @param aTmpNeighborMap Map from the global vertex index to the neighbors of the vertex (also global indices). This map is only required for the initialization
         * of the surface mesh and is created in the initialization call.
         */
        void initialize_neighbors(
                map< moris_index, Vector< moris_index > > &aTmpNeighborMap );

        /**
         * @brief Initializes the given side set (calls the cluster initializer for each cluster in the side set).
         * @param aSideSet
         */
        void initialize_side_set( map< moris_index, Vector< moris_index > > &aTmpNeighborMap, Set const *aSideSet );

        /**
         * @brief Initializes the given cluster (calls the cell initializer for each cell in the cluster).
         * @param aTmpNeighborMap
         * @param aCluster
         */
        void initialize_cluster( map< moris_index, Vector< moris_index > > &aTmpNeighborMap, Cluster const *const &aCluster, moris_index aClusterIndex );

        /**
         * @brief Initializes all information for the given cell.
         * @param aCell
         */
        void initialize_cell( map< moris_index, Vector< moris_index > > &aTmpNeighborMap, const Cell *aCell, int aCellOrdinal, moris_index aClusterIndex );

        /**
         * @brief Initializes all information for the given vertex.
         * @param aTmpNeighborMap
         * @param aCurrentLocalCellIndex
         * @param aSideVertices
         * @param aVertex
         */
        void initialize_vertex(
                map< moris_index, Vector< moris_index > > &aTmpNeighborMap,
                moris_index                                aCurrentLocalCellIndex,
                Vector< Vertex const * >                  &aSideVertices,
                Vertex const                              *aVertex );

        void initialize_vertex_coordinates();

        void initialize_facet_measure();

        void initialize_vertex_normals();

        // data
        Integration_Mesh_DataBase_IG const *mIGMesh;

        /**
         * @brief Map from the global vertex index to the index of the vertex in the surface mesh
         * @example If a (global) mesh consists of vertices from 0 to n but the surface mesh only consists of some of
         * those vertices (e.g. vertex 3, 18, 5, 20, ...), the map will provide mappings from 3 to 0, 18 to 1, 5 to 2 and 20
         * to 3 (and so on). It is not ensured that the indices will map to vertices in ascending order!
         */
        moris::map< moris_index, moris_index > mGlobalToLocalVertexIndex;

        /**
         * @brief The value at the n-th (local) index is the global index of the vertex in the global mesh.
         * global index of the vertex in the global mesh.
         * @example If a (global) mesh consists of vertices from 0 to n but the surface mesh only consists of some of
         * those vertices (e.g. vertex 3, 18, 5, 20, ...), the list will contain the vertices 3, 18, 5, 20, ... in that
         * order.
         */
        Vector< moris_index > mLocalToGlobalVertexIndex;

        /**
         * @brief List of neighboring vertices for each vertex in the surface mesh. The indices are the indices of the
         * vertices in the surface mesh, not the global indices!
         */
        Vector< Vector< moris_index > > mVertexNeighbors;

        /**
         * @brief Map from the global cell index to the index of the cell in the surface mesh
         * @example If a (global) mesh consists of cells from 0 to n but the surface mesh only consists of some of
         * those cells (e.g. cells 3, 18, 5, 20, ...), the map will provide mappings from 3 to 0, 18 to 1, 5 to 2 and 20
         * to 3 (and so on). It is not ensured that the indices will map to cells in ascending order!
         */
        moris::map< moris_index, moris_index > mGlobalToLocalCellIndex;

        /**
         * @brief The value at the n-th (local) index is the global index of the cell in the global mesh.
         * @example If a (global) mesh consists of cells from 0 to n but the surface mesh only consists of some of
         * those cells (e.g. cells 3, 18, 5, 20, ...), the list will contain the cells 3, 18, 5, 20, ... in that
         * order.
         */
        Vector< moris_index > mLocalToGlobalCellIndex;

        /**
         * @brief Stores the side ordinal for each cell which is used to determine the surface facet of the cell.
         * @example For a cell with three sides, that coincides with the surface mesh on side ordinal 1, the value at the
         * index of the cell in the surface mesh will be 1.
         */
        Vector< moris_index > mCellSideOrdinals;

        /**
         * @brief List of vertices that are part of the cell with the given index. The indices are the
         * indices of the vertices in the surface mesh, not the global indices!
         */
        Vector< Vector< moris_index > > mCellToVertexIndices;

        /**
         * @brief List of cell indices that the vertex with the given index is part of. The indices are the indices of
         * the cell in the surface mesh, not the global indices!
         */
        Vector< Vector< moris_index > > mVertexToCellIndices;

        Vector< mtk::Side_Set const * > mSideSets;

        Vector< Vector< moris_index > > mSideSetToClusterIndices;

        Vector< Vector< moris_index > > mClusterToCellIndices;

        Vector< moris_index > mCellToClusterIndices;

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
