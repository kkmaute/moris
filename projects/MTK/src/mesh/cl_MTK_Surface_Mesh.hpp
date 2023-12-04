/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Surface_Mesh.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_SURFACE_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_SURFACE_MESH_HPP_

#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "typedefs.hpp"
#include "cl_MTK_Side_Set.hpp"
#include <ostream>

namespace moris::mtk
{

    /**
     * @brief This class is used to extract (or better: to provide a view of) a surface mesh from a given integration mesh.
     * This class does not store any vertex or cell data but only provides the necessary information to access the data in the integration mesh.
     * The surface mesh will use local indices that only refer to the vertices/facets in the surface mesh.
     * A corresponding mapping between the global and local indices is provided.
     */
    class Surface_Mesh
    {
      private:
        std::shared_ptr< Integration_Mesh_DataBase_IG > mIGMesh;

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
        moris::Cell< moris_index > mLocalToGlobalVertexIndex;

        /**
         * @brief List of neighboring vertices for each vertex in the surface mesh. The indices are the indices of the
         * vertices in the surface mesh, not the global indices!
         */
        moris::Cell< moris::Cell< moris_index > > mVertexNeighbors;

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
        moris::Cell< moris_index > mLocalToGlobalCellIndex;

        /**
         * @brief Stores the side ordinal for each cell which is used to determine the surface facet of the cell.
         * @example For a cell with three sides, that coincides with the surface mesh on side ordinal 1, the value at the
         * index of the cell in the surface mesh will be 1.
         */
        moris::Cell< moris_index > mCellSideOrdinals;

        /**
         * @brief List of vertices that are part of the cell with the given index. The indices are the
         * indices of the vertices in the surface mesh, not the global indices!
         */
        moris::Cell< moris::Cell< moris_index > > mCellToVertexIndices;

        /**
         * @brief List of cell indices that the vertex with the given index is part of. The indices are the indices of
         * the cell in the surface mesh, not the global indices!
         */
        moris::Cell< moris::Cell< moris_index > > mVertexToCellIndices;

        moris::Cell< mtk::Side_Set * > mSideSets;

        moris::Cell< moris::Cell< moris_index > > mSideSetToClusterIndices;

        moris::Matrix< DDRMat > mDisplacements;

        /**
         * @brief Is used to store the normals of the facet. This member variable should not be accessed directly but via the getter function get_facet_normals().
         * It is only initialized when the getter function is called for the first time. Subsequent calls will return the already computed value.
         */
        Matrix< DDRMat >
                mFacetNormals = Matrix< DDRMat >( 0, 0 );

        /**
         * @brief Is used to store the measure of the facet. This member variable should not be accessed directly but via the getter function get_facet_measure().
         * It is only initialized when the getter function is called for the first time. Subsequent calls will return the already computed value.
         */
        Matrix< DDRMat > mVertexNormals = Matrix< DDRMat >( 0, 0 );

        /**
         * @brief Is used to store the measure of the facet. This member variable should not be accessed directly but via the getter function get_facet_measure().
         * It is only initialized when the getter function is called for the first time. Subsequent calls will return the already computed value.
         */
        Matrix< DDRMat > mFacetMeasure = Matrix< DDRMat >( 0, 0 );


      public:
        Surface_Mesh( Integration_Mesh *aIGMesh, const moris::Cell< std::string > &aSideSetNames );

        explicit Surface_Mesh( moris::Cell< mtk::Side_Set * > aSideSets );

        Surface_Mesh( Integration_Mesh_DataBase_IG *aIGMesh, const moris::Cell< std::string > &aSideSetNames );

        /**
         * @brief Returns the coordinates of all vertices in the surface mesh.
         * @return A (n x d) matrix where n is the number of vertices in the surface mesh and d is the dimension of the mesh.
         */
        [[nodiscard]] Matrix< DDRMat > get_vertex_coordinates() const;

        void set_displacement( Matrix< DDRMat > aDisplacements );

        /**
         * @brief Returns the indices of all neighboring vertices for each vertex in the surface mesh.
         * @return A list of lists. The outer list contains the neighbor-lists for each vertex. The inner list contains the indices of the neighbors.
         */
        [[nodiscard]] moris::Cell< moris::Cell< moris_index > > get_vertex_neighbors() const;

        /**
         * @brief Returns the facet normals for each facet in the surface mesh.
         * @return A (n x d) matrix where n is the number of facets in the surface mesh and d is the dimension of the mesh (holding the normal components).
         */
        [[nodiscard]] Matrix< DDRMat > get_facet_normals();

        /**
         * @brief Returns the facet measure (length/area) for each facet in the surface mesh.
         * @return A (n x 1) matrix where n is the number of facets in the surface mesh.
         */
        [[nodiscard]] Matrix< DDRMat > get_facet_measure();

        /**
         * @brief Returns the averaged vertex normals for each vertex in the surface mesh.
         * @details The vertex normals are averaged over all facets that are connected to the vertex, weighted by the respective facet measure.
         * @return A (n x d) matrix where n is the number of vertices in the surface mesh and d is the dimension of the mesh (holding the normal components).
         */
        [[nodiscard]] Matrix< DDRMat > get_vertex_normals();

        /**
         * @brief Returns the global index of a vertex with the given local index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aLocalVertexIndex The local index of the vertex in the surface mesh.
         * @return The global index of the vertex in the whole mesh.
         */
        [[nodiscard]] moris_index get_global_vertex_index( moris_index aLocalVertexIndex );

        /**
         * @brief Returns the global index of a cell with the given local index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aLocalCellIndex The local index of the cell in the surface mesh.
         * @return The global index of the cell in the whole mesh.
         */
        [[nodiscard]] moris_index get_global_cell_index( moris_index aLocalCellIndex );

        /**
         * @brief Returns the local index of a vertex with the given global index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aGlobalVertexIndex The global index of the vertex in the whole mesh.
         * @return The local index of the vertex in the surface mesh.
         */
        [[nodiscard]] moris_index get_local_vertex_index( moris_index aGlobalVertexIndex );

        /**
         * @brief Returns the local index of a cell with the given global index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aGlobalCellIndex The global index of the cell in the whole mesh.
         * @return The local index of the cell in the surface mesh.
         */
        [[nodiscard]] moris_index get_local_cell_index( moris_index aGlobalCellIndex );


      private:
        /**
         * @brief This function initializes the surface mesh from an IG Mesh and side set names. It is called by the constructor.
         * @details The surface mesh is initialized only once. During the initialization, the mapping between the global
         * to the local indices is created. The local indices are the indices of the vertices/facets in the surface mesh (starting at 0).
         */
        void initialize_from_side_sets( const moris::Cell< mtk::Side_Set * > &aSideSets );


        /**
         * @brief Assigns the neighbors to each vertex in the surface mesh after the first pass over all vertices.
         * This has to be done after all vertices have already been added to the list of vertices as the local indices of neighbors might otherwise not be known.
         * @param aTmpNeighborMap Map from the global vertex index to the neighbors of the vertex (also global indices). This map is only required for the initialization
         * of the surface mesh and is created in the initialization call.
         */
         void initialize_neighbors(
                 map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap );

         /**
          * @brief Initializes the given side set (calls the cluster initializer for each cluster in the side set).
          * @param aSideSet
          */
         void initialize_side_set( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap, Set const *aSideSet );

         /**
          * @brief Initializes the given cluster (calls the cell initializer for each cell in the cluster).
          * @param aTmpNeighborMap
          * @param aCluster
          */
         void initialize_cluster(
                 map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap,
                 Cluster const *const                           &aCluster );

         /**
          * @brief Initializes all information for the given cell.
          * @param aCell
          */
         void initialize_cell(
                 map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap,
                 const Cell                                     *aCell,
                 int                                             aCellOrdinal );

         /**
          * @brief Initializes all information for the given vertex.
          * @param aTmpNeighborMap
          * @param aCurrentLocalCellIndex
          * @param aSideVertices
          * @param aVertex
          */
         void initialize_vertex(
                 map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap,
                 moris_index                                     aCurrentLocalCellIndex,
                 moris::Cell< Vertex const * >                  &aSideVertices,
                 Vertex const                                   *aVertex );

         /**
          * @brief Get all cells that are part of the surface mesh.
          * @return A map from the global cell index to the cell.
          */
         map< moris_index, Cell_DataBase const * > get_cells();
    };


}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_CL_MTK_SURFACE_MESH_HPP_ */