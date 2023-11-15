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
#include <ostream>

namespace moris::mtk
{
    //    struct Surface_Vertex {
    //        mtk::Vertex const *mVertex;
    //        moris::Cell<Surface_Vertex *> mNeighbors;
    //        Matrix<moris::DDRMat> mCoords;
    //    };
    //
    //    struct Surface_Facet {
    //        moris::Cell<Surface_Vertex *> mVertices;
    //        moris::Cell<Surface_Facet *> mNeighbors;
    //        moris::Matrix<moris::DDRMat> mNormal{1, 3, 0.0f};
    //    };


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
        moris::Cell<moris_index> mCellSideOrdinals;

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

      public:
        Surface_Mesh( Integration_Mesh *aIGMesh, const moris::Cell< std::string > &aSideSetNames );

        Surface_Mesh( Integration_Mesh_DataBase_IG *aIGMesh, const moris::Cell< std::string > &aSideSetNames );

        [[nodiscard]] Matrix< DDRMat > get_vertex_coordinates() const;

        [[nodiscard]] moris::Cell< moris::Cell< moris_index > > get_vertex_neighbors() const;

        [[nodiscard]] Matrix< DDRMat > get_facet_normals() const;

        [[nodiscard]] Matrix< DDRMat > get_vertex_normals() const;

      private:
        /**
         * @brief This function initializes the surface mesh. It is called by the constructor.
         * @details The surface mesh is initialized only once. During the initialization, the mapping between the global
         * to the local indices is created. The local indices are the indices of the vertices/facets in the surface mesh (starting at 0).
         * @param aSideSetNames All side sets for which the surface mesh should be extracted.
         */
        void initialize_surface_mesh( moris::Cell< std::string > const &aSideSetNames );

        /**
         * @brief Small helper function that is used in the initializer
         * @param aCurrentCell
         */
        void initialize_cell( const Cell *aCurrentCell, int aCurrentCellOrdinal );

        /**
         * @brief Small helper function that is used in the initializer
         * @param aTmpNeighborMap
         * @param aCell
         * @param aLocalCellIndex
         * @param aCellOrdinal
         */
        void initialize_vertices(
                map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap,
                const Cell                                     *aCell,
                moris_index                                     aLocalCellIndex,
                int                                             aCellOrdinal );

        /**
         * @brief Assigns the neighbors to each vertex in the surface mesh after the first pass over all vertices.
         * This has to be done after all vertices have already been added to the list of vertices as the local indices of neighbors might otherwise not be known.
         * @param aTmpNeighborMap Map from the global vertex index to the neighbors of the vertex (also global indices). This map is only required for the initialization
         * of the surface mesh and is created in the initialization call.
         */
        void initialize_neighbors( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap );
    };


}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_CL_MTK_SURFACE_MESH_HPP_ */