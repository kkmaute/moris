/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Integration_Mesh_Cleanup.hpp
 *
 */

#include "cl_Cell.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "typedefs.hpp"
#include "cl_GEN_Geometry_Engine.hpp"

#ifndef SRC_cl_XTK_Integration_Mesh_Cleanup_HPP_
#define SRC_cl_XTK_Integration_Mesh_Cleanup_HPP_


using namespace moris;

namespace xtk
{
    // class Cut_Integration_Mesh;

    class Integration_Mesh_Cleanup
    {
        // ----------------------------------------------------------------------------

      protected:
        Cut_Integration_Mesh*     mInputMesh;
        Facet_Based_Connectivity* mInputFacetConnectivity;

        // ----------------------------------------------------------------------------

        moris::moris_index                mMergeNum;     // Number of vertices to merge
        moris::Cell< moris::moris_index > mMergeInds;    // Indices of vertices of merge
        moris::Cell< moris::moris_index > mMergeIndsV2;
        moris::Cell< moris::moris_index > mCellMergeInds;    // Indices of cells merged
        moris::Cell< moris::moris_index > mFacetMergeInds;

        moris::Cell< bool_t > mChildMeshBoundary;    // Flag if vertex index is on child mesh boundary
        moris::Cell< bool_t > mBlkPhaseBoundary;     // Flag if vertex index is on material boundary
        moris::Cell< bool_t > mOwnedVertices;        // Flag if vertex is owned by processor

        moris::Cell< moris_id > mNotOwnedIgCellIds;

        moris::uint mNumVerts;    // Total number of vertices


        moris::Cell< std::shared_ptr< IG_Vertex_Group > > mVertexGroups;
        moris::Cell< std::shared_ptr< IG_Cell_Group > >   mCellGroups;

        moris::Cell< moris::Cell< moris::moris_index > > mVertIndToCells;    // reference from a vertex index to all connected cells
        moris::Cell< moris::Cell< moris::moris_index > > mVertIndToVerts;    // reference from a vertex index to all connected vertices

        moris::Cell< moris::moris_index > mFlats;

        moris::uint mNumVertGroups;
        moris::uint mNumCellGroups;

        moris::real tMergeTimeTime;
        moris::real tFacetTime;

        // ----------------------------------------------------------------------------

      public:

        // ----------------------------------------------------------------------------

        Integration_Mesh_Cleanup() = default;

        // ----------------------------------------------------------------------------

        /**
         * @brief Construct a new class Integration_Mesh_Cleanup
         *
         * @param aXTKMesh Cut_Integration_Mesh class to alter
         * @param aFacetConnectivity Facet_Based_Connectivity class to alter
         *
         */
        Integration_Mesh_Cleanup(
                Cut_Integration_Mesh*     aXTKMesh,
                Facet_Based_Connectivity* aFacetConnectivity );


        // ----------------------------------------------------------------------------

        virtual ~Integration_Mesh_Cleanup();

        // ----------------------------------------------------------------------------

        /**
         * @brief set member variable mVertIndToCells which contains list of attached cells for each vertex index
         */
        void
        set_vertex_cell_connectivity();

        // ----------------------------------------------------------------------------

        /**
         * @brief set member variable mVertIndToVerts which contains list of attached vertices for each vertex index
         */
        void
        set_vertex_vertex_connectivity();

        // ----------------------------------------------------------------------------

        void
        adjust_vertex_vertex_connectivity( moris_index aVert1, moris_index aVert2 );

        // ----------------------------------------------------------------------------

        /**
         * @brief set member variable mBlkPhaseBoundary true/false if vertex is on bulk phase boundary
         */
        void
        set_blk_phase_vertex_flags();

        // ----------------------------------------------------------------------------

        /**
         * @brief set member variable mChildMeshBoundary true/false if vertex is on child mesh boundary
         */
        void
        set_child_mesh_boundary_flag();

        // ----------------------------------------------------------------------------

        /**
         * @brief compare vertex ids between processors
         * @param aNotOwnedCellVerts not owned cell vertices
         */
        void
        communicate_merged_cells( moris::Matrix< IdMat > aNotOwnedCellVerts );

        // ----------------------------------------------------------------------------

        /**
         * @brief set vertex index map between GEN and FEM
         * @param aGeometryEngine link between XTK and FEM
         */
        void
        make_GenMeshMap( moris::ge::Geometry_Engine* aGeometryEngine );

        // ----------------------------------------------------------------------------

        /**
         * @brief sets flag for if vertex is owned by processor
         */
        void
        set_owned_verts_flag();

        // ----------------------------------------------------------------------------

        /**
         * @brief set member variable mMergeNum, which is the number of vertices to be merged.
         * @brief set member variable mMergeInds, which is a list of indices of the vertices to be merged
         */
        void
        num_merges();

        // ----------------------------------------------------------------------------

        /**
         * @brief returns index of vertex that Vert1 should be merged with, using
         * delauney condition, and if it results in a flat or inverted cell
         * @param Vert1 Vertex index to merge
         */
        moris_index
        get_attached_vertex( moris_index aVert1 );

        // ----------------------------------------------------------------------------

        /**
         * @brief performs mesh simplification
         * @param aActiveIgCells Integration Cells to be altered
         * @param aGeometryEngine Geometry engine to set a GEN-MESH vertex index map
         */
        void
        perform(
                moris::Cell< moris::mtk::Cell* >& aActiveIgCells,
                moris::ge::Geometry_Engine*       aGeometryEngine );

        // ----------------------------------------------------------------------------

        /**
         * @brief return delauny fitness value of merge (minimum angle of transformed cells)
         *
         * @param aVert1 Vertex index to merge to Vert2
         * @param aVert2 Stationary vertex index being merged
         */
        double_t
        delauny(
                moris_index aVert1,
                moris_index aVert2 );

        // ----------------------------------------------------------------------------

        /**
         * @brief return false if merge inverts a triangle or tet
         *
         * @param aVert1 Vertex index to merge to Vert2
         * @param aVert2 Stationary vertex index being merged
         */
        bool_t
        calcInvert(
                moris_index aVert1,
                moris_index aVert2 );


        // ----------------------------------------------------------------------------

        /**
         * @brief return minimum angle of a triangle spanned by 3 coordinates
         *
         * @param C1 Coordinate 1
         * @param C2 Coordinate 2
         * @param C3 Coordinate 3
         */
        double_t
        minAngle2D(
                moris::Matrix< DDRMat > C1,
                moris::Matrix< DDRMat > C2,
                moris::Matrix< DDRMat > C3 );

        // ----------------------------------------------------------------------------

        /**
         * @brief return minimum solid angle of a tetrahedron spanned by 4 coordinates
         *
         * @param C1 Coordinate 1
         * @param C2 Coordinate 2
         * @param C3 Coordinate 3
         * @param C4 Coordinate 4
         */
        double_t
        minAngle3D(
                moris::Matrix< DDRMat > C1,
                moris::Matrix< DDRMat > C2,
                moris::Matrix< DDRMat > C3,
                moris::Matrix< DDRMat > C4 );

        // ----------------------------------------------------------------------------

        /**
         * @brief return solid angle between OA, OB, OC vectors
         *
         * @param node0 Point coordinates
         * @param nodeA Point coordinates
         * @param nodeB Point coordinates
         * @param nodeC Point coordinates
         */
        double_t
        solidAngle(
                moris::Matrix< DDRMat > nodeO,
                moris::Matrix< DDRMat > nodeA,
                moris::Matrix< DDRMat > nodeB,
                moris::Matrix< DDRMat > nodeC );

        // ----------------------------------------------------------------------------

        /**
         * @brief return angle between vectors A and B
         *
         * @param nodeA vector
         * @param nodeB vector
         */
        double_t
        angle(
                moris::Matrix< DDRMat > nodeA,
                moris::Matrix< DDRMat > nodeB );

        // ----------------------------------------------------------------------------

        /**
         * @brief acos calculation avoiding an input above 1 or below -1
         *
         * @param x angle to calculate acos of
         */
        double_t
        safeAcos( double_t x );

        // ----------------------------------------------------------------------------

        /**
         * @brief calculate area or volume of triangle or tetrahedron
         *
         * @param aCoords rows of coordinates spanning the shape/volume
         */
        double_t
        n_area( moris::Matrix< DDRMat > aCoords );

        // ----------------------------------------------------------------------------

        /**
         * @brief merge two vertices
         *
         * @param Vert1 Vertex index to merge
         * @param Vert2 Stationary vertex index in merge
         * @param aActiveIgCells group of cells to be altered
         */
        void
        merge( moris_index Vert1, moris_index Vert2, moris::Cell< moris::mtk::Cell* >& aActiveIgCells );

        // ----------------------------------------------------------------------------

        /**
         * @brief merges a list of vertices
         */
        void
        merge_list( moris::Cell< moris::mtk::Cell* >& aActiveIgCells );

        // ----------------------------------------------------------------------------

        /**
         * @brief corrects the indices of cells, vertices and facets
         */
        void
        shift_indices( moris::Cell< moris::mtk::Cell* >& aActiveIgCells );

        // ----------------------------------------------------------------------------

        /**
         * @brief returns index list of flat tris/tets
         */
        moris::Cell< moris_index >
        check_flats();

        // ----------------------------------------------------------------------------

        /**
         * @brief checks for coinciding vertices
         */
        void
        check_coinc_verts( moris::Cell< moris_index > aFlats, moris::Cell< moris::mtk::Cell* >& aActiveIgCells );

        // ----------------------------------------------------------------------------

        /**
         * @brief returns distance between two coordinates
         * @param aV1 Coordinates 1
         * @param aV2 Coordinates 2
         */
        double_t
        dist( moris::Matrix< DDRMat > aV1, moris::Matrix< DDRMat > aV2 );

        // ----------------------------------------------------------------------------

    };

    // ----------------------------------------------------------------------------

}    // namespace xtk

#endif    // SRC_cl_XTK_Integration_Mesh_Cleanup_HPP_