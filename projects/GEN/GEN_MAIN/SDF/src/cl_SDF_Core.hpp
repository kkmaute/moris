/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Core.hpp
 *
 */

#pragma once

#include "typedefs.hpp"    // COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SDF_Object.hpp"
#include "cl_SDF_Mesh.hpp"
#include "cl_SDF_Parameters.hpp"
#include "cl_SDF_Data.hpp"

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    class Core
    {
        Mesh& mMesh;
        Data& mData;

        uint mCandidateSearchDepth        = 1;
        real mCandidateSearchDepthEpsilon = 0.01;
        bool mVerbose;

        //-------------------------------------------------------------------------------

      public:
        //-------------------------------------------------------------------------------

        Core( Mesh&   aMesh,
                Data& aData,
                bool  aVerbose = false );

        //-------------------------------------------------------------------------------

        ~Core(){};

        //-------------------------------------------------------------------------------

        void
        set_candidate_search_depth( const uint aCandidateSearchDepth )
        {
            mCandidateSearchDepth = aCandidateSearchDepth;
        }

        //-------------------------------------------------------------------------------

        void
        set_candidate_search_epsilon( const real aCandidateSearchEpsilon )
        {
            mCandidateSearchDepthEpsilon = aCandidateSearchEpsilon;
        }

        //-------------------------------------------------------------------------------

        void
        calculate_raycast();

        //-------------------------------------------------------------------------------

        void
        calculate_raycast(
                Matrix< IndexMat >& aElementsAtSurface );

        //-------------------------------------------------------------------------------
        void
        calculate_raycast(
                Matrix< IndexMat >& aElementsAtSurface,
                Matrix< IndexMat >& aElementsInVolume );

        //-------------------------------------------------------------------------------

        void
        calculate_raycast_and_sdf( Matrix< DDRMat >& aSDF );

        //-------------------------------------------------------------------------------

        void
        calculate_raycast_and_sdf(
                Matrix< DDRMat >&   aSDF,
                Matrix< IndexMat >& aElementsAtSurface,
                Matrix< IndexMat >& aElementsInVolume );

        //-------------------------------------------------------------------------------

        void
        save_to_vtk( const std::string& aFilePath );

        //-------------------------------------------------------------------------------

        void
        save_unsure_to_vtk( const std::string& aFilePath );

        //-------------------------------------------------------------------------------

      private:
        //-------------------------------------------------------------------------------

        void
        voxelize_3D( const uint aAxis );

        //-------------------------------------------------------------------------------

        void
        voxelize_2D( const uint aAxis );

        //-------------------------------------------------------------------------------

        moris::Cell< Vertex* >
        set_candidate_list();

        //-------------------------------------------------------------------------------

        void
        calculate_udf( moris::Cell< Vertex* >& aCandidateList );

        //-------------------------------------------------------------------------------

        /**
         * Kehrwoche :
         * make sure that each vertex is really associated
         * to its closest triangle
         */
        void
        sweep();

        //-------------------------------------------------------------------------------

        void
        fill_sdf_with_values( Matrix< DDRMat >& aSDF );

        //-------------------------------------------------------------------------------

        /**
         * @brief Takes a point in 3D space and determines which triangles
         * in the positive and negative x direction the point could possibly intersect.
         * These triangles are added to mData.mCandidateFacets. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
         *
         * @param aPoint Point in space that lies within the bounding Y-Z plane of the candidate triangles
         */
        void
        preselect_triangles_x( const Matrix< DDRMat >& aPoint );

        //-------------------------------------------------------------------------------

        /**
         * @brief Takes a point in 3D space and determines which triangles
         * in the positive and negative y direction the point could possibly intersect.
         * These triangles are added to mData.mCandidateFacets. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
         *
         * @param aPoint Point in space that lies within the bounding X-Z plane of the candidate triangles
         */
        void
        preselect_triangles_y( const Matrix< DDRMat >& aPoint );

        //-------------------------------------------------------------------------------

        /**
         * @brief Takes a point in 3D space and determines which triangles
         * in the positive and negative z direction the point could possibly intersect.
         * These triangles are added to mData.mCandidateFacets. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
         *
         * @param aPoint Point in space that lies within the bounding X-Y plane of the candidate triangles
         */
        void
        preselect_triangles_z( const Matrix< DDRMat >& aPoint );

        //-------------------------------------------------------------------------------

        /**
         * @brief Takes a point in 2D space and determines which lines
         * in the positive aAxis direction the point could possibly intersect. The ray is being cast in the !aAxis direction
         * These lines are added to mData.mIntersectedFacets. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
         * Since preselection is sufficient for 2D, some of the lines are already marked as intersected. These are added to mData.mCandidateFacets
         *
         * @param aAxis 0 for lines in the x direction of aPoint, 1 for lines in the y direction of aPoint
         * @param aPoint Location in space that lies within the bounding aAxis coordinates of the candidate triangles
         */
        void
        preselect_lines(
                const uint              aAxis,
                const Matrix< DDRMat >& aPoint );

        //-------------------------------------------------------------------------------


        void
        intersect_triangles(
                const uint              aAxis,
                const Matrix< DDRMat >& aPoint );

        //-------------------------------------------------------------------------------

        /**
         * Takes all of potential facets (in mData.mIntersectedFacets) and computes the coordinate axis intersection location
         *  with the ray originating from aPoint
         *
         * @param aAxis coordinate axis to shoot ray in. 0 = x, 1 = y, 2 = z
         * @param aPoint spatial location of the origin of the ray
         */
        void
        intersect_ray_with_facets(
                const uint              aAxis,
                const Matrix< DDRMat >& aPoint );

        //-------------------------------------------------------------------------------

        void
        check_if_node_is_inside_triangles(
                const uint aAxis,
                const uint aNodeIndex );

        //-------------------------------------------------------------------------------

        void
        check_if_node_is_inside_lines(
                const uint aAxis,
                const uint aNodeIndex );
                
        //-------------------------------------------------------------------------------

        void
        calculate_candidate_points_and_buffer_diagonal();

        //-------------------------------------------------------------------------------

        void
        get_nodes_withing_bounding_box_of_triangle(
                Facet*                  aFacet,
                moris::Cell< Vertex* >& aNodes,
                moris::Cell< Vertex* >& aCandList );

        //-------------------------------------------------------------------------------

        /**
         * performs a floodfill in order to fix unsure nodes
         */
        void
        force_unsure_nodes_outside();

        //----------------------------------------------------------------------f---------

        void
        random_rotation();

        //-------------------------------------------------------------------------------

        void
        undo_rotation();
    };

    //-------------------------------------------------------------------------------
} /* namespace moris::sdf */
