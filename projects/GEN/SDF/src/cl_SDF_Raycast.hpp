/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_SDF_Raycaster.cpp
 *
 */

#pragma once

#include <fstream>

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "SDF_Tools.hpp"

#include "cl_SDF_Object.hpp"

namespace moris::sdf
{
    class Raycast
    {
      private:
        Object&               mObject;               // closed collection of facets that rays can hit               
        uint                  mDimension;            // number of coordinate axes for both the object and the cast point
        Matrix< DDRMat >      mPoint;                // global point coordinates, can change due to coordinate rotation
        Matrix< DDRMat >      mOriginalPoint;        // global point coordinates, needed to revert rotations that occur
        Matrix< DDRMat >      mFacetMinCoords;       // row = facet, column = min mDimension coordinate of the vertices of the facet
        Matrix< DDRMat >      mFacetMaxCoords;       // row = facet, column = max mDimension coordinate of the vertices of the facet
        moris::Cell< uint >   mCandidateFacets;      // index of facets that lie within the bounding box, but might not actually be intersected
        moris::Cell< Facet* > mIntersectedFacets;    // intersection has been found, coordinate of intersection still to be computed
        moris::Cell< real >   mCoordsK;              // intersection coordinates for all facets in mIntersectedFacets
        uint                  mPointIsInside;        // tribool for if the point is inside the object. 0 = outside, 1 = inside, 2 = unsure

      public:
        /**
         * Constructor for Raycast. Sets member data from input. Determines min and max facet coords and sets them accordingly. 
         * 
         * @param aObject closed collection of facets and vertices. generated from .obj or .stl file
         */
        Raycast( Object& aObject );

        /**
         * Conducts a raycast to determine if aPoint is inside or outside of the object. Casts rays in each coordinate direction until the point is resolved.
         * Also determines intersection locations along the coordinate axes of aPoint with each facet.
         *
         * @param aPoint coordinate point in which the ray will originate
         */
        void
        raycast_point( const Matrix< DDRMat >& aPoint );

        //-------------------------------------------------------------------------------

        /**
         * Performs a raycast in the aAxis direction and determines whether the point is inside or outside.
         * Only to be used for 3D object files.
         *
         * @param aAxis Coordinate axis along which the ray will be cast
         */
        void
        voxelize_3D( const uint aAxis );

        //-------------------------------------------------------------------------------

        /**
         * Performs a raycast in the aAxis direction and determines whether the point is inside or outside.
         * Only to be used for 2D object files.
         *
         * @param aAxis Coordinate axis along which the ray will be cast
         */
        void
        voxelize_2D( const uint aAxis );


        /**
         * @brief Takes a point in 3D space and determines which triangles
         * in the positive and negative x direction the point could possibly intersect.
         * These triangles are added to mCandidateFacets. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
         *
         * @param aPoint Point in space that lies within the bounding Y-Z plane of the candidate triangles
         */
        void preselect_triangles_x();

        //-------------------------------------------------------------------------------

        /**
         * @brief Takes a point in 3D space and determines which triangles
         * in the positive and negative y direction the point could possibly intersect.
         * These triangles are added to mCandidateFacets. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
         *
         * @param aPoint Point in space that lies within the bounding X-Z plane of the candidate triangles
         */
        void
        preselect_triangles_y();

        //-------------------------------------------------------------------------------

        /**
         * @brief Takes a point in 3D space and determines which triangles
         * in the positive and negative z direction the point could possibly intersect.
         * These triangles are added to mCandidateFacets. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
         *
         * @param aPoint Point in space that lies within the bounding X-Y plane of the candidate triangles
         */
        void
        preselect_triangles_z();

        //-------------------------------------------------------------------------------

        /**
         * @brief Takes a point in 2D space and determines which lines
         * in the positive aAxis direction the point could possibly intersect. The ray is being cast in the !aAxis direction
         * These lines are added to mIntersectedFacets. Helps speed up the raycast by avoiding checking intersections with unrelated facets.
         * Since preselection is sufficient for 2D, some of the lines are already marked as intersected. These are added to mCandidateFacets
         *
         * @param aAxis 0 for lines in the x direction of aPoint, 1 for lines in the y direction of aPoint
         * @param aPoint Location in space that lies within the bounding aAxis coordinates of the candidate triangles
         */
        void
        preselect_lines( uint aAxis );

        //-------------------------------------------------------------------------------


        /**
         * checks if the cast point will cast a ray that will actually intersect mCandidateFacets. If so, the facet is added to mIntersectedFacets
         * 
         * @param aAxis coordinate axis in which to cast from
         */
        void
        intersect_triangles( uint aAxis );

        //-------------------------------------------------------------------------------

        /**
         * Takes all of potential facets (in mIntersectedFacets) and computes the coordinate axis intersection location
         *  with the ray originating from aPoint
         *
         * @param aAxis coordinate axis to shoot ray in. 0 = x, 1 = y, 2 = z
         * @param aPoint spatial location of the origin of the ray
         */
        void
        intersect_ray_with_facets( uint aAxis );

        //-------------------------------------------------------------------------------

        /**
         * determines if mPoint is within the closed collection of triangles. sets mPointIsInside accordingly
         * NOTE: Only to be used with 3D raycasts
         * 
         * @param aAxis coordinate axis in which the ray was cast from
         */
        void
        check_if_node_is_inside_triangles( uint aAxis );

        //-------------------------------------------------------------------------------

        /**
         * determines if mPoint is within the closed collection of lines. sets mPointIsInside accordingly
         * NOTE: Only to be used with 2D raycasts
         * 
         * @param aAxis coordinate axis in which the ray was cast from
         */
        void
        check_if_node_is_inside_lines( uint aAxis );

        //-------------------------------------------------------------------------------

        /**
         * Generates a random rotation angle (and axis in 3D) and builds a rotation matrix.
         * Rotates both mObject and mPoint to raycast in another random direction in the event the node is unsure.
         * 
         */
        void
        random_rotation();

        //-------------------------------------------------------------------------------

        /**
         * Resets mPoint and mObject to the original coordinate frame in the event that they were rotated
         * 
         */
        void
        undo_rotation();

        //-------------------------------------------------------------------------------
        // Accessor functions
        //-------------------------------------------------------------------------------

        void
        set_point( const Matrix< DDRMat >& aPoint )
        {
            mPoint = aPoint;
        }
        
        moris::Cell< uint >
        get_candidate_facets()
        {
            return mCandidateFacets;
        }

        moris::Cell< Facet* >
        get_intersected_facets()
        {
            return mIntersectedFacets;
        }

        moris::Cell< real >
        get_intersection_coordinates()
        {
            return mCoordsK;
        }

        uint
        is_point_inside()
        {
            return mPointIsInside;
        }
    };    // class Raycast
}    // namespace moris::sdf