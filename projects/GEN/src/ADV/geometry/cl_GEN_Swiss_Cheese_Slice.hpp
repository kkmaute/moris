/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Swiss_Cheese_Slice.hpp
 *
 */

#pragma once

#include "cl_GEN_Combined_Field.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"

namespace moris::ge
{
    class Swiss_Cheese_Slice : public Combined_Field
    {
    public:

        /**
         * Constructor for a swiss cheese slice, with number of holes specified
         *
         * @param aLeftBound Left (-x) bound on hole center creation
         * @param aRightBound Right (+x) bound on hole center creation
         * @param aBottomBound Bottom (-y) bound on hole center creation
         * @param aTopBound Top (+y) bound on hole center creation
         * @param aNumXHoles Number of holes in the x direction
         * @param aNumYHoles Number of holes in the y direction
         * @param aXSemidiameter Superellipse semi-diameter in the x direction
         * @param aYSemidiameter Superellipse semi-diameter in the y direction
         * @param aExponent Superellipse exponent
         * @param aExponent Superellipse scaling
         * @param aRegularization Superellipse regularization
         * @param aShift Superellipse shift (near zero)
         * @param aOffset Offset to be applied on subsequent rows in the y direction
         * @param aParameters Additional parameters
         */
        Swiss_Cheese_Slice(
                real             aLeftBound,
                real             aRightBound,
                real             aBottomBound,
                real             aTopBound,
                uint             aNumXHoles,
                uint             aNumYHoles,
                real             aXSemidiameter,
                real             aYSemidiameter,
                real             aExponent = 2.0,
                real             aScaling = 1.0,
                real             aRegularization = 1e-8,
                real             aShift = 0.0,
                real             aOffset = 0.0,
              Level_Set_Parameters aParameters = Level_Set_Parameters());

        /**
         * Constructor for a swiss cheese slice, with hole spacing specified
         *
         * @param aLeftBound Left (-x) bound on hole center creation
         * @param aRightBound Right (+x) bound on hole center creation
         * @param aBottomBound Bottom (-y) bound on hole center creation
         * @param aTopBound Top (+y) bound on hole center creation
         * @param aTargetXSpacing Targeted spacing between hole centers in the x direction
         * @param aTargetYSpacing Targeted spacing between hole centers in the y direction
         * @param aXSemidiameter Superellipse semi-diameter in the x direction
         * @param aYSemidiameter Superellipse semi-diameter in the y direction
         * @param aExponent Superellipse exponent
         * @param scaling  Superellipse scaling
         * @param aRegularization Superellipse regularization
         * @param aShift Superellipse shift (near zero)
         * @param aOffset Offset to be applied on subsequent rows in the y direction
         * @param aParameters Additional parameters
         */
        Swiss_Cheese_Slice(
                real             aLeftBound,
                real             aRightBound,
                real             aBottomBound,
                real             aTopBound,
                real             aTargetXSpacing,
                real             aTargetYSpacing,
                real             aXSemidiameter,
                real             aYSemidiameter,
                real             aExponent,
                real             aScaling,
                real             aRegularization,
                real             aShift,
                real             aOffset = 0.0,
                bool             aAllowLessThanTargetSpacing = false,
                Level_Set_Parameters aParameters = Level_Set_Parameters());

    private:

        /**
         * Standard private function for creating holes to eliminate redundant code.
         *
         * @param aXOrigin x-coordinate of the lower left hole center
         * @param aYOrigin y-coordinate of the lower left hole center
         * @param aNumXHoles Number of holes in the x direction
         * @param aNumYHoles Number of holes in the y direction
         * @param aXDelta Hole center delta in the x direction
         * @param aYDelta Hole center delta in the y direction
         * @param aOffset Offset to be applied on subsequent rows in the y direction
         * @param aXSemidiameter Superellipse semi-diameter in the x direction
         * @param aYSemidiameter Superellipse semi-diameter in the y direction
         * @param aExponent Superellipse exponent
         * @param aExponent Superellipse scaling
         * @param aRegularization Superellipse regularization
         * @param aShift Superellipse shift (near zero)
         * @param aNumRefinements The number of refinement steps to use for this field
         * @param aRefinementMeshIndices Indices of meshes to perform refinement on
         * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = {} refinement)
         * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
         * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
         * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
         */
        void create_holes(
                real           aXOrigin,
                real           aYOrigin,
                uint           aNumXHoles,
                uint           aNumYHoles,
                real           aXDelta,
                real           aYDelta,
                real           aOffset,
                real           aXSemidiameter,
                real           aYSemidiameter,
                real           aExponent,
                real           aScaling,
                real           aRegularization,
                real           aShift,
              Level_Set_Parameters aParameters);
    };
}
