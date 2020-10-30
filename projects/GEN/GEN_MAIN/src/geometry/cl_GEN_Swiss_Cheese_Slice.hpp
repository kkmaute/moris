#ifndef MORIS_CL_GEN_SWISS_CHEESE_SLICE_HPP
#define MORIS_CL_GEN_SWISS_CHEESE_SLICE_HPP

#include "cl_GEN_Multigeometry.hpp"

namespace moris
{
    namespace ge
    {
        class Swiss_Cheese_Slice : public Multigeometry
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
                 * @param aName Name of this field for identification
                 * @param aNumRefinements The number of refinement steps to use for this geometry
                 * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
                 * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
                 * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
                 * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
                 */
                Swiss_Cheese_Slice(
                        real        aLeftBound,
                        real        aRightBound,
                        real        aBottomBound,
                        real        aTopBound,
                        uint        aNumXHoles,
                        uint        aNumYHoles,
                        real        aXSemidiameter,
                        real        aYSemidiameter,
                        real        aExponent = 2.0,
                        real        aScaling = 1.0,
                        real        aRegularization = 1e-8,
                        real        aShift = 0.0,
                        real        aOffset = 0.0,
                        std::string aName = "",
                        Matrix<DDSMat>  aNumRefinements = {{}},
                        Matrix<DDSMat>  aNumPatterns = {{}},
                        sint        aRefinementFunctionIndex = -1,
                        sint        aBSplineMeshIndex = -1,
                        real        aBSplineLowerBound = -1.0,
                        real        aBSplineUpperBound = 1.0);

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
                 * @param aName Name of this field for identification
                 * @param aNumRefinements The number of refinement steps to use for this geometry
                 * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
                 * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
                 * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
                 * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
                 */
                Swiss_Cheese_Slice(
                        real        aLeftBound,
                        real        aRightBound,
                        real        aBottomBound,
                        real        aTopBound,
                        real        aTargetXSpacing,
                        real        aTargetYSpacing,
                        real        aXSemidiameter,
                        real        aYSemidiameter,
                        real        aExponent,
                        real        aScaling,
                        real        aRegularization,
                        real        aShift,
                        real        aOffset = 0.0,
                        bool        aAllowLessThanTargetSpacing = false,
                        std::string aName = "",
                        Matrix<DDSMat>  aNumRefinements = {{}},
                        Matrix<DDSMat>  aNumPatterns = {{}},
                        sint        aRefinementFunctionIndex = -1,
                        sint        aBSplineMeshIndex = -1,
                        real        aBSplineLowerBound = -1.0,
                        real        aBSplineUpperBound = 1.0);

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
                 * @param aNumRefinements The number of refinement steps to use for this geometry
                 * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
                 * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
                 * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
                 * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
                 */
                void create_holes(
                        real aXOrigin,
                        real aYOrigin,
                        uint aNumXHoles,
                        uint aNumYHoles,
                        real aXDelta,
                        real aYDelta,
                        real aOffset,
                        real aXSemidiameter,
                        real aYSemidiameter,
                        real aExponent,
                        real aScaling,
                        real aRegularization,
                        real aShift,
                        Matrix<DDSMat>  aNumRefinements,
                        Matrix<DDSMat>  aNumPatterns,
                        sint aRefinementFunctionIndex,
                        sint aBSplineMeshIndex,
                        real aBSplineLowerBound,
                        real aBSplineUpperBound);
        };
    }
}

#endif //MORIS_CL_GEN_SWISS_CHEESE_SLICE_HPP
