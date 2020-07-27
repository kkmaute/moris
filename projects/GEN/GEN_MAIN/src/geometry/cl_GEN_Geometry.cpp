#include "cl_GEN_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Geometry::Geometry(sint aNumRefinements,
                           sint aRefinementFunctionIndex,
                           sint aBSplineMeshIndex,
                           real aLevelSetLowerBound,
                           real aLevelSetUpperBound)
                : mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex),
                  mBSplineMeshIndex(aBSplineMeshIndex),
                  mLevelSetLowerBound(aLevelSetLowerBound),
                  mLevelSetUpperBound(aLevelSetUpperBound)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Geometry::get_num_refinements()
        {
            return mNumRefinements;
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Geometry::get_refinement_function_index()
        {
            return mRefinementFunctionIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Geometry::get_bspline_mesh_index()
        {
            return mBSplineMeshIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Geometry::get_level_set_lower_bound()
        {
            return mLevelSetLowerBound;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Geometry::get_level_set_upper_bound()
        {
            return mLevelSetUpperBound;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Geometry::conversion_to_level_set()
        {
            return (mBSplineMeshIndex >= 0);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
