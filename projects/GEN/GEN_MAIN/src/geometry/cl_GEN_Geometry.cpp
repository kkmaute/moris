#include "cl_GEN_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Geometry::Geometry(sint aNumRefinements, sint aRefinementFunctionIndex)
                : mNumRefinements(aNumRefinements),
                  mRefinementFunctionIndex(aRefinementFunctionIndex)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Geometry::sensitivities_available()
        {
            return true;
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

    }
}
