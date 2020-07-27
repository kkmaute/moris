#ifndef MORIS_CL_GEN_GEOMETRY_HPP
#define MORIS_CL_GEN_GEOMETRY_HPP

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry : virtual public Field
        {
        private:
            sint mNumRefinements;
            sint mRefinementFunctionIndex;
            sint mBSplineMeshIndex;
            real mLevelSetLowerBound;
            real mLevelSetUpperBound;

        public:

            /**
             * Constructor for a geometry, which sets the refinement info and a B-spline mesh index.
             *
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for level set discretization (-1 = no B-splines)
             * @param aLevelSetLowerBound The lower bound for a B-spline level set field describing this geometry
             * @param aLevelSetUpperBound The upper bound for a B-spline level set field describing this geometry
             */
            Geometry(sint aNumRefinements,
                     sint aRefinementFunctionIndex,
                     sint aBSplineMeshIndex,
                     real aLevelSetLowerBound,
                     real aLevelSetUpperBound);

            /**
             * This function will return true when called less than the number of refinements set for this geometry,
             * and false otherwise. This is to determine for a given refinement call if this geometry needs refinement.
             *
             * @return if to perform an additional refinement with this geometry
             */
            sint get_num_refinements();

            /**
             * Gets the index of a user-defined refinement function used within HMR.
             *
             * @return User-defined refinement function index
             */
            sint get_refinement_function_index();

            /**
             * Gets the index of a B-Spline mesh for creating a B-Spline level set discretization for this geometry.
             *
             * @return B-Spline mesh index, or -1 if not creating level set
             */
            sint get_bspline_mesh_index();

            /**
             * Gets the lower bound for the B-spline level set field.
             *
             * @return Lower bound
             */
            real get_level_set_lower_bound();

            /**
             * Get the upper bound for the B-spline level set field.
             *
             * @return Upper bound
             */
            real get_level_set_upper_bound();

            /**
             * Function for determining if this geometry is to be used for seeding a B-spline level set field.
             *
             * @return Logic for level set creation
             */
            virtual bool conversion_to_level_set();

        };
    }
}

#endif /*MORIS_CL_GEN_GEOMETRY_HPP*/
