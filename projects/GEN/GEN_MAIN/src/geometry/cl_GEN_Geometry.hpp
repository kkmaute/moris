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

        public:

            /**
             * Constructor for a geometry, which sets the refinement info
             *
             * @param aNumRefinements The number of refinement steps to use for this geometry
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             */
            Geometry(sint aNumRefinements, sint aRefinementFunctionIndex);

            /**
             * Lets the geometry engine know if sensitivities are available, otherwise it will perform finite
             * differencing instead for intersection locations
             *
             * @return If sensitivities are implemented or not (true in base class)
             */
            virtual bool sensitivities_available();

            /**
             * This function will return true when called less than the number of refinements set for this geometry,
             * and false otherwise. This is to determine for a given refinement call if this geometry needs refinement.
             *
             * @return if to perform an additional refinement with this geometry
             */
            sint get_num_refinements();

            /**
             * Gets the index of a user-defined refinement function used within HMR
             *
             * @return User-defined refinement function index
             */
            sint get_refinement_function_index();

        };
    }
}

#endif /*MORIS_CL_GEN_GEOMETRY_HPP*/
