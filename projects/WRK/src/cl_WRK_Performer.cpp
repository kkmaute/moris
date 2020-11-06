#include "cl_WRK_Performer.hpp"

namespace moris
{
    namespace wrk
    {
        //--------------------------------------------------------------------------------------------------------------

        uint Performer::get_num_refinement_fields()
        {
            return 0; // No fields by default
        }

        //--------------------------------------------------------------------------------------------------------------

        real Performer::get_field_value(
                uint aFieldIndex,
                uint aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return 0.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        sint Performer::get_refinement_function_index(
                uint aFieldIndex,
                uint aRefinementIndex)
        {
            return -1; // Default refinement
        }

        //--------------------------------------------------------------------------------------------------------------
    }
}
