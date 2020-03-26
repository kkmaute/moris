// Project header files
#include "cl_OPT_Algorithm.hpp" // OPT/src

namespace moris
{
    namespace opt
    {
        Algorithm::Algorithm()
        {
        }

        // -------------------------------------------------------------------------------------------------------------

        Algorithm::~Algorithm()
        {
        }

        // -------------------------------------------------------------------------------------------------------------

        void Algorithm::initialize()
        {
            mActive.set_size( mProblem->get_num_constraints(), 1, 0.0 ); // set the size of the vector containing active constraints flag
        }
    }
}
