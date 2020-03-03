// Project header files
#include "ios.hpp"
#include "cl_OPT_Algorithm_API.hpp" // OPT/src
#include "cl_OPT_Algorithm_GCMMA.hpp" // OPT/src
#include "cl_OPT_Algorithm_LBFGS.hpp" // OPT/src
#include "cl_OPT_Algorithm_SQP.hpp" // OPT/src
#include "cl_OPT_Algorithm_Sweep.hpp" // OPT/src

// -----------------------------------------------------------------------------
extern moris::Logger gLogger;

namespace moris
{
    namespace opt
    {
        Algorithm_API::Algorithm_API(const std::string aAlgId )
            : mAlgorithm(nullptr)
        {
            if (aAlgId == "GCMMA")      // Globally convergent method of moving asymptotes
            {
                mAlgorithm = new OptAlgGCMMA();
            }
            else if (aAlgId == "SQP")   // Sparse sequential quadratic programming
            {
                mAlgorithm = new Algorithm_SQP();
            }
            else if (aAlgId == "LBFGS") // Limited memory Broyden-Fletcher-Goldfarb-Shanno algorithm
            {
                // mAlgorithm = new Algorithm_LBFGS();
                MORIS_ERROR(false, "LBFGS currently is not working.");
            }
            else if (aAlgId == "SWEEP") // Performs sweep of desired optimization variables
            {
                mAlgorithm = new Algorithm_Sweep();
            }
            else
            {
                MORIS_LOG_ERROR ( "Algorithm not yet implemented.");
                assert::error( "In cl_Opt_Alg_API.cpp" );
            }
        }

        // ---------------------------------------------------------------------

        Algorithm_API::~Algorithm_API()
        {
            delete mAlgorithm;
        }
    } // namespace opt
}     // namespace moris
