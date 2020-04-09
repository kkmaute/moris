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
        Algorithm_API::Algorithm_API(ParameterList aParameterList)
            : mAlgorithm(nullptr)
        {
            std::string tAlgorithmName = aParameterList.get<std::string>("algorithm");
            if (tAlgorithmName == "gcmma")      // Globally convergent method of moving asymptotes
            {
                mAlgorithm = new OptAlgGCMMA();
            }
            else if (tAlgorithmName == "sqp")   // Sparse sequential quadratic programming
            {
                mAlgorithm = new Algorithm_SQP();
            }
            else if (tAlgorithmName == "lbfgs") // Limited memory Broyden-Fletcher-Goldfarb-Shanno algorithm
            {
                // mAlgorithm = new Algorithm_LBFGS();
                MORIS_ERROR(false, "LBFGS currently is not working.");
            }
            else if (tAlgorithmName == "sweep") // Performs sweep of desired optimization variables
            {
                mAlgorithm = new Algorithm_Sweep();
            }
            else
            {
                MORIS_LOG_ERROR ( "Algorithm not yet implemented.");
                assert::error( "In cl_Opt_Alg_API.cpp" );
            }
            mAlgorithm->mParameterList = aParameterList;
        }

        // ---------------------------------------------------------------------

        Algorithm_API::~Algorithm_API()
        {
            delete mAlgorithm;
        }
    } // namespace opt
}     // namespace moris
