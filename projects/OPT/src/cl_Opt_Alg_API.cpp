// Project header files
#include "ios.hpp"
#include "cl_Opt_Alg_API.hpp" // OPT/src
#include "cl_Opt_Alg_GCMMA.hpp" // OPT/src
#include "cl_Opt_Alg_LBFGS.hpp" // OPT/src
#include "cl_Opt_Alg_SQP.hpp" // OPT/src
#include "cl_Opt_Alg_Sweep.hpp" // OPT/src

// -----------------------------------------------------------------------------
extern moris::Logger gLogger;

namespace moris
{
    namespace opt
    {
        OptAlgAPI::OptAlgAPI( const std::string aAlgId )
            : mOptAlg(nullptr)
        {
            if ( aAlgId == "GCMMA" )      // Globally convergent method of moving asymptotes
            {
                mOptAlg = new OptAlgGCMMA();
            }
            else if ( aAlgId == "SQP" )   // Sparse sequential quadratic programming
            {
                mOptAlg = new OptAlgSQP();
            }
            else if ( aAlgId == "SWEEP" ) // Performs sweep of desired optimization variables
            {
                mOptAlg = new OptAlgSweep();
            }
            else
            {
                MORIS_LOG_ERROR ( "Algorithm not yet implemented.");
                assert::error( "In cl_Opt_Alg_API.cpp" );
            }
        }

        // ---------------------------------------------------------------------

        OptAlgAPI::~OptAlgAPI()
        {
            delete mOptAlg;
        }
    } // namespace opt
}     // namespace moris
